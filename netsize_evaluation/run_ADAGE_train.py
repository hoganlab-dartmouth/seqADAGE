'''
This script is used to run the ADAGE_train.py script to build a batch of
ADAGE models with specific training parameters.

NOTE: This step is very time-consuming if each model is built sequentially.
This script just provides an example of how to build many ADAGE models.
We highly recommend configuring this script so that the construction
of each ADAGE model is done as one job and distributed across a
computing cluster.

Usage:
    run_ADAGE_train.py <input_file> <output_path> <model_sizes>
    <corruption_levels> <initialization_seeds> [--reshuffle=<shuffle_seeds>]
    [--subsample=<subsample_sizes>] [--fileCheck=<fileCheck>]
    run_ADAGE_train.py -h | --help

Options:
    -h --help                       Show this screen.
    <input_file>                    file path to the input data
    <output_path>                   file path to the output folder
    <model_sizes>                   model sizes, numbers separated by comma
    <corruption_levels>             corruption levels, numbers separated by
                                    comma
    <initialization_seeds>          range of initialization random seeds,
                                    two ints separated by comma
    --reshuffle=<shuffle_seeds>     range of reshuffling random seeds
                                    two ints separated by comma [default: 0,0]
    --subsample=<subsample_sizes>   sample sizes to subsample from the input
                                    data, numbers separated by comma
                                    [default: None]
    --fileCheck=<fileCheck>         check whether all models already exist in
                                    the designated output folder and
                                    only build non-existing models
                                    [default: False]
'''

import os
import subprocess
from docopt import docopt


def file_check(file_to_check, run_command):
    '''
    This function checks whether a file already exist. If not, it will
    carry out commands in run_command.

    Inputs:
    file_to_check: the file to check existence, should be full file path
    run_command: commands to run in the terminal if the file does not exist
    '''

    if not os.path.exists(file_to_check):
        print(file_to_check + " does not exist, resubmit it.")
        subprocess.call(run_command, shell=True)
    else:
        print(file_to_check + " already exists.")


def run_ADAGE_train(input_file, output_path, model_sizes, corruption_levels,
                    initialization_seeds, shuffle_seeds,
                    subsample_sizes=None, fileCheck=False):
    '''
    This function runs a batch of ADAGE_train.py script in the terminal
    with the specified run parameters.

    Inputs:
    input_file: character, file path to the input data file for ADAGE_train.py
    output_path: character, the output folder to store ADAGE_train.py outputs
    model_sizes: character, a series of model sizes to use, numbers separated
                 by comma
    corruption_levels: character, a series of corruption levels to use, numbers
                       separated by comma
    initialization_seeds: character, range of initialization random seeds, two
                          numbers separated by comma
    shuffle_seeds: character, range of reshuffling random seeds, two numbers
                   separated by comma
    subsample_sizes: a series of sample size to use in the subsampling analysis,
                     numbers separated by comma
    fileCheck: boolean, whether to check the model existence before build

    '''

    model_sizes = model_sizes.strip().split(",")
    corruption_levels = corruption_levels.strip().split(",")
    init_start, init_end = map(int, initialization_seeds.strip().split(","))
    initialization_seeds = map(str, range(init_start, init_end))
    shuffle_start, shuffle_end = map(int, shuffle_seeds.strip().split(","))
    shuffle_seeds = map(str, range(shuffle_start, shuffle_end))
    if subsample_sizes is not None:
        subsample_sizes = subsample_sizes.strip().split(",")

    # following parameters are default to a set of values suitable for training
    # gene expression data examined previously
    skip_col = '0'
    batch_size = '10'
    epoch_size = '500'
    learn_rate = '0.01'

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    for o in model_sizes:
        if not os.path.exists(os.path.join(output_path, o)):
            os.mkdir(os.path.join(output_path, o))
        new_output_path = os.path.join(output_path, o)
        for c in corruption_levels:
            for l in initialization_seeds:
                run_args = [input_file, skip_col, o, batch_size, epoch_size,
                            c, learn_rate, new_output_path, l]
                name_args = [os.path.basename(input_file.replace('.pcl', '')),
                             o, batch_size, epoch_size, c, learn_rate, l]

                # to keep it simple, we will use the same random seed both as
                # initialization seed and reshuffling seed in most cases,
                # and only differentiate them in the subsampling analysis when
                # shuffle_seeds does not equal to 0.
                if len(shuffle_seeds) == 0:

                    # modify this run_command to submit it as a job in
                    # a computing cluster
                    run_command = ("python ./ADAGE_train.py {a[0]} {a[1]} "
                                   "{a[2]} {a[3]} {a[4]} {a[5]} {a[6]} {a[7]}"
                                   " --seed1 {a[8]} --seed2 {a[8]}"
                                   ).format(a=run_args)

                    # if using fileCheck, it will first check whether the model
                    # already exists, and will only build it when not. If
                    # fileCheck is not used, it will build the model and the old
                    # model will be overwritten.
                    if fileCheck:
                        output_filename = ("{a[0]}_{a[1]}_batch{a[2]}_"
                                           "epoch{a[3]}_corrupt{a[4]}_lr{a[5]}"
                                           "_seed1_{a[6]}_seed2_{a[6]}_network_"
                                           "SdA.txt").format(a=name_args)
                        output_file = os.path.join(new_output_path,
                                                   output_filename)
                        file_check(output_file, run_command)
                    else:
                        subprocess.call(run_command, shell=True)

                else:
                    for p in shuffle_seeds:
                        for s in subsample_sizes:

                            # modify this run_command to submit it as a job in
                            # a computing cluster
                            run_args = run_args.extend([p, s])
                            run_command = ("python ./ADAGE_train.py {a[0]} "
                                           "{a[1]} {a[2]} {a[3]} {a[4]} {a[5]} "
                                           "{a[6]} {a[7]} --seed1 {a[9]} "
                                           "--seed2 {a[8]} --subsample "
                                           "{a[10]}").format(a=run_args)

                            if fileCheck:
                                name_args = name_args.extend([p, s])
                                output_filename = ("{a[0]}_{a[1]}_batch{a[2]}_"
                                                   "epoch{a[3]}_corrupt{a[4]}_"
                                                   "lr{a[5]}_seed1_{a[7]}_seed2"
                                                   "_{a[6]}_subsize_{a[8]}_"
                                                   "network_SdA.txt"
                                                   ).format(a=name_args)
                                output_file = os.path.join(new_output_path,
                                                           output_filename)
                                file_check(output_file, run_command)
                            else:
                                subprocess.call(run_command, shell=True)


if __name__ == '__main__':

    arguments = docopt(__doc__, version=None)
    input_file = arguments['<input_file>']
    output_path = arguments['<output_path>']
    model_sizes = arguments['<model_sizes>']
    corruption_levels = arguments['<corruption_levels>']
    initialization_seeds = arguments['<initialization_seeds>']
    shuffle_seeds = arguments['--reshuffle']
    subsample_sizes = arguments['--subsample']
    fileCheck = arguments['--fileCheck']

    # convert fileCheck string to boolean
    if fileCheck == "False" or fileCheck == "F" or fileCheck == "FALSE":
        fileCheck = False
    elif fileCheck == "True" or fileCheck == "T" or fileCheck == "TRUE":
        fileCheck = True
    else:
        raise ValueError("Cannot covert {} to a boolean".format(fileCheck))

    run_ADAGE_train(input_file, output_path, model_sizes, corruption_levels,
                    initialization_seeds, shuffle_seeds, subsample_sizes,
                    fileCheck)
