'''
This script calculates reconstruction costs of all models in a target
folder. If the test_size is set to 0 (default), it will calculate the mean
reconstruction cost of all samples in the data compendium on a specific model.
If the test_size does not equal to 0, the script considers the provided models
being built for subsampling analysis, and will try to determine the train size
from a model's file name and use the user provided test_size to calculate train
and test reconstruction costs. The test_size should not exceed the total sample
size minus the max sample size used in the subsampling analysis.


Usage:
    run_ADAGE_cost.py <parent_directory> <data_file> <out_file>
    <sample_size> [--testsize=<test_size>]
    run_ADAGE_cost.py -h | --help

Options:
    -h --help               Show this screen.
    <parent_directory>      the parent directory that stores all models
                            with different network sizes
    <data_file>             file path to the data compendium file
    <out_file>              file path to the output file, a tab-delimited file
                            that stores the train and test cost for each model
    <sample_size>           the number of samples in the compendium, int
    --testsize=<test_size>  the number of samples in the test set, int
                            [default: 0].

'''

import os
import csv
from docopt import docopt
import ADAGE_cost


if __name__ == '__main__':

    arguments = docopt(__doc__, version=None)
    parent_directory = arguments['<parent_directory>']
    data_file = arguments['<data_file>']
    out_file = arguments['<out_file>']
    sample_size = arguments['<sample_size>']
    test_size = arguments["--testsize"]

    net_structure = [name for name in os.listdir(parent_directory) if
                     os.path.isdir(os.path.join(parent_directory, name))]

    with open(out_file, 'wb') as out_fh:
        writer = csv.writer(out_fh, delimiter='\t')
        if test_size != "0":
            writer.writerow(["NetworkStructure", "TrainSize", "TestSize",
                             "CorruptLevel", "Seed1", "Seed2", "TrainCost",
                             "TestCost"])
        else:
            writer.writerow(["NetworkStructure", "CorruptLevel", "Seed1",
                             "Seed2", "TrainCost"])
        for i in net_structure:
            directory = os.path.join(parent_directory, i)
            for root, dirs, filenames in os.walk(directory):
                for f in filenames:
                    if 'network' in f:
                        net_file = os.path.join(root, f)
                        print "processing " + net_file
                        name_list = f.split('_')
                        netsize = i
                        corruption = name_list[4]
                        seed1 = name_list[7]
                        seed2 = name_list[9]

                        if test_size != "0":
                            # for models built for subsampling analysis, its
                            # subsampling size is specified in the file name
                            train_size = name_list[11]
                        else:
                            train_size = sample_size

                        costs = ADAGE_cost.get_cost(data_file=data_file,
                                                    network_file=net_file,
                                                    net_size=int(netsize),
                                                    skip_col=0,
                                                    permute_seed=int(seed2),
                                                    train_size=int(train_size),
                                                    test_size=int(test_size))
                        if test_size != "0":
                            writer.writerow([netsize, train_size, test_size,
                                             corruption, seed1, seed2,
                                             str(costs[0]), str(costs[1])])
                        else:
                            writer.writerow([netsize, corruption, seed1, seed2,
                                             str(costs[0])])
