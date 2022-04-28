"""
This script provides functions to calculate reconstruction errors for an
ADAGE model.
If running in command line, it outputs the train and test reconstruction errors
for an model using certain-sized samples in the dataset as training set and
certain-sized samples in the datasets as testing set. Training and testing sets
are selected randomly.
If you want to get the reconstruction error of a model using all samples in
the dataset, set the trainsize to the number of samples in the dataset and set
testsize to 0.

Usage:
    ADAGE_cost.py <data-file> <skip-col> <network-file> <net-structure>
    [--seed2=<SEED2>] [--trainsize=<train-size>]
    [--testsize=<test-size>]
    ADAGE_cost.py -h | --help

Options:
    -h --help                       Show this screen.
    <data-file>                     File path of the test sample file
    <skip-col>                      int, the number of column to be skipped
                                    between the first gene ID column and the
                                    first experimental column
    <network-file>                  Existing network that are used in the test
    <net-structure>                 The net structure used by the existing
                                    network. A list of int separated by comma
    --seed2=<SEED2>                 Random seed for permuting order of samples
                                    [default: 123]
    --trainsize=<train-size>        The size of the training set [default: 500]
    --testsize=<test-size>          The size of the testing set [default: 200]

"""

from docopt import docopt
import numpy
import ADAGE_network
import sys
sys.path.insert(0, '../data_collection/')  # for importing PCLfile
from pcl import PCLfile


def sigmoid(x):
    '''sigmoid function'''
    y = 1/(1 + numpy.exp(-x))
    return y


def recon_error(x, z):
    '''
    this function calculates the reconstruction error between input layer x
    and reconstruction layer z for one sample.
    '''
    L = -sum(x*numpy.log(z)+(1-x)*numpy.log(1-z))
    return L


def compute_cost(data_set, network_parameters):
    '''
    this function computes the mean reconstruction error when testing a dataset
    on a model specified by its network parameters.

    input:
        data_set    a numpy array that stores one sample per row
        network_parameters  the network parameters got from read_network_file
                            function in ADAGE_network.py. It's a tuple
                            (weight_matrix, hidden bias, visible bias),
                            each in the numpy array format.
    return:
        the averaged reconstruction error of all samples
    '''
    W = numpy.matrix(network_parameters[0])
    W_prime = W.T
    bhid = network_parameters[1]
    bvis = network_parameters[2]
    x = numpy.matrix(data_set)
    y = sigmoid(x * W + bhid)
    z = sigmoid(y * W_prime + bvis)
    L = []
    # calculate the reconstruction error for each sample
    for sample_idx in xrange(x.shape[0]):
        x_row = x[sample_idx, :]
        # convert the matrix format into a flattened array
        x_row = numpy.array(x_row).flatten()
        z_row = z[sample_idx, :]
        z_row = numpy.array(z_row).flatten()
        cost_row = recon_error(x_row, z_row)
        L.append(cost_row)
    cost = numpy.mean(L)
    return cost


def cost_per_node_per_sample(data_set, network_parameters):
    '''
    This function calculates a reconstruction cost for each node-sample
    combination.

    input:
        data_set    a numpy array that stores one sample per row
        network_parameters  the network parameters got from read_network_file
                            function in ADAGE_network.py. It's a tuple
                            (weight_matrix, hidden bias, visible bias),
                            each in the numpy array format.
    return:
        a list of lists, the inside list stores reconstruction error for one
        sample across all nodes, the outside list stores the reconstruction
        error lists across samples

    '''
    W = network_parameters[0]
    bhid = network_parameters[1]
    bvis = network_parameters[2]
    net_size = W.shape[1]
    sample_size = data_set.shape[0]

    cost_matrix = []
    for sample_idx in xrange(sample_size):
        cost_per_sample = []
        for node_idx in xrange(net_size):
            w_vector = W[:, node_idx]
            # reshape the w vector to make it into N*1 matrix form
            w_matrix = numpy.reshape(w_vector, (W.shape[0], 1))
            bhid_value = bhid[node_idx]
            bvis_value = bvis[node_idx]
            network_parameters = (w_matrix, bhid_value, bvis_value)
            x_vector = data_set[sample_idx, :]
            # reshape the x vector to make it into 1*M matrix form
            x_matrix = numpy.reshape(x_vector, (1, data_set.shape[1]))
            cost_per_node_sample = compute_cost(x_matrix,
                                                network_parameters)
            cost_per_sample.append(cost_per_node_sample)
        cost_matrix.append(cost_per_sample)
    return cost_matrix


def get_cost(data_file, network_file, net_size, skip_col=0, permute_seed=123,
             train_size=500, test_size=200):
    '''
    This function calculates the train and test reconstruction errors
    for an model using certain-sized samples in the dataset as training set and
    certain-sized samples in the datasets as testing set. Training and testing
    sets are selected randomly.

    input:
    data_file           File path of the test sample file
    network_file        Existing network that are used in the test
    net_size            int, number of hidden nodes in the model
    skip_col            int, the number of column to be skipped
                        between the first gene ID column and the
                        first experimental column [default: 0]
    permute_seed        int, random seed for permuting order of samples
                        [default: 123]
    train_size          int, the size of the training set [default: 500]
    test_size           int, the size of the testing set [default: 200]

    return:
    A list [train_cost, test_cost] that stores the reconstruction errors
    of the training set and testing set. If train_size equals the compendium
    size, only [train_cost] is returned.

    '''

    # read in the data compendium
    datasets = PCLfile(data_file, skip_col)
    compendium_size = len(datasets.sample_list)

    # get_train_test_sample function first permutes sample order in datasets
    # and then returns first train_size samples as train_set and the next
    # test_size samples as test_set. train_sample and test_sample stores
    # the corresponding sample labels.
    train_set, train_sample, test_set, test_sample = \
        datasets.get_train_test_sample(seed=permute_seed,
                                       train_size=train_size,
                                       test_size=test_size)
    gene_size = train_set.shape[1]  # number of genes

    # read in weight matrix and bias vectors
    network_parameters = ADAGE_network.read_network_file(network_file,
                                                         gene_size,
                                                         net_size)
    # calculate train and test costs
    train_cost = compute_cost(train_set, network_parameters)
    if train_size != compendium_size:
        test_cost = compute_cost(test_set, network_parameters)
        return [train_cost, test_cost]
    else:
        return [train_cost]


if __name__ == '__main__':

    arguments = docopt(__doc__, version=None)
    data_file = arguments['<data-file>']
    skip_col = int(arguments['<skip-col>'])
    network_file = arguments['<network-file>']
    net_size = int(arguments['<net-structure>'])
    random_seed_2 = int(arguments['--seed2'])
    train_size = int(arguments['--trainsize'])
    test_size = int(arguments['--testsize'])

    costs = get_cost(data_file=data_file, network_file=network_file,
                     net_size=net_size, skip_col=skip_col,
                     permute_seed=random_seed_2, train_size=train_size,
                     test_size=test_size)

    if len(costs) == 2:
        print ("Training cost: " + str(costs[0]) + "\tTesting cost: " +
               str(costs[1]))
    else:
        print ("Training cost: " + str(costs[0]))
