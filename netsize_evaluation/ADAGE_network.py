'''
This script provides a helper function to read in ADAGE network from a network
file. Use it with import.
'''

import numpy


def read_network_file(network_file, input_size, net_size):
    '''
    A function to read in network file. The network file stores weight matrix,
    hidden bias vector, and visible bias vector.

    input:
        network_file    file path to the network file
        input_size      the number of features in the input layer
        net_size        the number of nodes in the hidden layer
    output:
        A tuple (weight_matrix, hidden bias, visible bias), each in the numpy
        array format
    '''
    with open(network_file, 'r') as network_fh:
        network_fh.next()  # skip the layer count line
        network_fh.next()  # skip 'weight matrix' line

        # Get the weight matrix
        W = []
        input_count = 0
        for line in network_fh:
            line = line.strip().split('\t')
            W.append(line)
            input_count += 1
            # when it reaches the input size
            if input_count == input_size:
                break
        W = numpy.array(W, dtype=float)

        network_fh.next()  # skip 'hidden bias vector' line

        # Get the bias vector of hidden layer
        h_bias = []
        output_count = 0
        for line in network_fh:
            line = line.strip()
            h_bias.append(line)
            output_count += 1
            # when it reaches the number of nodes in hidden layer
            if output_count == net_size:
                break
        h_bias = numpy.array(h_bias, dtype=float)

        network_fh.next()  # skip 'visible bias vector' line

        # Get the bias vector of visible output layer
        v_bias = []
        input_count = 0
        for line in network_fh:
            line = line.strip()
            v_bias.append(line)
            input_count += 1
            if input_count == input_size:
                break
        v_bias = numpy.array(v_bias, dtype=float)

        parameters = (W, h_bias, v_bias)
        return parameters
