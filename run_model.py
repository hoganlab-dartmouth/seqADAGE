"""

Georgia Doing 2020
Denoising Autoencoder with P.a. RNAseq

Usage:

	run_model.py <dataset> [--seed=<SEED>]

Options:
	-h --help			Show this screen
	<dataset>			File patht to RNAseq compnedium
	--seed = <SEED>		Random seed for training [default: 1]

Output:
	
	Loss and validation loss
"""


import argparse 
import numpy as np
import csv
import pandas as pd

from keras.layers import Input, Dense
from keras.models import Model
from keras import optimizers
from keras import regularizers
from tensorflow.keras import initializers
from keras import layers
import keras as keras





def run_model(input_file,seed):

	all_comp = np.loadtxt(open(input_file, "rb"), delimiter=',', skiprows = 1)
	gene_num = np.size(all_comp, 0)

	x_train = all_comp.transpose()
	x_train = x_train.reshape((len(x_train), np.prod(x_train.shape[1:])))
	noise_factor = 0.1
	x_train_noisy = x_train + noise_factor * np.random.normal(loc=0.0, scale=1.0, size=x_train.shape) 
	x_train_noisy = np.clip(x_train_noisy, 0., 1.)

	# this is the size of our encoded representations
	encoding_dim = 50  
	print(encoding_dim)

	# this is our input placeholder
	input_img = Input(shape=(gene_num,))

	# "encoded" is the encoded representation of the input
	init_s = keras.initializers.random_normal(mean=0.0, stddev=0.05, seed=seed)
	encoded = Dense(encoding_dim, activation='sigmoid',
					kernel_initializer = init_s,
					bias_initializer='zeros')(input_img)
    #activity_regularizer = regularizers.l1(10e-5)

	# "decoded" is the lossy reconstruction of the input
	decoded = Dense(gene_num, activation='sigmoid',
			kernel_initializer = init_s,
			bias_initializer='zeros')(encoded)

	# this model maps an input to its reconstruction
	autoencoder = Model(input_img, decoded)

	# this model maps an input to its encoded representation
	encoder = Model(input_img, encoded)

	# create a placeholder for an encoded (32-dimensional) input
	encoded_input = Input(shape=(encoding_dim,))
	# retrieve the last layer of the autoencoder model
	decoder_layer = autoencoder.layers[-1]
	# create the decoder model
	decoder = Model(encoded_input, decoder_layer(encoded_input))

	autoencoder.compile(optimizer="adadelta", loss='binary_crossentropy') 
	history = autoencoder.fit(x_train_noisy, x_train,
	                epochs=5,
	                batch_size=10,
	                shuffle=True,
	                validation_split = 0.1,          
	                verbose=1)

	weights = autoencoder.get_weights()

	np.savetxt('./outputs/' + input_file[:-4]+  '_en_weights.csv', np.matrix(weights[0]), fmt = '%s', delimiter=',')
	np.savetxt('./outputs/' + input_file[:-4]+ '_en_bias', np.matrix(weights[1]), fmt = '%s', delimiter=',')
	np.savetxt('./outputs/' + input_file[:-4]+ '_de_weights.csv', np.matrix(weights[2]), fmt = '%s', delimiter=',')
	np.savetxt('./outputs/' + input_file[:-4]+ '_de_bias.csv', np.matrix(weights[3]), fmt = '%s', delimiter=',')
	np.savetxt('./outputs/' + input_file[:-4]+ '_loss.csv', np.matrix(history.history['loss']), fmt = '%s', delimiter=',')
	np.savetxt('./outputs/' + input_file[:-4]+ '_val_loss.csv', np.matrix(history.history['val_loss']), fmt = '%s', delimiter=',')


if __name__ == '__main__':
		parser = argparse.ArgumentParser(description='Get training set.')
		parser.add_argument('filename',type=str, nargs=1, help='filpath to training set.')

		args=parser.parse_args()
		print(args.filename[0])

		run_model(args.filename[0])


