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
	
	**Weights, Loss and validation loss saved as files
"""

import os
os.environ['KERAS_BACKEND'] = 'tensorflow'
import keras as keras
import argparse 
import numpy as np
import csv
import pandas as pd

from keras import optimizers, regularizers, layers
from keras.layers import Input, Dense
from keras.models import Model, Sequential
#from tensorflow.keras import initializers
import TiedWeightsEncoder as tw
import Adage as ad


from keras import initializers


def run_model(input_file, seed=123, epochs=50, kl1=0, kl2=0,
			  act='sigmoid', init='glorot_uniform', tied=True, batch_size=10):
	"""

	"""
	print(keras.backend.backend())
	#all_comp = np.loadtxt(open(input_file, "rb"), delimiter=',', skiprows = 1)
	all_comp = pd.read_csv(input_file, index_col=0)
	 # this is the size of our input
	gene_num = np.size(all_comp, 0)
	encoding_dim = 300

	x_train, x_train_noisy = prep_data(all_comp, seed)

	if(tied):
		autoencoder = linked_ae(encoding_dim, gene_num, act, init,
								seed, kl1, kl2)
	else:
		autoencoder = unlinked_ae(encoding_dim, gene_num, act, init,
								  seed, kl1, kl2)
	autoencoder, history = train_model(autoencoder, x_train, 
		                               x_train_noisy, epochs, seed, batch_size)
	weights = autoencoder.get_weights()

	file_desc = (input_file[:-4] + '_seed:' + str(seed)	
				 + "_kl1:" + str(kl1)	
				 + "_kl2:" + str(kl2) 	
				 +  "_act:" + act 		
				 + '_init:' + init 		
				 + '_ep:' + str(epochs) 	
				 + '_tied:' + str(tied)	
				 + '_batch:' + str(batch_size))
				

	write_data(file_desc, weights, history)

	adage = ad.Adage(autoencoder, history, all_comp)
	return adage


def unlinked_ae(encoding_dim, gene_num, act, init,seed, kl1, kl2):
	input_img = Input(shape=(gene_num,))
	init_s = keras.initializers.glorot_normal(seed=seed)
	encoded = Dense(encoding_dim, input_shape=(gene_num, ), 
					activation=act, #sigmoid
					#kernel_initializer = init_s,
					kernel_initializer = init,
					kernel_regularizer = regularizers.l1_l2(l1=kl1, l2=kl2),
    				activity_regularizer = regularizers.l1(0))(input_img) 

	decoded = Dense(gene_num, activation='sigmoid',
    				activity_regularizer = regularizers.l1(0))(encoded)

	# this model maps an input to its reconstruction
	autoencoder = Model(input_img, decoded)
	return(autoencoder)


def linked_ae(encoding_dim, gene_num, act, init,seed, kl1, kl2):
	input_img = Input(shape=(gene_num,))
	init_s = keras.initializers.glorot_normal(seed=seed)
	encoded = Dense(encoding_dim, input_shape=(gene_num, ), activation=act, 
					#kernel_initializer = init_s,
					kernel_initializer = init,
					kernel_regularizer = regularizers.l1_l2(l1=kl1, l2=kl2),
    				activity_regularizer = regularizers.l1(0))


	decoder = tw.TiedWeightsEncoder(input_shape=(encoding_dim,), 
					output_dim=gene_num,encoded=encoded, activation=None)

	autoencoder = Sequential()
	autoencoder.add(encoded)
	autoencoder.add(decoder)
	return(autoencoder)


def prep_data(all_comp, seed):
	all_comp_np = all_comp.to_numpy()
	# this is the size of our input
	gene_num = np.size(all_comp_np, 0)

	x_train = all_comp_np.transpose()
	x_train = x_train.reshape((len(x_train), 
							   np.prod(x_train.shape[1:])))
	noise_factor = 0.1
	x_train_noisy = x_train + (noise_factor 
		            * np.random.normal(loc=0.0,scale=1.0, size=x_train.shape)) 
	x_train_noisy = np.clip(x_train_noisy, 0., 1.)

	return(x_train, x_train_noisy)

def train_model(autoencoder, x_train, x_train_noisy, epochs, seed, batch_size):

	np.random.seed(seed)
	train_idxs = np.random.choice(x_train.shape[0], 
							      int(x_train.shape[0]*0.9), replace=False)

	x_train_train = x_train[train_idxs,:]
	x_train_test = x_train[~np.in1d(range(x_train.shape[0]),train_idxs),:]

	x_train_noise_train = x_train_noisy[train_idxs,:]
	x_train_noise_test = x_train_noisy[~np.in1d(range(x_train.shape[0]),
		                                              train_idxs),:]


	optim = optimizers.Adadelta() # lr=0.001, rho=0.95, epsilon=1e-07
	optim = optimizers.SGD() # lr=0.001, rho=0.95, epsilon=1e-07
	autoencoder.compile(optimizer=optim, loss='binary_crossentropy') # binary_crossentropy mse

	history = autoencoder.fit(x_train_noise_train, x_train_train,
	              	epochs=epochs,
	                batch_size=batch_size,
	                shuffle=True,
	                #validation_split = 0.1,          
	                verbose=1,
	                validation_data=(np.array(x_train_noise_test),
	                				 np.array(x_train_test)))
	return(autoencoder, history)


def write_data(file_desc, weights, history):
	"""
	Save logs and output for a model in an outputs foolder
	"""
	np.savetxt('./outputs/' + file_desc + '_en_weights.csv', 
		np.matrix(weights[0]), fmt = '%s', delimiter=',')
	np.savetxt('./outputs/' + file_desc + '_en_bias', 
		np.matrix(weights[1]), fmt = '%s', delimiter=',')
	np.savetxt('./outputs/' + file_desc + '_de_weights.csv', 
		np.matrix(weights[2]), fmt = '%s', delimiter=',')
	np.savetxt('./outputs/' + file_desc + '_de_bias.csv', 
		np.matrix(weights[3]), fmt = '%s', delimiter=',')
	np.savetxt('./outputs/' + file_desc + '_loss.csv', 
		np.matrix(history.history['loss']), fmt = '%s', delimiter=',')
	np.savetxt('./outputs/' + file_desc + '_val_loss.csv', 
		np.matrix(history.history['val_loss']), fmt = '%s', delimiter=',')


if __name__ == '__main__':
		parser = argparse.ArgumentParser(description='Get training set.')
		parser.add_argument('filename',type=str, nargs=1, 
			help='filpath to training set.')
		args=parser.parse_args()
		run_model(args.filename[0])


