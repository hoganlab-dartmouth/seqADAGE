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

#import os
#os.environ['KERAS_BACKEND'] = 'tensorflow'
#import keras as keras

import tensorflow as tf
from tensorflow import keras

import argparse
import numpy as np
import csv
import pandas as pd
import math

from tensorflow.keras import optimizers, regularizers, layers, initializers, models
from tensorflow.keras.callbacks import TensorBoard, ModelCheckpoint, EarlyStopping, ReduceLROnPlateau
#from keras.layers import Input, Dense
#from keras.models import Model, Sequential
#from tensorflow.keras import initializers
import TiedWeightsEncoder as tw
import Adage as ad
import tensorflow as tf



from keras import initializers


def run_count_autoencoder(input_file, seed=123, epochs=50, kl1=0, kl2=0, lr = 0.01,
			  act='sigmoid', init='glorot_uniform', tied=True, batch_size=10, v=1):
	"""

	"""
	print(keras.backend.backend())
	#all_comp = np.loadtxt(open(input_file, "rb"), delimiter=',', skiprows = 1)
	all_comp = pd.read_csv(input_file, index_col=0)
	#print(all_comp.dtypes)
	 # this is the size of our input
	gene_num = np.size(all_comp, 0)
	encoding_dim = 300

	x_train, x_train_noisy = prep_data(all_comp, seed)
	kl1_tf = tf.cast(kl1, dtype=tf.float64)
	if(tied):
		autoencoder = linked_ae(encoding_dim, gene_num, act, init,
								seed, kl1, kl2)
	else:
		autoencoder = unlinked_ae(encoding_dim, gene_num, act, init,
								  seed, kl1_tf, kl2)
	autoencoder, history = train_model(autoencoder, x_train,
		                               x_train_noisy, epochs, seed, batch_size,lr,v)
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
	input_img = layers.Input(shape=(gene_num,))
	init_s = initializers.glorot_normal(seed=seed)
	encoded = layers.Dense(encoding_dim, input_shape=(gene_num, ),
					activation=act, #sigmoid
					#kernel_initializer = init_s,
					kernel_initializer = init,
					kernel_regularizer = regularizers.l1_l2(l1=kl1, l2=kl2),
    				activity_regularizer = regularizers.l1(0))(input_img)

	decoded = layers.Dense(gene_num, activation='sigmoid',
    				activity_regularizer = regularizers.l1(0))(encoded)

	# this model maps an input to its reconstruction
	autoencoder = models.Model(input_img, decoded)
	return(autoencoder)


def linked_ae(encoding_dim, gene_num, act, init,seed, kl1, kl2):
	input_img = layers.Input(shape=(gene_num,))
	init_s = initializers.glorot_normal(seed=seed)
	encoded = layers.Dense(encoding_dim, input_shape=(gene_num, ), activation=act,
					#kernel_initializer = init_s,
					kernel_initializer = init,
					kernel_regularizer = regularizers.l1_l2(l1=kl1, l2=kl2),
    				activity_regularizer = regularizers.l1(0))


	decoder = tw.TiedWeightsEncoder(input_shape=(encoding_dim,),
					output_dim=gene_num,encoded=encoded, activation="sigmoid")

	autoencoder = keras.Sequential()
	autoencoder.add(encoded)
	autoencoder.add(decoder)
	return(autoencoder)


def prep_data(all_comp, seed):

	all_comp_np = all_comp.values.astype("float64")
	print(np.shape(all_comp_np))
	# this is the size of our input
	gene_num = np.size(all_comp_np, 0)

	x_train = all_comp_np.transpose()
	#x_train = x_train.reshape((len(x_train),
	#						   np.prod(x_train.shape[1:])))
	noise_factor = 0.1
	x_train_noisy = x_train + (noise_factor
		            * np.random.normal(loc=0.0,scale=1.0, size=x_train.shape))
	x_train_noisy = np.clip(x_train_noisy, 0., 1.)

	return(x_train, x_train_noisy)


# def nbinom_pmf_tf(x,n,p):
#     coeff = tf.math.lgamma(n + x) - tf.math.lgamma(x + 1) - tf.math.lgamma(n)
#     return tf.cast(tf.exp(coeff + n * tf.math.log(p) + x * tf.math.log(1 - p)),dtype=tf.float32)


# def loss_neg_bin_tf_differentiable(y_pred, y_true):
# 	print(tf.cast(y_pred,dtype=tf.float32).shape)
# 	print(tf.cast(y_true,dtype=tf.float32).shape)
# 	result = tf.map_fn(lambda x: -nbinom_pmf_tf(x[1]
# 												, x[0][0]
#                                                 , tf.minimum(tf.constant(0.99,dtype=tf.float32),x[0][1]))
# 					   ,(tf.cast(y_pred, dtype=tf.float32),tf.cast(y_true,dtype=tf.float32))
#                        ,dtype=tf.float32)
# 	print(result.shape)
# 	result = tf.reduce_sum(result,axis=1)
# 	print(result.shape)
# 	print(result)
# 	return result

def nbinom_pmf_tf(x,n,p):
    coeff = tf.math.lgamma(n + x) - tf.math.lgamma(x + 1) - tf.math.lgamma(n)
    return tf.cast(tf.exp(coeff + n * tf.math.log(p) + x * tf.math.log(1 - p)),dtype=tf.float32)

def nbl(y_true,n,p):
	nll = (tf.math.lgamma(n)
		+ tf.math.lgamma(y_true +1)
		- tf.math.lgamma(n + y_true)
		- n * tf.math.log(p)
		- y_true * tf.math.log(1-p)
		)
	return nll

def nbl2(y_true,y_pred):
	y_true = tf.cast(y_true, tf.float32)
	y_pred = tf.cast(y_pred, tf.float32)
	t1 = tf.math.lgamma(1e6 + 1e-10) + tf.math.lgamma(y_true + 1.0) - tf.math.lgamma(y_true+1e6+1e-10)
	t2 = (1e6+y_true) * tf.math.log(1.0 + (y_pred/(1e6+1e-10))) + (y_true * (tf.math.log(1e6+1e-10) - tf.math.log(y_pred + 1e-10)))
	final = t1+t2
	final = tf.reduce_mean(final)
	return final

def _nan2inf(x):
    return tf.where(tf.math.is_nan(x), tf.zeros_like(x)+np.inf, x)

def sigmoid(x):
	print(x)
	return 1 / (1 + math.exp(-x))

def zinbl2(y_true,y_pred):
	y_true = tf.cast(y_true, tf.float32)
	y_pred = tf.cast(y_pred, tf.float32)
	eps = 1e-10
	theta = 1e6
	#print(pi)

	#pi = y_pred / y_pred + math.exp
	pi = 0.5
	ridge_lambda = 0
	t1 = tf.math.lgamma(1e6 + 1e-10) + tf.math.lgamma(y_true + 1.0) - tf.math.lgamma(y_true+1e6+1e-10)
	t2 = (1e6+y_true) * tf.math.log(1.0 + (y_pred/(1e6+1e-10))) + (y_true * (tf.math.log(1e6+1e-10) - tf.math.log(y_pred + 1e-10)))

	nb_case = t1+t2 - tf.math.log(1.0 - pi+eps)

	zero_nb = tf.pow(theta/(theta+y_pred+eps), theta)
	zero_case = -tf.math.log(pi + (( 1.0-pi)*zero_nb)+eps)
	result = tf.where(tf.less(y_true,1e-8), zero_case, nb_case)
	ridge = ridge_lambda*tf.square(pi)
	result += ridge


	final = tf.reduce_mean(result)
	return _nan2inf(final)

def loss_neg_bin_tf_differentiable2(y_pred, y_true):
    result = tf.map_fn(lambda x: nbl2(x[1] # -nbinom_pmf_tf   nbl
                       , x[0])
                       ,(y_pred,y_true)
                       ,dtype=tf.float32)
    result = tf.reduce_sum(result,axis=1)
    return result

def loss_neg_bin_tf_differentiable(y_pred, y_true):
    result = tf.map_fn(lambda x: nbl2(x[1] # -nbinom_pmf_tf   nbl
                                                , x[0][0]
                                                , tf.minimum(tf.constant(0.99,dtype=tf.float32),x[0][1]))
                       ,(y_pred,y_true)
                       ,dtype=tf.float32)
    result = tf.reduce_sum(result,axis=1)
    return result


def train_model(autoencoder, x_train, x_train_noisy, epochs, seed, batch_size, lr,v):

	np.random.seed(seed)
	train_idxs = np.random.choice(x_train.shape[0],
							      int(x_train.shape[0]*0.9), replace=False)

	x_train_train = x_train[train_idxs,:]
	x_train_test = x_train[~np.in1d(range(x_train.shape[0]),train_idxs),:]

	x_train_noise_train = x_train_noisy[train_idxs,:]
	x_train_noise_test = x_train_noisy[~np.in1d(range(x_train.shape[0]),
		                                              train_idxs),:]
	#print(x_train_noise_test.shape)
	#print(x_train_test.shape)


	optim = optimizers.Adadelta(lr = 0.001) # lr=0.001, rho=0.95, epsilon=1e-07
	optim = optimizers.RMSprop(lr = lr, clipvalue=5.0)
	#autoencoder.compile(optimizer=optim, loss='binary_crossentropy') # binary_crossentropy mse
	autoencoder.compile(optimizer=optim, loss=zinbl2) # binary_crossentropy mse zinbl2

	callbacks = []

	lr_cb = ReduceLROnPlateau(monitor='val_loss', patience = 10, verbose=True)
	callbacks.append(lr_cb)

	history = autoencoder.fit(x_train_noisy, x_train,
	              	epochs=epochs,
	                batch_size=batch_size,
	                shuffle=True,
	                validation_split = 0.1,
	                callbacks=callbacks,
	                verbose=v#,
	                #validation_data=(np.array(x_train_noise_test),
	                #				 np.array(x_train_test))
	                )
	return(autoencoder, history)


def write_data(file_desc, weights, history):
	"""
	Save logs and output for a model in an outputs foolder
	"""
	print(np.shape(weights))
	np.savetxt('../outputs/' + file_desc + '_en_weights_dca.csv',
		np.matrix(weights[0]), fmt = '%s', delimiter=',')
	np.savetxt('../outputs/' + file_desc + '_en_bias_dca',
		np.matrix(weights[1]), fmt = '%s', delimiter=',')
	#np.savetxt('../outputs/' + file_desc + '_de_weights.csv',
	#	np.matrix(weights[2]), fmt = '%s', delimiter=',')
	#np.savetxt('../outputs/' + file_desc + '_de_bias.csv',
#		np.matrix(weights[3]), fmt = '%s', delimiter=',')
	np.savetxt('../outputs/' + file_desc + '_loss_dca.csv',
		np.matrix(history.history['loss']), fmt = '%s', delimiter=',')
	np.savetxt('../outputs/' + file_desc + '_val_loss_dca.csv',
		np.matrix(history.history['val_loss']), fmt = '%s', delimiter=',')


if __name__ == '__main__':
		parser = argparse.ArgumentParser(description='Get training set.')
		parser.add_argument('filename',type=str, nargs=1,
			help='filpath to training set.')
		args=parser.parse_args()
		run_count_autoencoder(args.filename[0])
