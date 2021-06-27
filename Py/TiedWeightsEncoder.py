"""
Georgia Doing 2020
TiedWeightsEncoder class

A class for a tied-weight autoencoder.

Usage
	constructor

"""
#import os
##os.environ['KERAS_BACKEND'] = 'theano'
#from keras.layers import Input, Dense, Activation, Layer
#from keras import activations
#import keras.backend as K


import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import optimizers, regularizers, layers, initializers, models, activations
from tensorflow.keras.layers import Layer


class TiedWeightsEncoder(Layer):

	"""
	Encoder - Decoder tied weights
	"""

	def __init__(self, output_dim, encoded, activation="sigmoid", **kwargs):
		self.output_dim = output_dim
		self.encoder = encoded
		self.activation = activations.get(activation)
		super(TiedWeightsEncoder, self).__init__(**kwargs)

	def build(self, input_shape):
		self.kernel = self.encoder.weights
		super(TiedWeightsEncoder, self).build(input_shape)

	def call(self, x):
		output = tf.matmul(x - self.encoder.weights[1], tf.transpose(self.encoder.weights[0]), name = 'dot')
		#print(output)
		if self.activation is not None:
			output = self.activation(output)
		return output

	def compute_output_shape(self, input_shape):
		return input_shape[0], self.output_dim
