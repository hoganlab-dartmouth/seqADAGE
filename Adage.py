"""
Georgia Doing 2020
ADAGE class

A class for a tied-weight autoencoder.

Usage
	constructor

"""


from keras.models import Model, Sequential
import numpy as np
import pandas
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

class Adage(object):

	"""
	
	"""

	def __init__(self, autoencoder, history, train_comp):
		self.autoencoder = autoencoder
		#super(Adage, self).__init__(**kwargs)
		self.history = history
		self.weights = autoencoder.get_weights()[0]
		self.loss = history.history['loss']
		self.val_loss = history.history['val_loss']
		self.hwg_cutoff = 2.5
		self.hw_genes = np.abs(self.weights) > (np.std(self.weights) * 2.5)
		self.kegg_ps = []
		self.go_ps = []
		self.compendium = train_comp
		self.activities = np.dot( self.compendium.T, self.weights)
		



	def set_hwg_cutoff(self, x):
	 	self.hwg_cutoff = x
	 	self.hw_genes = np.abs(self.weights) > (np.std(self.weights) * x)
	 	return self.hw_genes

	def calc_enrich(self, path_file):
		kegg = pandas.read_csv(path_file, 
                                  header = None, sep = '\t')
		temp_kegg_en = [1]*self.weights.shape[1]
		for i in range(self.weights.shape[1]):
			path_en = []
			sig_genes = self.compendium.index[np.where(self.hw_genes[:,i])]
			for j in range(kegg.shape[0]):
				path_genes = kegg[2][j].split(';')
				x = len(list(set(sig_genes) & set(path_genes)))
				n = self.weights.shape[0]
				k = len(sig_genes)
				m = len(path_genes)
				p = hypergeom.logsf(x,n,k,m)
				path_en.append(-p)
			#path_en_c = multipletests(path_en)[1]
			temp_kegg_en[i] = path_en
		kegg_df = pandas.DataFrame(temp_kegg_en, columns = kegg[0])
		return(kegg_df.replace([np.inf, -np.inf, 'nan'], 0))	

	def set_kegg(self, path_file):
		self.kegg_ps = self.calc_enrich(path_file)
		return(self.kegg_ps)

	def set_go(self, path_file):
		self.go_ps = self.calc_enrich(path_file)
		return(self.go_ps)

	def set_act(self, comp):
	 	self.activities = np.dot(comp, self.weights)
	 	return self.activities

	def calc_act(self, gene_exp):
	 	return np.dot(gene_exp, self.weights)