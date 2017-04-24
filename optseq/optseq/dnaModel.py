import numpy as np 
from keras.models import Sequential
from keras.layers.core import Dense, Dropout, Activation, Flatten
import random
import pandas as pd
from sklearn.cross_validation import train_test_split
from sklearn import preprocessing
import scipy.stats as stats

np.random.seed(1337)

class dnaModel(object):
	""" A CNN model for DNA sequences """
	def __init__(self, input_file):
		#The input_file should contain only two columns: sequence + output
		#this should create a CNN/model object which can be trained,tested optimally (hyperparam opt)
		
		self.model, self.X_train, self.Y_train, self.X_test, self.Y_test = self.__create_model(input_file)

	def __create_model(self, input_file):
		xtrain, ytrain, xtest, ytest = [], [], [], []

		raw_data = pd.read_excel(input_file, header=0, parse_cols="A,B")
		print raw_data.columns

		#Cleaning expression column
		df = raw_data[np.isfinite(raw_data[u' expression'])]

		cnn = Sequential()
		return cnn, xtrain, ytrain, xtest, ytest

	def __oneHotEncoder(seq):
		base_dict = {u'A':[1,0,0,0],u'C':[0,1,0,0],u'G':[0,0,1,0],u'T':[0,0,0,1]}
		return np.array([base_dict[x] for x in seq])
		
	def __oneHotDecoder(encseq):
		dec_seq = ""
		for x in encseq:
			if (x == np.array([1,0,0,0])).all():
				dec_seq += u'A'
			elif (x == np.array([0,1,0,0])).all():
				dec_seq += u'C'
			elif (x == np.array([0,0,1,0])).all():
				dec_seq += u'G'
			elif (x == np.array([0,0,0,1])).all():
				dec_seq += u'T'
		return dec_seq

	def __train(self):
		return 0

	def compile(self):
		return 0

	def design(self):
		return 0

	def save(self):
		return 0