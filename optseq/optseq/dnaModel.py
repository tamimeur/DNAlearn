import numpy as np 
from keras.models import Sequential
from keras.layers.core import Dense, Dropout, Activation, Flatten
import random
import pandas as pd
from sklearn.cross_validation import train_test_split
from sklearn import preprocessing
import scipy.stats as stats
import logging

np.random.seed(1337)

class dnaModel(object):
	""" A CNN model for DNA sequences """
	def __init__(self, input_file):
		#The input_file should contain only two columns: sequence + output
		#this should create a CNN/model object which can be trained,tested optimally (hyperparam opt)
		#logging.basicConfig(filename='optseq.log',level=logging.DEBUG)
		#logging.info("test")
		self.model, self.X_train, self.Y_train, self.X_test, self.Y_test, self.batch, self.epoch = \
				self.__create_model(input_file)
		

	###inits/bulids a new model and optimizes hyperparams
	###trains/fits and saves hyperparams of training to dnaModel members

	def __create_model(self, input_file):
		xtrain, ytrain, xtest, ytest = [], [], [], []

		raw_data = pd.read_excel(input_file, header=0, parse_cols="A,B")
		print raw_data.columns

		###########  Cleaning output column  ###################
		df = raw_data[np.isfinite(raw_data[u' expression'])]

		############# Format NN inputs #############
		seq_len = len(df['sequence'][0])
		X_data = np.empty([len(df),150,4])
		indx = 0

		Y_data = np.array(df[[u' expression']])

		for seq in df[u'sequence']:
			X_data[indx] = self.__oneHotEncoder(seq)
			indx += 1

		# print "CNN input \n", X_data
		# print "CNN output \n", Y_data

		########## RANDOM TEST/TRAIN SPLIT #########
		xtrain, xtest, ytrain, ytest = train_test_split(X_data, Y_data, test_size=0.15, random_state=42)
		norm_train = preprocessing.StandardScaler().fit_transform(ytrain)

		########## Build CNN Model #############
		cnn = Sequential()
		cnn.add(Convolution1D(nb_filter=30,filter_length=6,input_dim=4,input_length=150,border_mode="same", activation='relu'))
		cnn.add(Dropout(0.1))
		cnn.add(Convolution1D(nb_filter=40,filter_length=6,input_dim=4,input_length=150,border_mode="same", activation='relu'))

		cnn.add(Flatten())

		cnn.add(Dense(40))
		cnn.add(Dropout(0.2))
		cnn.add(Activation('relu'))

		cnn.add(Dense(1))
		cnn.add(Activation('linear'))

		#compile the model
		adam = Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-08)
		#model.compile(loss='mean_squared_error', optimizer=adam)
		rms = RMSprop(lr=0.001, rho=0.9, epsilon=1e-08)
		cnn.compile(loss='mean_squared_logarithmic_error', optimizer=rms)
		#print 'Model compiled in {0} seconds'.format(time.time() - start_time)

		########### Train CNN model ############
		cnn.fit(X_train, norm_train, batch_size=128, nb_epoch=6, verbose=1)

		return cnn, xtrain, ytrain, xtest, ytest, 128, 6

	def __oneHotEncoder(self,seq):
		base_dict = {u'A':[1,0,0,0],u'C':[0,1,0,0],u'G':[0,0,1,0],u'T':[0,0,0,1]}
		return np.array([base_dict[x] for x in seq])

	def __oneHotDecoder(self,encseq):
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

	#### this just does a replicate fit with saved previous keras and dnaModel ... models
	def __train(self):
		return 0

	def compile(self):
		return 0

	def design(self):
		return 0

	def save(self):
		return 0