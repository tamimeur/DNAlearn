import numpy as np 
from keras.models import Sequential
from keras.layers.core import Dense, Dropout, Activation, Flatten
from keras.layers.convolutional import Convolution1D
from keras.optimizers import RMSprop, SGD, Adam
from keras.wrappers.scikit_learn import KerasRegressor
import random
import pandas as pd
from sklearn.model_selection import GridSearchCV
from sklearn.cross_validation import train_test_split
from sklearn import preprocessing
from sklearn.metrics import make_scorer
import scipy.stats as stats
import logging

np.random.seed(1337)

#this Needs to change
#seq_len = 150 

def create_model(learn_rate=0.001, seq_len=15):
	########## Build CNN Model #############
	cnn = Sequential()
	cnn.add(Convolution1D(nb_filter=30,filter_length=6,input_dim=4,input_length=seq_len,border_mode="same", activation='relu'))
	cnn.add(Dropout(0.1))
	cnn.add(Convolution1D(nb_filter=40,filter_length=6,input_dim=4,input_length=seq_len,border_mode="same", activation='relu'))

	cnn.add(Flatten())

	#this dense layer should be tuned to various percentages
	#of the first input layer, like [(seq_len*2)/3, seq_len/3]
	cnn.add(Dense(40))
	cnn.add(Dropout(0.2))
	cnn.add(Activation('relu'))

	cnn.add(Dense(1))
	cnn.add(Activation('linear'))

	#compile the model
	adam = Adam(lr=learn_rate, beta_1=0.9, beta_2=0.999, epsilon=1e-08)
	#model.compile(loss='mean_squared_error', optimizer=adam)
	rms = RMSprop(lr=learn_rate, rho=0.9, epsilon=1e-08)

	#optimizer, lr, and momeuntum should be tuned
	#need a custom metric (r^2 methinks)
	cnn.compile(loss='mean_squared_logarithmic_error', optimizer=rms)#, metrics=['accuracy'])
	#print 'Model compiled in {0} seconds'.format(time.time() - start_time)
	return cnn

def my_custom_r2_func(ground_truth, predictions):
	slope, intercept, r_value, p_value, std_err = stats.linregress(ground_truth.reshape(-1),predictions.reshape(-1))
	print "Evaluating a model meow: ", r_value**2
	return r_value**2
	# diff = np.abs(ground_truth - predictions).max()
	# return np.log(1 + diff)

class dnaModel(object):
	""" A CNN model for DNA sequences """
	def __init__(self, input_file):
		#The input_file should contain only two columns: sequence + output
		#this should create a CNN/model object which can be trained,tested optimally (hyperparam opt)
		#logging.basicConfig(filename='optseq.log',level=logging.DEBUG)
		#logging.info("test")
		self.model, self.X_train, self.Y_train, self.X_test, self.Y_test, self.batch, self.epoch = \
				self.__opt_model(input_file)
		

	###inits/bulids a new model and optimizes hyperparams
	###trains/fits and saves hyperparams of training to dnaModel members
	###currently a simple gridsearch, but should switch to a stochaistic 
	###optimization alg with mean performance as the cost fcn
	def __opt_model(self, input_file):
		xtrain, ytrain, xtest, ytest = [], [], [], []

		raw_data = pd.read_excel(input_file, header=0, parse_cols="A,B")
		print raw_data.columns

		###########  Cleaning output column  ###################
		df = raw_data[np.isfinite(raw_data[u' expression'])]

		############# Format NN inputs #############
		seq_len = len(df['sequence'][0])
		X_data = np.empty([len(df),seq_len,4])
		indx = 0

		Y_data = np.array(df[[u' expression']])

		for seq in df[u'sequence']:
			X_data[indx] = self.__oneHotEncoder(seq)
			indx += 1

		# print "CNN input \n", X_data
		# print "CNN output \n", Y_data

		########## RANDOM TEST/TRAIN SPLIT #########
		normed_out = preprocessing.StandardScaler().fit_transform(Y_data)
		xtrain, xtest, ytrain, ytest = train_test_split(X_data, normed_out, test_size=0.15, random_state=42)
		# xtrain, xtest, ytrain, ytest = train_test_split(X_data, Y_data, test_size=0.15, random_state=42)
		# norm_train = preprocessing.StandardScaler().fit_transform(ytrain)

		model = KerasRegressor(build_fn=create_model, nb_epoch=6, batch_size=128, verbose=1)
		# define the grid search parameters
		num_bases = [seq_len]
		learn_rate = [0.001, 0.01, 0.1]
		#momentum = [0.0, 0.2, 0.4, 0.6, 0.8, 0.9]
		param_grid = dict(learn_rate=learn_rate, seq_len = num_bases)#, momentum=momentum)
		#specify my own scorer for GridSearchCV that uses r2 instead of the estimator's scorer
		grid = GridSearchCV(estimator=model, param_grid=param_grid, scoring=make_scorer(my_custom_r2_func, greater_is_better=True), n_jobs=1)
		grid_result = grid.fit(xtrain, ytrain)
		# summarize results
		print("Best: %f using %s" % (grid_result.best_score_, grid_result.best_params_))
		means = grid_result.cv_results_['mean_test_score']
		stds = grid_result.cv_results_['std_test_score']
		params = grid_result.cv_results_['params']
		for mean, stdev, param in zip(means, stds, params):
		    print("%f (%f) with: %r" % (mean, stdev, param))

		########### Train CNN model ############
		#batch_size and epoch should be tuned
		#cnn.fit(xtrain, norm_train, batch_size=128, nb_epoch=6, verbose=1)

		return model, xtrain, ytrain, xtest, ytest, 128, 6

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