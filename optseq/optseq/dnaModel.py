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
		self.raw_data = pd.read_excel(input_file, header=0, parse_cols="A,B")
		self.X_train = []
		self.Y_train = []
		self.X_test = []
		self.Y_test = []
		self.layers = 3

		self.create_model()

	def create_model(self):
		df = self.raw_data[np.isfinite(self.raw_data[u' expression'])]
		return 0

	def train(self):
		return 0