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