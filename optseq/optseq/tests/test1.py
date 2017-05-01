from unittest import TestCase
import optseq.dnamodel as dm
import pandas as pd
import numpy as np
import os

class TestDnaModel(TestCase):
	def test_load_model(self):
		testpath = os.path.dirname(os.path.abspath(__file__))
		raw_data = pd.read_excel(testpath+'/test_input.xlsx', header=0)
		col_names = ['sequence']
		for i in range(1, len(raw_data.columns)):
			col_names.append('output'+str(i))
		raw_data.columns = col_names
		df = raw_data[np.isfinite(raw_data['output1'])]
		dnaCNN = dm.dnaModel(df, filename=testpath+'/test_model.h5')

#test onehotencoder and onehotdecoder


