# DNAlearn

DNAlearn is a recommender system for DNA sequence optimization. Users who want to model their high-througput DNA sequencing results can do so easily with DNAlearn. 

DNAlearn builds and tunes an optimal convolutional neural network (CNN) to model DNA sequencing data. Using this model, DNAlearn recommends an optimized sequence - predicted to maximize experimental output. 

## Installation
I recommend installing in a virtual environment, as DNAlearn will attempt to install a specific version of Keras(1.2) and various other packages.
	
	pip install dnalearn

## Usage
DNAlearn requires an input xlsx file of DNA sequencing data with the following structure:

Example:

| Sequences     | Output        |
| ------------- |:-------------:|
| actgactg ...  | 3      |
| actgactg ...  | 5      |



To generate a new model using the command line interface:

	dnalearn input.xlsx

To retrain a previously generated model using the command line interface:

	dnalearn input.xlsx -m model.h5

DNAlearn provides the option to tune the following hyperparameters:

	loss  (0='mean_squared_error'; 1='mean_squared_logarithmic_error'; 2='poisson')
	optimizer (0='rms'; 1='adam'; 2='sgd')
	lr (learning rate)
	dense1_neurons (number of neurons in the first dense layer)
	dense2_neurons (number of neurons in the second dense layer - 0 means no 2nd dense layer)
	filter_batch (the filter batch size for convolutional layers)
	filter_len (the length of the filters for convultional layers)

To specify hyperparameter values to vary, you must provide a config.json file using the following command:

	dnalearn input.xlsx -c config.json

Below is an example config file:

	{
	"optimizer": [0,1],
	"loss": [0,2],
	"lr": [0.001],
	"dense1_neurons": [50,90],
	"dense2_neurons": [0],
	"filter_batch": [20],
	"filter_len": [6,10]
	}

By combinatorially constructing models using these hyperparameters, DNAlearn will test (in this example) 16 different models and return the best one. Be careful - providing too many hyperparameters to choose from, and therefore a large amount of models to test, will take a large amount of computation time.



Have fun!
