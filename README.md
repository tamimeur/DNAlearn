# DNAlearn

DNAlearn is a recommender system for DNA sequence optimization. Users who want to model their high-througput DNA sequencing results can do so easily with DNAlearn. 

DNAlearn builds and tunes an optimal convolutional neural network (CNN) to model DNA sequencing data. Using this model, DNAlearn recommends an optimized sequence - predicted to maximize experimental output. 

## Installation
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
