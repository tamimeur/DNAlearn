import click
import os, sys
import errno
import pandas as pd
import numpy as np
import dnamodel as dm
#import logging


all_colors = 'black', 'red', 'green', 'yellow', 'blue', 'magenta', \
             'cyan', 'white'

@click.command()
@click.argument('input_file')
#@click.option('--max', type=click.File('rb'), help='Maximize the last column in the input file')
#@click.option('--min', '-c', type=click.File('rb'), help='Minimize the last column in the input file')


def main(input_file):
    """Automated, experimentally-driven design and optimization of DNA sequences."""
    greet = 'Hello'
    #click.echo(click.style('{0}, {1}'.format(greet, 'OptSeq users!'), fg='blue', blink=False))
    #click.echo(input_file)
    print input_file


    #there has to be an option for a saved, pre-existing model (for iterating)
    raw_data = pd.read_excel(input_file, header=0, parse_cols="A,B")
    print raw_data.columns
    df = raw_data[np.isfinite(raw_data[u' expression'])]

    dnaCNN = dm.dnaModel(df) 
    ### dnaCNN is now a model object with a CNN and training/testing data ###
    ### it compiles ### 
    ### and trains (default my current hyperparams, but this should ###
    ### eventually tune/optimize hyperparams) the CNN ### 
    
    # dnaCNN.design()
    ### when we execute design (default maximizes col B), ###
    ### it returns a set (default some percentage of the input data) ###
    ### of new designs to test (default cheap). ###

    ### I should get this to work with a database, so that ###
    ### models are stored, saved, and updated in the db, that ###
    ### that way we can iterate without having to dl and reload ###
    ### the dnaModel. ###

    dnaCNN.save()
    ### BUT FOR NOW ... the last thing it does is output a saved model file. ###

    dnaCNN.test()

    dnaCNN.design()
