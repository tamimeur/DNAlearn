import click
import os, sys
import errno
import logging
import dnamodel as dm

all_colors = 'black', 'red', 'green', 'yellow', 'blue', 'magenta', \
             'cyan', 'white'

@click.command()
#@click.option('--max', type=click.File('rb'), help='Maximize the last column in the input file')
#@click.option('--min', '-c', type=click.File('rb'), help='Minimize the last column in the input file')


def main(input_file):
    """Automated, experimentally-driven design of genetic circuits."""
    greet = 'Hello'
    click.echo(click.style('{0}, {1}'.format(greet, 'OptSeq users!'), fg='blue', blink=False))
    
    #logging.basicConfig(filename='optseq.log',level=logging.DEBUG)
    #logging.info('Welcome to OptSeq!')

    #there has to be an option for a saved, pre-existing model (for iterating)
    # dnaCNN = dm.dnaModel(input_file)
    # dnaCNN.compile()
    # dnaCNN.design()
