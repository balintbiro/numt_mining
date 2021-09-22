#import the required modules
import os
import pandas as pd
from subprocess import call


def create_lastal_db(organism_name):
    organism_dir=f'../data/{organism_name}/'
    genome_mask=pd.Series(os.listdir(organism_
