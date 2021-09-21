#importing the required modules
import os
import pandas as pd
import urllib.request   

#create repository for every organisms
def mkdir(organism_name):
    parent_dir='../results/'
    path=os.path.join(parent_dir+organism_name)
    os.mkdir(path)