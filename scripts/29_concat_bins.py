#!/usr/bin/env python
# coding: utf-8

# %% 

# Author: Kate Mortensen
# Last Modified: 5/15/2023

# In[2]:

import pandas as pd
import numpy as np
import argparse
import os
import glob

# In[3]:

# Initialize parser
parser = argparse.ArgumentParser()
 
# Adding optional argument
parser.add_argument("-mag_dir", "--mag_dir", help = "path to MAGs directory")
parser.add_argument("-outdir", "--outdir", help = "path to output directory")

# Read arguments from command line
args = parser.parse_args()
mag_dir = args.mag_dir
outdir = args.outdir

try: 
    os.mkdir(outdir) 
except OSError as error: 
    print(error)  


inputs = glob.glob(f"{mag_dir}/*")
with open(f"{outdir}/bins_concat.fasta", "w") as output:
  first = True
  for file in inputs:
    with open(file, "r") as inputfile:
        for line in inputfile:
            output.write(line )


'''node2count'''            
'''Dictionary of node IDs to occurrences to make sure there are no duplicate node IDs across mags.'''
'''node > count'''

node2count = {}
with open(f"{outdir}/bins_concat.fasta", "r") as fp:
    for line in fp:
        if line[0] == '>':
            node = line
            if node2count.get(node) == None:
                node2count[node] = 1
            else:
                new_count = node2count.get(node) + 1
                node2count[node] = new_count

'''assert no redundancy in node IDs across MAGs'''

for node in node2count.keys():
    assert node2count.get(node) == 1
    

    
    
