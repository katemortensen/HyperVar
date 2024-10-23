#!/usr/bin/env python
# coding: utf-8

# %% 

# Author: Kate Mortensen
# Date: 5-9-2023
# Purpose: data file for spacer count bar plot with respect to repeat guide, assembly type, assembly name 

# In[2]:

import pandas as pd
import numpy as np
import argparse
import os


# In[3]:

# Initialize parser
parser = argparse.ArgumentParser()
 
# Adding optional argument
parser.add_argument("-cdhit_dir", "--cdhit_dir", help = "path to directory where CD-HIT-EST output exists (ex: \'/usr/home/wkdir/CD-HIT-EST\')", required=True)
parser.add_argument("-repeat_guides", "--repeat_guides", help = "path to repeat guides fasta (/wkdir/CD-HIT-EST/top25_representative_repeats.fa)", required=True)
parser.add_argument("-outdir", "--outdir", help = "path to save output directory", required=True)

# Read arguments from command line
args = parser.parse_args()
cdhit_dir = args.cdhit_dir
repeat_guides = args.repeat_guides
outdir = args.outdir

try: 
    os.mkdir(outdir) 
except OSError as error: 
    print(error) 
    
# %% 
# cdhit_dir = '/home/kmorten/ADMB/DASTOOL_MAGS/CD-HIT-EST'
# # mags_dir = '/home/kmorten/ADMB/DASTOOL_MAGS/bins'
# # accessions = '/home/kmorten/ADMB/SRR_list.txt'
# repeat_guides = "/home/kmorten/ADMB/DASTOOL_MAGS/CD-HIT-EST/top25_representative_repeats.fa"
# outdir = "/home/kmorten/ADMB/DASTOOL_MAGS/results"
# nodeIDs_file = '/home/kmorten/ADMB/DASTOOL_MAGS/CONTIG_RELATIVE_ABUNDANCE/nodeIDs.txt'

try: 
    os.mkdir(outdir) 
except OSError as error: 
    print(error) 

# %% 
repr_reps = ['rNA', 'rALL']
with open(repeat_guides, 'r') as fp:
    lines = fp.readlines()
    for line in lines:
        if line.startswith(">"):
            repr_reps.append(line.strip(">").strip("\n").strip(" "))

# %% 
all_assemblies = []
with open(nodeIDs_file, 'r') as fp:
    next(fp)
    for line in fp:
        all_assemblies.append(line.split('\t')[0])
all_assemblies = set(all_assemblies)

# %% 
'''rep2div2spcount'''
'''repeat_guide > diversity_of_assembly_type(s) > spacer_count'''
'''Dictionary of Counts of Unique Spacers shared by different assembly types'''

rep2assembly2spcount = {}
for rep in repr_reps:
    path = f'{cdhit_dir}/cd-hit-est-{rep}_spacer.out.clstr'        
    with open(path, "r") as fp:
        lines = fp.readlines()
        