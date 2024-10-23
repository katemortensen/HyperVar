#!/usr/bin/env python
# coding: utf-8

# In[2]:
# Author: Kate Mortensen
# Last Modified: 9/4/2023
Purpose = "Report spacer counts for assemblies and assembly types."


# In[5]:

import pandas as pd
import numpy as np
import argparse
import os


# In[5]:
# Initialize parser
parser = argparse.ArgumentParser()
 
# Arguments
parser = argparse.ArgumentParser(description=Purpose)
parser.add_argument("-input_files", "--input_files", help = "list of paths to assembly2sp_summary.tsv files (ex: \'/usr/home/subwkdir/results/assembly2sp_count_summary.tsv\')", type=str, required=True)
parser.add_argument("-topmost_bins", "--topmost_bins", help = "top X number most abundant bins")
parser.add_argument("-outdir", "--outdir", help = "path to save output directory", required=True)
parser.add_argument("-output_file", "--output_file", help = "output file name", required=True)
args = parser.parse_args()

# input_files = args.input_files
input_files = [i for i in args.input_files.split(',')]
topmost_bins = args.topmost_bins
outdir = args.outdir
output_file = args.output_file 

try: 
    os.mkdir(outdir) 
except OSError as error: 
    print(error) 
    

# In[6]:

# input_files = ['/home/kmorten/ADMB//METAWRAP_MAGS/results/assembly2sp_summary.tsv',
#  '/home/kmorten/ADMB//DASTOOL_MAGS/results/assembly2sp_summary.tsv',
#  '/home/kmorten/ADMB//MAGSCOT_MAGS/results/assembly2sp_summary.tsv']
# topmost_bins = 5
# outdir = '/home/kmorten/ADMB/results'
# output_file = 'assembly2sp_summary_all.tsv'

# try: 
#     os.mkdir(outdir) 
# except OSError as error: 
#     print(error)  


# #### Spacer Count

# %% 

lines_seen_so_far = set()
line_count = 0
file_out = open(f'{outdir}/{output_file}', 'w')
for input_file in input_files:
    tool = input_file.split('_')[0].split('/')[-1]
    file_in = open(f'{input_file}', 'r')
    for line in file_in:
        if line not in lines_seen_so_far:
            lines_seen_so_far.add(line)       
            if line_count == 0:
                new_line = line.strip('\n') + '\ttool\n'
                file_out.write(new_line)      
            elif line.split('\t')[0] == 'coassembly':
                new_line = line.strip('\n') + '\tMETASPADES\n'
                file_out.write(new_line)      
            elif 'Individual' in line:
                new_line = line.strip('\n') + '\tMETASPADES\n'
                file_out.write(new_line)      
            else: 
                new_line = line.strip('\n') + f'\t{tool}\n'
                file_out.write(new_line)
        line_count += 1
    file_in.close()
file_out.close()
file_out_lines = open(f'{outdir}/{output_file}').readlines()
assert len(file_out_lines) == len(lines_seen_so_far) 
  

# %%
