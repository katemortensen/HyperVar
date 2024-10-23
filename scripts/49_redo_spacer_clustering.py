#!/usr/bin/env python
# coding: utf-8

# %% 

# Author: Kate Mortensen
# Date: 01-11-2024
# Ven Diagram: Spacers from individual assemblies only, co-assembly only, and in common with both individual assemblies and co-assembly
# Purpose: Cluster spacers from all assembly types and assembly methods for "spacer fidelity table"

# In[2]:

import pandas as pd
import numpy as np
import argparse
import os

# In[3]:

# Initialize parser
parser = argparse.ArgumentParser()

# argument for a list of strings
def list_of_strings(arg):
    return arg.split(',')
 
# Adding optional argument
parser.add_argument("-wkdir", "--wkdir", help = "path to project working directory.")
parser.add_argument("-accessions", "--accessions", help = "path to file with list of sample accession numbers.")
parser.add_argument("-mag_dirs", "--mag_dirs", help = "list of paths to MAG subwkdirs", type=list_of_strings)
parser.add_argument("-out_dir", "--out_dir", help = "path to save output (not including file name).")

# Read arguments from command line
args = parser.parse_args()
wkdir = args.wkdir
accessions = args.accessions
mag_dirs = args.mag_dirs
out_dir = args.out_dir

samples = open(accessions, "r").read().strip('\n').split('\n')
file_out = f'{out_dir}/all_spacers.tsv'

print("Output File: ", file_out)
print("MAG dirs: ", mag_dirs)

# comment this out before running pipeline

# %% 
# temp working variables 
wkdir = "/home/kmorten/ADMB"
accessions = f'{wkdir}/SRR_list.txt'
samples = open(accessions, "r").read().strip('\n').split('\n')
mag_dirs = [f'{wkdir}/METAWRAP_MAGS/',
            f'{wkdir}/DASTOOL_MAGS/', 
            f'{wkdir}/MAGSCOT_MAGS/']
out_dir = f'{wkdir}/SPACERS'
file_out = f'{out_dir}/all_spacers.tsv'

# %% 
'''Compile All Spacers'''

if os.path.exists(file_out):
    # delete the file
    os.remove(file_out)

count = 0
with open(file_out, 'a') as fp_out:
    # spacers from co-assembly
    sp_path = f'{wkdir}/CRISPRONE/rNA/coassembly/coassembly-spacer.seq'
    for line in open(sp_path):
        fp_out.write(line)
        count += 1
    # spacers from samples 
    for sample in samples :
        sp_path = f'{wkdir}/CRISPRONE/rNA/{sample}/{sample}-spacer.seq'
        for line in open(sp_path, 'r'):
            if line.startswith('>'):
                old_contig_name = line.strip('>')
                new_contig_name = f'>sample_{old_contig_name}'
                fp_out.write(new_contig_name)
            else:
                fp_out.write(line)
            count += 1
    # spacers from MAG method            
    for mag_dir in mag_dirs:
        method = mag_dir.strip('/').split('/')[-1].strip('/').split('_')[0]
        bins = []
        remove_from_bin_list = samples + ['coassembly', 'unbinned']
        for dir in os.listdir(f'{mag_dir}/CRISPRONE/rNA'):
            bin_tmp = dir.split('/')[-1]
            if bin_tmp not in remove_from_bin_list:
                bins.append(bin_tmp)
        # spacers from bins 
        for bin in bins: 
            sp_path = f'{mag_dir}/CRISPRONE/rNA/{bin}/{bin}-spacer.seq'
            for line in open(sp_path, 'r'):
                if line.startswith('>'):
                    old_contig_name = line.strip('>')
                    new_contig_name = f'>{method}_{old_contig_name}'
                    fp_out.write(new_contig_name)
                else:
                    fp_out.write(line)
                count += 1
                
lines = 0
with open(file_out, 'r') as fp:
    lines = sum(1 for line in fp)
    print('Total Number of lines:', lines)

assert count == lines                    
         
         
