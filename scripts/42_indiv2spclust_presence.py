#!/usr/bin/env python
# coding: utf-8

# Author: Kate Mortensen
# Date: 9-6-2023
# Purpose: To populate a df with spacer presence in individual assemblies for clustering 

# %% 
import pandas as pd
import numpy as np
import argparse

# %% 
# Initialize parser
parser = argparse.ArgumentParser()
 
# Adding optional argument
parser.add_argument("-input", "--input", help = "path to cd-hit-est-#_spacer.out.clstr", required=True)
parser.add_argument("-output", "--output", help = "path to save output including file name", required=True)

# Read arguments from command line
args = parser.parse_args()
input = args.input
rep = input.strip('_spacer.out.clstr').split('-')[-1]
output = args.output

# %%
input = "/home/kmorten/ADMB/CD-HIT-EST/cd-hit-est-rNA_spacer.out.clstr"
output = "/home/kmorten/ADMB/results/indiv2spclust_presence.tsv"
output = "/home/kmorten/ADMB/results/indiv2spclust_presence_matrix.tsv"

# %% 
'''collect all spacer clusters from Individual assemblies'''

assembly2clust = {}
all_spacers = []
pseudo_spacers = []
with open(input, 'r') as fp:
    for line in fp:
        if line[0] == '>':
            linep=line.strip('\n').strip('>').split(' ')
            pseudo_spacer=f'{linep[0]}_{linep[1]}'
            pseudo_spacers.append(pseudo_spacer)
        else:
            assembly=line.split('>')[1].split('_')[0]
            contig=line.split(f'{assembly}_')[1].split('...')[0]
            all_spacers.append(contig)
            if assembly2clust.get(assembly) == None:
                assembly2clust[assembly] = {pseudo_spacer:1}
            else:
                assembly2clust[assembly][pseudo_spacer] = 1
assert len(all_spacers) == len(set(all_spacers))
            

# %% 
with open(output, 'w') as fp_out:
    header = 'assembly\tpseudo_spacer\tpresence\n'
    fp_out.write(header)
    for assembly in assembly2clust.keys():
        for pseudo_spacer in pseudo_spacers:
            if assembly2clust[assembly].get(pseudo_spacer) == None:
                presence = 0
            else:
                presence = 1
            fp_out.write(f'{assembly}\t{pseudo_spacer}\t{presence}\n')
                

# %%
with open(output, 'w') as fp_out:
    header_list = ['pseudo_spacer']
    for i in assembly2clust.keys():
        header_list.append(i)
    header = ('\t').join(header_list)
    fp_out.write(f'{header}\n')
    for pseudo_spacer in pseudo_spacers:
        list_of_presence = []
        for assembly in assembly2clust.keys():
            if assembly2clust[assembly].get(pseudo_spacer) == None:
                presence = '0'
            else:
                presence = '1'
            list_of_presence.append(presence)
        new_row_list= [pseudo_spacer] + list_of_presence
        new_row = ('\t').join(new_row_list)
        fp_out.write(f'{new_row}\n')
                



# %%
