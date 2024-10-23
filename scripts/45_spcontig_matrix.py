#!/usr/bin/env python
# coding: utf-8

# Author: Kate Mortensen
# Date: 10/11/2023
# Purpose: cluster spacer containing contigs using t-SNE

# %% 
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from collections import OrderedDict
from operator import getitem

# %% 
# Initialize parser
parser = argparse.ArgumentParser()
 
# Adding optional argument
parser.add_argument("-accessions", "--accessions", help = "file with list of sample accessions", required=True)
parser.add_argument("-spacer_seq_file", "--spacer_seq_file", help = "path to concat_rNA_spacer.seq", required=True)
parser.add_argument("-spacer_counts_file", "--spacer_counts_file", help = "path to assembly2sp_count_sorted.tsv", required=True)
parser.add_argument("-abund_file", "--abund_file", help = "path to relative_contig_abundance_profiles_crt_<mapping_tool>.tsv", required=True)
parser.add_argument("-topmost_bin_count", "--topmost_bin_count", help = "the number bins with the highest spacer counts (bins ranked in descending order of spacer counts, selecting the highest rank 1st)", required=True)
parser.add_argument("-basename", "--basename", help = "basename of output file", required=True)
parser.add_argument("-outdir", "--outdir", help = "dir path to save outdir", required=True)

# Read arguments from command line
args = parser.parse_args()
accessions = args.accessions
spacer_seq_file = args.spacer_seq_file
spacer_counts_file = args.spacer_counts_file
abund_file = args.abund_file
topmost_bin_count = int(args.topmost_bin_count)
basename = args.basename
outdir = args.outdir
samples = open(accessions, "r").read().strip('\n').split('\n')

try: 
    os.mkdir(f'{outdir}')
except OSError as error: 
    print(error)   

# delete hardcoding before running script

# %%   
# mag_tool = 'METAWRAP'
# mag_tool = 'DASTOOL'
# mag_tool = 'MAGSCOT'


# mapping_tool = 'minimap'
# mapping_tool = 'salmon'

# topmost_bin_count = 5

# accessions = "/home/kmorten/ADMB/SRR_list.txt" 
# samples = open(accessions, "r").read().strip('\n').split('\n')
# spacer_seq_file = f'/home/kmorten/ADMB/{mag_tool}_MAGS/CD-HIT-EST/concat_rNA_spacer.seq'
# spacer_counts_file = f'/home/kmorten/ADMB/{mag_tool}_MAGS/results/assembly2sp_count_sorted.tsv'
# abund_file = f'/home/kmorten/ADMB/{mag_tool}_MAGS/CONTIG_RELATIVE_ABUNDANCE/relative_contig_abundance_profiles_crt_{mapping_tool}.tsv'
# outdir = f'/home/kmorten/ADMB/{mag_tool}_MAGS/results/'
# basename = f'spcontig_matrix_{mapping_tool}'

# %% 
df_spseq = pd.read_csv(f'{spacer_seq_file}', sep = '\t')
df_spcount = pd.read_csv(f'{spacer_counts_file}', sep = '\t')
df_spcount_sub = df_spcount.loc[(df_spcount['assembly_type']=='Co-assembled') & (df_spcount['assembly'] != 'coassembly')]
df_abund = pd.read_csv(f'{abund_file}', sep = '\t')

# %% 
'''bin2spacer_count'''
'''bin > spacers,spacer_count'''
bin2spacer_count = {}
with open(spacer_seq_file, "r") as fp: 
    for line in fp:
        if line[0] == '>':
            assembly = line.strip('>').split('_')[0]
            contig = line.strip('\n').split(assembly)[1].strip('_')
            if bin2spacer_count.get(assembly) == None:
                bin2spacer_count.update({assembly:{'spacers':[], 'spacer_count':0}})
            new_contig_list = [contig] + [i for i in bin2spacer_count[assembly].get('spacers')]
            bin2spacer_count[assembly]['spacers'] = new_contig_list
            new_spacer_count = 1 + bin2spacer_count[assembly].get('spacer_count')
            assert len(set(new_contig_list)) == new_spacer_count 
            bin2spacer_count[assembly]['spacer_count'] = new_spacer_count

# %% 
'''bins with the most spacers'''
top_bins = []
count = 1
for bin in list(df_spcount_sub['assembly']):
    if (bin != 'unbinned') and (count <= topmost_bin_count):
        top_bins.append(bin)
        count += 1

# %% 
'''cross check non-redunant spacer counts'''          
for bin in top_bins:
    count1 = bin2spacer_count[bin].get('spacer_count')
    count2 = int(df_spcount_sub[df_spcount_sub['assembly'] == bin]['spacer_count'])
    assert count1 == count2
 
# %% 
# df_spseq = pd.read_csv(f'{spacer_seq_file}', sep = '\t')
# df_spcount = pd.read_csv(f'{spacer_counts_file}', sep = '\t')
# df_abund = pd.read_csv(f'{abund_file}', sep = '\t')

# # %% 
# # bins with the most spacer (highest spacer counts)
# df_spcount_sub = df_spcount.loc[(df_spcount['assembly_type']=='Co-assembled') & (df_spcount['assembly'] != 'coassembly')]
# top_bins = list(df_spcount_sub['assembly'][0:topmost_bin_count])


# '''bin2spacer_count'''
# '''bin > spacers,spacer_count'''

# bin2spacer_count = {}
# with open(spacer_seq_file, "r") as fp: 
#     for line in fp:
#         if line[0] == '>':
#             assembly = line.strip('>').split('_')[0]
#             contig = line.strip('\n').split(assembly)[1].strip('_')
#             if bin2spacer_count.get(assembly) == None:
#                 bin2spacer_count.update({assembly:{'spacers':[], 'spacer_count':0}})
#             new_contig_list = [contig] + [i for i in bin2spacer_count[assembly].get('spacers')]
#             bin2spacer_count[assembly]['spacers'] = new_contig_list
#             new_spacer_count = 1 + bin2spacer_count[assembly].get('spacer_count')
#             assert len(set(new_contig_list)) == new_spacer_count 
#             bin2spacer_count[assembly]['spacer_count'] = new_spacer_count

# for bin in top_bins:
#     count1 = bin2spacer_count[bin].get('spacer_count')
#     count2 = int(df_spcount_sub[df_spcount_sub['assembly'] == bin]['spacer_count'])
#     assert count1 == count2

# %%
'''assembly2profile'''
'''Weighted & Non-weighted Relative Redistributed Read Counts Per Million '''
'''assembly->contig->sample->wrrcpm ,rrcpm'''

path_in = abund_file
assembly2profile = {}
header_accounted = 0
cols = []

file = open(path_in, 'r')
for line in file:
    linep = line.strip("\n").split("\t")
    if header_accounted == 0:
        cols = linep
        header_accounted += 1
    else:
        lineDict_tmp = dict(zip(cols,linep))
        assembly_type = lineDict_tmp.get('assembly_type')
        assembly = lineDict_tmp.get('assembly')
        contig = lineDict_tmp.get('contig')
        crt = lineDict_tmp.get('crt')
        sample2profile_tmp = {}
        for sample in samples:
            print(sample)
            # rrcpm = float(lineDict_tmp.get(sample))
            if lineDict_tmp.get(sample) == None:
                print(lineDict_tmp.keys())
                print(linep)
            rrcpm = float(lineDict_tmp.get(sample))
            sample2profile_tmp.update({sample:{'rrcpm':rrcpm}})
        # update assembly2profile with entire sample2profile_tmp dictionary
        if assembly2profile.get(assembly) == None:
            assembly2profile[assembly] = {contig:sample2profile_tmp}
        else:
            assembly2profile[assembly][contig] = sample2profile_tmp

# %% 
'matrix for all contig containing spacers in the topmost spacer abundant bins'

for bin in top_bins:
    # for bin in bin2spacer_count.keys():
    path_out = f'{outdir}/{basename}_{bin}.tsv'
    with open(path_out, "w") as fp_out:
        new_row = ['spacer_contigs'] + samples
        new_row_string = "\t".join(new_row).strip("\t") + "\n"
        fp_out.write(new_row_string)    
        spacer_contigs = []
        for sp in list(set(bin2spacer_count[bin].get('spacers'))):
            spacer = sp.split(':')[0]
            spacer_contigs.append(spacer)
        for spacer in list(set(spacer_contigs)):
            new_row = [spacer] + [str(assembly2profile[bin][spacer][sample].get('rrcpm')) for sample in samples]
            new_row_string = "\t".join(new_row).strip("\t") + "\n"
            fp_out.write(new_row_string)


