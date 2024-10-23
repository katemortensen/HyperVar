#!/usr/bin/env python
# coding: utf-8


# Author: Kate Mortensen
# Date: 01-15-2024
# Ven Diagram: Spacers from individual assemblies only, co-assembly only, and in common with both individual assemblies and co-assembly
# Purpose: After finding representative repeat-associated spacers, clustering these repeat-associated spacers, find all individual assemblies and bins that are associated with each representative spacer for clustered spacers. 

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
parser.add_argument("-input", "--input", help = "path to cd-hit-est-rNA_spacer.out.clstr.")
parser.add_argument("-accessions", "--accessions", help = "path to file with list of sample accession numbers.")
parser.add_argument("-mag_dirs", "--mag_dirs", help = "list of paths to MAG subwkdirs", type=list_of_strings)
parser.add_argument("-out_dir", "--out_dir", help = "path to save output (not including file name).")

# Read arguments from command line
args = parser.parse_args()
input = args.input
accessions = args.accessions
mag_dirs = args.mag_dirs
out_dir = args.out_dir

samples = open(accessions, "r").read().strip('\n').split('\n')
methods = [mag_dir.strip('/').split('/')[-1].strip('/').split('_')[0] for mag_dir in mag_dirs]
file_out = f'{out_dir}/sp2venn_all.tsv'

print("Output File: ", file_out)
print("MAG dirs: ", mag_dirs)

# comment this out before running pipeline


# %% 
# temp working variables 

# wkdir = "/home/kmorten/ADMB"
# mag_dirs = [f'{wkdir}/METAWRAP_MAGS/',
#             f'{wkdir}/DASTOOL_MAGS/', 
#             f'{wkdir}/MAGSCOT_MAGS/']

# input = f'{wkdir}/SPACERS/cd-hit-est-rNA_spacer.out.clstr'
# methods = [mag_dir.strip('/').split('/')[-1].strip('/').split('_')[0] for mag_dir in mag_dirs]
# out_dir = f'{wkdir}/SPACERS'
# file_out = f'{out_dir}/sp2venn_all.tsv'


# %% 

'''clustDict'''
'''cluster > methods > assemblies > contigs'''
'''Dictionary parsing CD-HIT-EST output of spacer clusters'''

clustDict = {}
count = -1
sp_count = 0

with open(input) as fp:
    lines = fp.readlines()
    for line in lines:
        if line.startswith(">") :
            count += 1
            clust_name = "clust_%s" % (count)
        if line[0] != '>':
            method = line.split(">")[1].split("_")[0]
            sp_count += 1
            if method == 'coassembly':
                assembly = method
            else:
                assembly = line.split(method)[1].split('_')[1]
            contig = line.strip('\n').split(assembly)[1].strip('_').split('...')[0]
            if clustDict.get(clust_name) == None:
                clustDict[clust_name] = {method:{assembly:[contig]}}
            elif clustDict[clust_name].get(method) == None:
                clustDict[clust_name][method] = {assembly:[contig]}
            elif clustDict[clust_name][method].get(assembly) == None:
                clustDict[clust_name][method][assembly] = [contig]
            elif clustDict[clust_name][method].get(assembly) != None:
                old_contig_list = clustDict[clust_name][method].get(assembly)
                new_contig_list = old_contig_list + [contig]
                clustDict[clust_name][method][assembly] = new_contig_list


# %% 

'''NRdict'''
'''Non-redundant Spacer Dictionary'''
'''I,C,CI and M,IM,CM,CIM > clust counts, mag methods > clust_counts'''
'''representative spacer (aka non-redundant (NR) aka cluster) > I,C,M,IC,IM,CM,CIM > method > assemblies > contigs'''
'''Dictionary of Spacer Clusters and the Assemblies for the spacers in a cluster originate'''

fidelities = ['I','C','M','CI','IM','CM','CIM']

NRdict = {}
for clust in clustDict.keys():
    methods_tmp = list(set(clustDict.get(clust).keys()))
    lfidelity = []
    if 'sample' in methods_tmp:
        lfidelity.append('I') # I = individual sample
        methods_tmp.remove('sample')
    if 'coassembly' in methods_tmp: 
        lfidelity.append('C') # C = coassembly
        methods_tmp.remove('coassembly')
    if len(methods_tmp) > 0:
        lfidelity.append('M') # M = mags
    lfidelity.sort()
    fidelity = ("").join(lfidelity)
    if 'M' not in fidelity:
        if NRdict.get(fidelity) == None:
            NRdict[fidelity] = 1
        else:
            new_count = NRdict.get(fidelity) + 1
            NRdict[fidelity] = new_count
    elif 'M' in fidelity:
        if NRdict.get(fidelity) == None:
            NRdict[fidelity] = dict(zip(methods, [0]*len(methods)))
        for m in methods:
            if m in methods_tmp:
                new_count = NRdict[fidelity].get(m) + 1
                NRdict[fidelity][m] = new_count


for fidelity in NRdict.keys():
    if 'M' not in fidelity:
        tmp_count = NRdict.get(fidelity)
        print(f'{fidelity} NR spacers: {tmp_count}')
    elif 'M' in fidelity:
        for m in NRdict[fidelity].keys():
            tmp_count = NRdict[fidelity].get(m)
            print(f'{fidelity}-{m} NR spacers: {tmp_count}')
    

# %% 

'''dataframe of non-redundant spacer counts'''

df_nr = pd.DataFrame({'Fidelity' : ['I','C','M','IC','IM','CM','ICM']})
col_count = 1

for m in methods:
    fcounts = []
    for fidelity in list(df_nr['Fidelity']):
        f = ''.join(sorted(fidelity))
        tmp_count = 0
        if 'M' not in f:
            tmp_count = NRdict.get(f)
            fcounts.append(tmp_count)
        elif 'M' in f:
            tmp_count = NRdict[f].get(m)
            fcounts.append(tmp_count)
        if tmp_count == None:
            print(m, ' ', f)
    new_col_name = 'NR_SP_' + m # non-redundant spacers (NR_SP)
    df_nr.insert(col_count, new_col_name, fcounts, True)
    col_count += 1

mlist = []
for f in df_nr['Fidelity'] :
    if 'M' in f:
        mlist.append(f)


df_sub = df_nr[df_nr['Fidelity'].isin(mlist)]
tot = df_sub.loc[ : ,df_sub.columns != 'Fidelity'].sum()
df_nr.loc[len(df_nr.index)] = ['Total M'] + list(tot)
df_nr.to_csv(file_out, sep = '\t', index=False)

df_nr

