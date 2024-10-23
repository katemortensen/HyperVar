#!/usr/bin/env python
# coding: utf-8

# Author: Kate Mortensen
# Date: 11-2-2023
# Purpose: Plot spacer counts over time (where time is represented by time)

# %% 
# Importing the library
import argparse
import pandas as pd
import os
import matplotlib.pyplot as plt

# %% 
# Initialize parser
parser = argparse.ArgumentParser()
 
# Adding optional argument
parser.add_argument("-accessions", "--accessions", help = "path to directory of accessions", required=True)
parser.add_argument("-spcount_file", "--spcount_file", help = "path to assembly2sp_count_sorted_all.tsv file", required=True)
parser.add_argument("-read_count_file", "--read_count_file", help = "path to sample_reads_count.txt file", required=True)
parser.add_argument("-outdir", "--outdir", help = "dir path to save output", required=True)

# Read arguments from command line
args = parser.parse_args()
accessions = args.accessions
spacer_counts_file = args.spcount_file
read_counts_file = args.read_count_file
outdir = args.outdir

try: 
    os.mkdir(outdir)
except OSError as error: 
    print(error)   

# delete hardcoding before running script

# %% 

# accessions = "/home/kmorten/ADMB/SRR_list.txt"
# spacer_counts_file = f'/home/kmorten/ADMB/results/assembly2sp_count_sorted_all.tsv'
# read_counts_file = f'/home/kmorten/ADMB/CONTIG_RELATIVE_ABUNDANCE/sample_reads_count.txt'
# outdir = "/home/kmorten/ADMB/results/"


# %% 
samples = open(accessions, "r").read().strip('\n').split('\n')

df_spcount = pd.read_csv(f'{spacer_counts_file}', sep = '\t')
df_spcount_sub = df_spcount.loc[df_spcount['assembly_type']=='Individual']

assert len(df_spcount_sub) == len(samples)
for i in samples:
    assert i in list(df_spcount_sub['assembly'])

# %% 
df_samples = pd.DataFrame(samples, columns = ['assembly'])
df_merge = df_samples.merge(df_spcount_sub[['assembly','spacer_count']], on = 'assembly', how = 'left') 


# %% 
x = list(df_merge['assembly'])
y = list(df_merge['spacer_count'])
plt.figure(figsize=(7, 5)) 
plt.bar(x, y, color ='gray', width = 0.4)
# plt.bar(x, y, color ='black', width = 0.4)
# lr = Ridge()
# Ridge(alpha=1.0, copy_X=True, fit_intercept=True, max_iter=None,
#    normalize=False, random_state=None, solver='auto', tol=0.001)
# plt.plot(x, lr.coef_*x+lr.intercept_, color='orange')
plt.xticks(rotation=45)
plt.xlabel("Samples (by time)", fontsize=15)
plt.ylabel("Raw Spacer Count", fontsize=15)
plt.title("Spacer Counts Over Time", fontsize=20)
plt.savefig(f"{outdir}/python_figs/spacer_counts_by_sample.pdf", bbox_inches='tight', format="pdf")  
plt.savefig(f"{outdir}/python_figs/spacer_counts_by_sample.png", bbox_inches='tight', format="png")  
plt.show()


# %%
df_reads = pd.read_csv(f'{read_counts_file}', sep = '\t')
df_reads.rename(columns = {"sample": "assembly"}, inplace = True) 
df_merge2 = df_samples.merge(df_reads[['assembly','total_reads']], on = 'assembly', how = 'left') 


x = list(df_merge2['assembly'])
y = list(df_merge2['total_reads'])
plt.figure(figsize=(7, 5)) 
plt.bar(x, y, color ='gray', width = 0.4)
plt.xticks(rotation=45)
plt.xlabel("Samples (by time)", fontsize=15)
plt.ylabel("Read Count", fontsize=15)
plt.title("Read Counts Over Time", fontsize=20)
plt.savefig(f"{outdir}/python_figs/read_counts_by_sample.pdf", bbox_inches='tight', format="pdf")  
plt.savefig(f"{outdir}/python_figs/read_counts_by_sample.png", bbox_inches='tight', format="png")  
plt.show()
# %%
