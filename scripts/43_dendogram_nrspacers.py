#!/usr/bin/env python
# coding: utf-8

# Author: Kate Mortensen
# Date: 9-6-2023
# Purpose: To cluster spacer presence in individual assemblies in a bidirectional heirchical fashion 

#%% 

# Installations
# pip install scikit-learn
# pip install sunbird

# %% 
# Importing the library
import argparse
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

# %% 
# Initialize parser
parser = argparse.ArgumentParser()
 
# Adding optional argument
parser.add_argument("-input", "--input", help = "path to indiv2spclust_presence.tsv", required=True)
parser.add_argument("-n", "--n", help = "number of spacers to randomly sample", type=int, required=True)
parser.add_argument("-outdir", "--outdir", help = "dir path to save output", required=True)

# Read arguments from command line
args = parser.parse_args()
input = args.input
n = int(args.n)
outdir = args.outdir

try: 
    os.mkdir(outdir)
except OSError as error: 
    print(error)   

# delete hardcoding before running script

# %% 
# input = "/home/kmorten/ADMB/results/indiv2spclust_presence_matrix.tsv"
# outdir = "/home/kmorten/ADMB/results/python_figs"
# n = 1000

# %% 
df = pd.read_csv(input, sep='\t')
samples = list(df.columns)[1:]
col_pseudo_spacers = list(df['pseudo_spacer'])

# %% 

# swap rows and columns 
dfT = df[samples].T
# randomly sample n pseudo_clusters
dfT = dfT.sample(n=n, axis='columns')
# sort by sample 
dfT = dfT.sort_index()

# %%
cmap = sns.cm.rocket_r
# # Define colors
# colors = ["white", "black"]
# cmap = LinearSegmentedColormap.from_list('Custom', colors, len(colors))

#%%
fg = sns.clustermap(
    dfT,
    figsize=(10, 5),
    row_cluster=False, 
    xticklabels=False, 
    cmap = cmap,
    cbar_pos = None
)

fg.fig.suptitle('Heirarchical Clustering of Non-redundant Spacers', 
                fontsize = 18, y = 1.05)
ax = fg.ax_heatmap 
ax.set_xlabel("Spacers", fontsize=15)
ax.set_ylabel("Samples Across Time", fontsize=15)

# %% 
fg.savefig(f'{outdir}/dendogram_nrspacers.pdf', format="pdf", bbox_inches="tight")
fg.savefig(f'{outdir}/dendogram_nrspacers.png', format="png", bbox_inches="tight")

# %% 

# plt.savefig(f'{outdir}/dendogram_nrspacers.pdf', format="pdf", bbox_inches="tight")
# plt.savefig(f'{outdir}/dendogram_nrspacers.png', format="png", bbox_inches="tight")

# plt.show()
# %%

