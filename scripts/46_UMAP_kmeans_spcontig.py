#!/usr/bin/env python
# coding: utf-8

# Author: Kate Mortensen
# Date: 10/11/2023
# Last Modified: 11/3/2023
# Purpose: cluster spacer containing contigs using k-means and visualize using UMAP


# %%
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from kneed import DataGenerator, KneeLocator
from sklearn.decomposition import PCA
import umap.umap_ as umap


# %% 
# Initialize parser
parser = argparse.ArgumentParser()
 
# Adding optional argument
parser.add_argument("-spacer_seq_file", "--spacer_seq_file", help = "path to concat_rNA_spacer.seq", required=True)
parser.add_argument("-spacer_counts_file", "--spacer_counts_file", help = "path to assembly2sp_count_sorted.tsv", required=True)
parser.add_argument("-abund_file", "--abund_file", help = "path to relative_contig_abundance_profiles_crt_<mapping_tool>.tsv", required=True)
parser.add_argument("-topmost_bin_count", "--topmost_bin_count", help = "the number bins with the highest spacer counts (bins ranked in descending order of spacer counts, selecting the highest rank 1st)", required=True)
parser.add_argument("-avg_abund_file", "--avg_abund_file", help = "with RRCPM abundances of contigs averaged over samples", required=True)
parser.add_argument("-mapping_tool", "--mapping_tool", help = "salmon or minimap etc.", required=True)
parser.add_argument("-outdir", "--outdir", help = "dir path to save outdir", required=True)

# Read arguments from command line
args = parser.parse_args()
spacer_seq_file = args.spacer_seq_file
spacer_counts_file = args.spacer_counts_file
abund_file = args.abund_file
topmost_bin_count = int(args.topmost_bin_count)
avg_abund_file = args.avg_abund_file 
mapping_tool = args.mapping_tool
outdir = args.outdir

try: 
    os.mkdir(f'{outdir}')
except OSError as error: 
    print(error)   

# delete hardcoding before running script

# %%  

# mag_tool = 'METAWRAP'
# # mag_tool = 'DASTOOL'
# # # mag_tool = 'MAGSCOT'

# mapping_tool = 'minimap'
# # mapping_tool = 'salmon'

# topmost_bin_count = 5
# #accessions = "/home/kmorten/ADMB/SRR_list.txt" 
# #samples = open(accessions, "r").read().strip('\n').split('\n')
# spacer_seq_file = f'/home/kmorten/ADMB/{mag_tool}_MAGS/CD-HIT-EST/concat_rNA_spacer.seq'
# spacer_counts_file = f'/home/kmorten/ADMB/{mag_tool}_MAGS/results/assembly2sp_count_sorted.tsv'
# abund_file = f'/home/kmorten/ADMB/{mag_tool}_MAGS/CONTIG_RELATIVE_ABUNDANCE/relative_contig_abundance_profiles_crt_{mapping_tool}.tsv'
# avg_abund_file = f'/home/kmorten/ADMB/{mag_tool}_MAGS/CONTIG_RELATIVE_ABUNDANCE/contig_abund_{mapping_tool}.tsv' 
# outdir = f'/home/kmorten/ADMB/{mag_tool}_MAGS/results/'

# %% 
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

'''bins with the most spacers'''
top_bins = []
count = 1
for bin in list(df_spcount_sub['assembly']):
    if (bin != 'unbinned') and (count <= topmost_bin_count):
        top_bins.append(bin)
        count += 1

'''cross check non-redunant spacer counts'''          
for bin in top_bins:
    count1 = bin2spacer_count[bin].get('spacer_count')
    count2 = int(df_spcount_sub[df_spcount_sub['assembly'] == bin]['spacer_count'])
    assert count1 == count2
    
# %% 
for bin in top_bins:
    
    spcontig_matrix = f'{outdir}/spcontig_matrix_{mapping_tool}_{bin}.tsv'
    
    if os.path.exists(spcontig_matrix) == False:
        print('matrix file not found at path:', spcontig_matrix)
    
    df = pd.read_csv(f'{spcontig_matrix}', sep = '\t')
    assert len(df['spacer_contigs']) == len(set(df['spacer_contigs']))
    
    df_avg_abund = pd.read_csv(f'{avg_abund_file}', sep = '\t')
    df.rename(columns = {"spacer_contigs": "contig"}, inplace = True) 
    df_merge = df.merge(df_avg_abund[['contig', 'arrcpm']], on = 'contig', how = 'left') 
    
    if len(df['contig']) <= 1:
        print(f'mag ({bin}) needs at least 2 contigs with spacers.')
    else :
        try :
            
            # K-means
            scaled_df = StandardScaler().fit_transform(df.drop("contig", axis=1))

            #initialize kmeans parameters
            kmeans_kwargs = {
            "init": "random",
            "n_init": 10,
            "random_state": 1 
            }
            
            #create list to hold SSE values for each k
            sse = []
            for k in range(1, 11):
                kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
                kmeans.fit(scaled_df)
                sse.append(kmeans.inertia_)

            # elbow plot 
            kneedle = KneeLocator(range(1, 11), sse, S=1.0, curve="convex", direction="decreasing")
            # kneedle.plot_knee_normalized()
            kneedle.plot_knee(figsize=(7,5))
            plt.title(f'Knee Point - {bin}', fontsize = 20)
            plt.xlabel("Number of Clusters", fontsize = 15)
            plt.ylabel("SSE", fontsize = 15)
            plt.legend(["data", "knee/elbow"], fontsize="15", loc ="upper right")
            plt.savefig(f"{outdir}/python_figs/elbow_UMAP_kmeans_spcontig_{mapping_tool}_{bin}.pdf", bbox_inches='tight', format="pdf")  
            plt.savefig(f"{outdir}/python_figs/elbow_UMAP_kmeans_spcontig_{mapping_tool}_{bin}.png", bbox_inches='tight', format="png")  
            plt.show()
            knee = round(kneedle.knee, 3)      
            
            # instantiate the k-means class, using optimal number of clusters
            # kmeans = KMeans( n_clusters=knee, init="random", n_init=10, random_state=1)
            kmeans = KMeans(n_clusters = knee, init = 'k-means++', random_state = 1)

            #fit k-means algorithm to data
            kmeans.fit(scaled_df)

            #view cluster assignments for each observation
            kmeans.labels_
            # https://www.machinelearningplus.com/predictive-modeling/k-means-clustering/
            y_pred = kmeans.fit_predict(scaled_df)
            df_merge['labels'] = kmeans.labels_
            
            # UMAP
            df = df.set_index('contig')
            embedding = umap.UMAP(random_state=42, metric="cosine").fit_transform(df.values)
            embedding.shape
                      
            fig, ax = plt.subplots(figsize=(7,5))
            g = sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], 
                        c=df_merge['arrcpm'],cmap="Reds",
                        edgecolor="black", style=y_pred, s=100,
                        legend="full", ax=ax)
            g.legend(title="Cluster",loc='upper right', 
                     bbox_to_anchor=(1.2, 1.0), ncol=1)

            norm = plt.Normalize(df_merge['arrcpm'].min(), df_merge['arrcpm'].max())
            sm = plt.cm.ScalarMappable(cmap="Reds", norm=norm)
            sm.set_array([])

            cax = fig.add_axes([ax.get_position().x1+0.05, ax.get_position().y0, 0.02, ax.get_position().height / 2])
            ax.figure.colorbar(sm, cax=cax, label="Abundance")

            ax.set_title(f'UMAP - {bin}', fontsize=20)
            # ax.suptitle("Abundance Profiles of Spacer-containg Contigs", fontsize = 15)
            ax.set_xlabel('UMAP 1', fontsize=15)
            ax.set_ylabel('UMAP 2', fontsize=15)
            fig.savefig(f"{outdir}/python_figs/UMAP_kmeans_spcontig_{mapping_tool}_{bin}.pdf", bbox_inches='tight', format="pdf")  
            fig.savefig(f"{outdir}/python_figs/UMAP_kmeans_spcontig_{mapping_tool}_{bin}.png", bbox_inches='tight', format="png")  
                        
        except Exception as e:
            print("error message: ",e)





# %%
