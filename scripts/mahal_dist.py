

# %%   
import argparse
import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from kneed import DataGenerator, KneeLocator
from sklearn.decomposition import PCA
import scipy as stats 
from scipy.stats import chi2 
import umap.umap_ as umap
import matplotlib.cm as cm


# %% 

mag_tool = 'METAWRAP'
# mag_tool = 'DASTOOL'
# # mag_tool = 'MAGSCOT'

mapping_tool = 'minimap'
# mapping_tool = 'salmon'

topmost_bin_count = 5
accessions = "/home/kmorten/ADMB/SRR_list.txt" 
spacer_seq_file = f'/home/kmorten/ADMB/{mag_tool}_MAGS/CD-HIT-EST/concat_rNA_spacer.seq'
spacer_counts_file = f'/home/kmorten/ADMB/{mag_tool}_MAGS/results/assembly2sp_count_sorted.tsv'
abund_file = f'/home/kmorten/ADMB/{mag_tool}_MAGS/CONTIG_RELATIVE_ABUNDANCE/relative_contig_abundance_profiles_crt_{mapping_tool}.tsv'
avg_abund_file = f'/home/kmorten/ADMB/{mag_tool}_MAGS/CONTIG_RELATIVE_ABUNDANCE/contig_abund_{mapping_tool}.tsv' 
outdir = f'/home/kmorten/ADMB/{mag_tool}_MAGS/results/'

# wabund_file = f'/home/kmorten/ADMB/{mag_tool}_MAGS/CONTIG_RELATIVE_ABUNDANCE/weighted_relative_contig_abundance_profiles_crt_{mapping_tool}.tsv'
# A50_file = f'/home/kmorten/ADMB/{mag_tool}_MAGS/CONTIG_RELATIVE_ABUNDANCE/A50_cosine_distance_{mapping_tool}.tsv'

# %% 
samples = open(accessions, "r").read().strip('\n').split('\n')

df_spcount = pd.read_csv(f'{spacer_counts_file}', sep = '\t')
df_spcount_sub = df_spcount.loc[(df_spcount['assembly_type']=='Co-assembled') & (df_spcount['assembly'] != 'coassembly')]
df_abund = pd.read_csv(f'{abund_file}', sep = '\t')

# df_wabund = pd.read_csv(f'{wabund_file}', sep = '\t')
# df_A50 = pd.read_csv(f'{A50_file}', sep = '\t')


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


# Mahalanobis function
def calculateMahalanobis(y=None, data=None, cov=None): 
  
    y_mu = y - np.mean(data) 
    if not cov: 
        cov = np.cov(data.values.T) 
    inv_covmat = np.linalg.inv(cov) 
    left = np.dot(y_mu, inv_covmat) 
    mahal = np.dot(left, y_mu.T) 
    return mahal.diagonal()

# %% 


for bin in top_bins:
# for bin in [top_bins[0]]:
    
    spcontig_matrix = f'{outdir}/spcontig_matrix_{mapping_tool}_{bin}.tsv'    
    if os.path.exists(spcontig_matrix) == False:
        print('matrix file not found at path:', spcontig_matrix)    
    df_spcontig = pd.read_csv(f'{spcontig_matrix}', sep = '\t')
    assert len(df_spcontig['spacer_contigs']) == len(set(df_spcontig['spacer_contigs']))
    
    # add abundance profiles for all contigs in bin
    df_abund_tmp = df_abund[df_abund['assembly']==bin].copy()
    
    # add average abundance for all contigs in bin
    df_avg_abund = pd.read_csv(f'{avg_abund_file}', sep = '\t')
    df_abund_tmp = df_abund_tmp.merge(df_avg_abund[['contig', 'arrcpm']], on = 'contig', how = 'left') 
    
    # add spacer presence
    df_abund_tmp['spacer_presence'] = len(df_abund_tmp)*[0]
    for contig in df_spcontig['spacer_contigs']:
        df_abund_tmp.loc[df_abund_tmp['contig']==contig,'spacer_presence']=1
    assert len(df_abund_tmp[df_abund_tmp['spacer_presence']==1]) == len(df_spcontig)
        

    if len(df_abund_tmp['contig']) <= 1:
        print(f'mag ({bin}) needs at least 2 contigs with spacers.')
    else :
        try :
            # df_mahal = df_abund_tmp.drop(["assembly_type", "assembly", "contig","crt","spacer_presence"], axis=1).copy()
            df_mahal = df_abund_tmp.copy()
            # Creating a new column in the dataframe that holds 
            # the Mahalanobis distance for each row 
            # df_mahal['Mahalanobis'] = calculateMahalanobis(y=df_mahal, data=df_mahal[list(df_mahal.columns)]) 
            df_mahal['Mahalanobis'] = calculateMahalanobis(y=df_mahal[samples], data=df_mahal[samples]) 

            # calculate p-value for each mahalanobis distance 
            df_mahal['p-val'] = 1 - chi2.cdf(df_mahal['Mahalanobis'], len(samples)-1) 
            
            # add contig and spacer_presence columns back
            df_mahal['spacer_presence']= df_abund_tmp['spacer_presence'] 

            # add outlier boolean 
            df_mahal['outlier'] = len(df_mahal)*[0]
            for contig in df_mahal['contig']:
                df_mahal.loc[df_mahal['p-val'] <= 0.001,'outlier']=1                        

            outliers_count = 0
            for i in df_mahal['p-val']:
                if i <= 0.001:
                    outliers_count += 1
            percent_outliers = round((outliers_count/len(df_mahal))*100, 2)

            spcontig_outliers_count = 0
            for i in df_mahal[df_mahal['spacer_presence']==1]['p-val']:
                if i <= 0.001:
                    spcontig_outliers_count += 1
            total_spcontigs = len(df_mahal[df_mahal['spacer_presence']==1])
            sp_percent_outliers = round((spcontig_outliers_count/total_spcontigs)*100, 2)
           
            nonspcontig_outliers_count = 0
            for i in df_mahal[df_mahal['spacer_presence']==0]['p-val']:
                if i <= 0.001:
                    nonspcontig_outliers_count += 1  
            total_nonspcontigs = len(df_mahal[df_mahal['spacer_presence']==0])
            nonsp_percent_outliers = round((nonspcontig_outliers_count/total_nonspcontigs)*100, 2)
           
            print(bin)
            print('total contigs: ', len(df_mahal))
            print('raw percent outliers: ', percent_outliers, '%')
            print('percent spcontigs: ', round((total_spcontigs/len(df_mahal))*100, 2), '%')
            print('percentage of outliers that are spcontigs: ', round((spcontig_outliers_count/outliers_count)*100, 2), '%')
            print('percent of outliers that are nonspcontigs: ',round((nonspcontig_outliers_count/outliers_count)*100, 2), '%')
            print('percent of spcontigs that are outliers: ', sp_percent_outliers, '%')
            print('percent of non_spcontigs that are outliers: ', nonsp_percent_outliers, '%', '\n')

            
            '''K-means'''
            
            scaled_df = StandardScaler().fit_transform(df_mahal.drop(["assembly_type", "assembly", "contig","crt","spacer_presence"], axis=1))

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
            # kneedle.plot_knee(figsize=(7,5))
            # plt.title(f'Knee Point - {bin}', fontsize = 20)
            # plt.xlabel("Number of Clusters", fontsize = 15)
            # plt.ylabel("SSE", fontsize = 15)
            # plt.legend(["data", "knee/elbow"], fontsize="15", loc ="upper right")
            # # plt.savefig(f"{outdir}/python_figs/elbow_UMAP_kmeans_{mapping_tool}_{bin}.pdf", bbox_inches='tight', format="pdf")  
            # # plt.savefig(f"{outdir}/python_figs/elbow_UMAP_kmeans_{mapping_tool}_{bin}.png", bbox_inches='tight', format="png")  
            # # plt.show()
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
            df_mahal['labels'] = kmeans.labels_

            '''UMAP'''            
            df_mahal = df_mahal.set_index('contig')
            embedding = umap.UMAP(random_state=42, metric='cosine').fit_transform(df_mahal[samples])
            embedding.shape

            '''PLOT 1, COLOR = ABUNDANCE, SHAPE = SPACER PRESENCE'''
            fig, ax = plt.subplots(figsize=(7,5))
            g = sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], 
                        c=df_mahal['arrcpm'],cmap="Reds",
                        edgecolor="black", style=df_mahal['spacer_presence'], s=100,
                        legend="full", ax=ax)
            g.legend(title="Spacer Presence",loc='upper right',bbox_to_anchor=(1.2, 1.0), ncol=1)

            norm = plt.Normalize(df_mahal['arrcpm'].min(), df_mahal['arrcpm'].max())
            sm = plt.cm.ScalarMappable(cmap="Reds", norm=norm)
            sm.set_array([])

            cax = fig.add_axes([ax.get_position().x1+0.05, ax.get_position().y0, 0.02, ax.get_position().height / 2])
            ax.figure.colorbar(sm, cax=cax, label="Abundance")

            ax.set_title(f'UMAP - {bin}', fontsize=20)
            # ax.suptitle("Abundance Profiles of Contigs", fontsize = 15)
            ax.set_xlabel('UMAP 1', fontsize=15)
            ax.set_ylabel('UMAP 2', fontsize=15)
            # fig.savefig(f"{outdir}/python_figs/UMAP_kmeans_abund_{mapping_tool}_{bin}.pdf", bbox_inches='tight', format="pdf")  
            # fig.savefig(f"{outdir}/python_figs/UMAP_kmeans_abund_{mapping_tool}_{bin}.png", bbox_inches='tight', format="png")  
            plt.show()
            
            '''PLOT 2, COLOR = PVAL, SHAPE = SPACER PRESENCE'''
            fig, ax = plt.subplots(figsize=(7,5))
            g = sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], 
                        c=df_mahal['p-val'],cmap="Blues_r",
                        edgecolor="black", style=df_mahal['spacer_presence'], s=100,
                        legend="full", ax=ax)
            g.legend(title="Spacer Presence",loc='upper right',bbox_to_anchor=(1.2, 1.0), ncol=1)

            norm = plt.Normalize(df_mahal['p-val'].min(), df_mahal['p-val'].max())
            sm = plt.cm.ScalarMappable(cmap="Blues_r", norm=norm)
            sm.set_array([])

            cax = fig.add_axes([ax.get_position().x1+0.05, ax.get_position().y0, 0.02, ax.get_position().height / 2])
            ax.figure.colorbar(sm, cax=cax, label="p-val")

            ax.set_title(f'UMAP Cos - {bin}', fontsize=20)
            # ax.suptitle("Abundance Profiles of Contigs", fontsize = 15)
            ax.set_xlabel('UMAP 1', fontsize=15)
            ax.set_ylabel('UMAP 2', fontsize=15)
            # fig.savefig(f"{outdir}/python_figs/UMAP_kmeans_pval_{mapping_tool}_{bin}.pdf", bbox_inches='tight', format="pdf")  
            # fig.savefig(f"{outdir}/python_figs/UMAP_kmeans_pval_{mapping_tool}_{bin}.png", bbox_inches='tight', format="png")  
            plt.show()
            
            
            '''PLOT 3, COLOR = MAHAL DISTANCE'''
            # fig, ax = plt.subplots(figsize=(7,5))
            # g = sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], 
            #             c=df_mahal['Mahalanobis'],cmap=cm.gray,
            #             edgecolor="black", style=y_pred, s=100,
            #             legend="full", ax=ax)
            # g.legend(title="Cluster",loc='upper right',bbox_to_anchor=(1.2, 1.0), ncol=1)

            # norm = plt.Normalize(df_mahal['Mahalanobis'].min(), df_mahal['Mahalanobis'].max())
            # sm = plt.cm.ScalarMappable(cmap="Reds", norm=norm)
            # sm.set_array([])

            # cax = fig.add_axes([ax.get_position().x1+0.05, ax.get_position().y0, 0.02, ax.get_position().height / 2])
            # ax.figure.colorbar(sm, cax=cax, label="Mahalanobis")

            # ax.set_title(f'UMAP - {bin}', fontsize=20)
            # ax.suptitle("Abundance Profiles of Contigs", fontsize = 15)
            # ax.set_xlabel('UMAP 1', fontsize=15)
            # ax.set_ylabel('UMAP 2', fontsize=15)
            # fig.savefig(f"{outdir}/python_figs/UMAP_kmeans_pval_{mapping_tool}_{bin}.pdf", bbox_inches='tight', format="pdf")  
            # fig.savefig(f"{outdir}/python_figs/UMAP_kmeans_pval_{mapping_tool}_{bin}.png", bbox_inches='tight', format="png")  
            
            '''PLOT 4, COLOR = OUTLIER, SHAPE=SPACER PRESENCE'''
            # fig, ax = plt.subplots(figsize=(7,5))
            # g = sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], 
            #             hue=df_mahal['outlier'],
            #             edgecolor="black", style=df_mahal['spacer_presence'],
            #             s=100, legend="full", ax=ax)
            # g.legend(title="Spacer Presence",loc='upper right',bbox_to_anchor=(1.2, 1.0), ncol=1)

            # ax.set_title(f'UMAP - {bin}', fontsize=20)
            # # ax.suptitle("Abundance Profiles of Contigs", fontsize = 15)
            # ax.set_xlabel('UMAP 1', fontsize=15)
            # ax.set_ylabel('UMAP 2', fontsize=15)
            # # fig.savefig(f"{outdir}/python_figs/UMAP_kmeans_mahaloutlier_{mapping_tool}_{bin}.pdf", bbox_inches='tight', format="pdf")  
            # # fig.savefig(f"{outdir}/python_figs/UMAP_kmeans_mahaloutlier_{mapping_tool}_{bin}.png", bbox_inches='tight', format="png")  
            # plt.show()
            
            '''PLOT 5, COLOR = SPACER PRESENCE'''
            # fig, ax = plt.subplots(figsize=(7,5))
            # g = sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], 
            #             hue=df_mahal['spacer_presence'],
            #             edgecolor="black", style=y_pred, s=100,
            #             legend="full", ax=ax)
            # g.legend(title="Cluster",loc='upper right',bbox_to_anchor=(1.2, 1.0), ncol=1)

            # ax.set_title(f'UMAP - {bin}', fontsize=20)
            # ax.suptitle("Abundance Profiles of Contigs", fontsize = 15)
            # ax.set_xlabel('UMAP 1', fontsize=15)
            # ax.set_ylabel('UMAP 2', fontsize=15)
            # # fig.savefig(f"{outdir}/python_figs/UMAP_kmeans_presence_{mapping_tool}_{bin}.pdf", bbox_inches='tight', format="pdf")  
            # # fig.savefig(f"{outdir}/python_figs/UMAP_kmeans_presence_{mapping_tool}_{bin}.png", bbox_inches='tight', format="png")  
                      
        except Exception as e:
            print("error message: ",e)

# %%
