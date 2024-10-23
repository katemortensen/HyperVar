#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Author: Kate Mortensen
# Date: 1/17/2023
Purpose = "Purpose: to calculate pearson correlation coefficients between relative abundance profiles of contigs."


# In[16]:
import argparse
import pandas as pd 
import os
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt


# %%
# Initialize parser
parser = argparse.ArgumentParser()
 
# Adding optional argument
parser.add_argument("-topmost_bin_count", "--topmost_bin_count", help = "quantity of topmost spacer abundant bins for comparison", required=True)
parser.add_argument("-abund_dir", "--abund_dir", help = "path to directory where CONTIG_RELATIVE_ABUNDANCE output exists (ex: \'/usr/home/wkdir/CONTIG_RELATIVE_ABUNDANCE\')", required=True)
parser.add_argument("-outdir", "--outdir", help = "path to results directory", required=True)
parser.add_argument("-mapping_tool", "--mapping_tool", help = "(salmon or minimap)", required = True )

args = parser.parse_args()

topmost_bin_count = int(args.topmost_bin_count)
abund_dir = args.abund_dir
outdir = args.outdir
mapping_tool = args.mapping_tool

try: 
    os.mkdir(outdir) 
except OSError as error: 
    print(error) 


print(f'Comparing subset of {topmost_bin_count} most spacer abundant bins')
print("output directory: ", outdir )

# stop and comment out hardcoded params

# # %%
# mapping_tool = 'minimap'
# topmost_bin_count = 5 

# abund_dir = '/home/kmorten/ADMB/DASTOOL_MAGS/CONTIG_RELATIVE_ABUNDANCE'
# outdir = '/home/kmorten/ADMB/DASTOOL_MAGS/results'
# spcount_file = '/home/kmorten/ADMB/DASTOOL_MAGS/results/assembly2sp_count_sorted.tsv'

# abund_dir = '/home/kmorten/ADMB/METAWRAP_MAGS/CONTIG_RELATIVE_ABUNDANCE'
# outdir = '/home/kmorten/ADMB/METAWRAP_MAGS/results'
# spcount_file = '/home/kmorten/ADMB/METAWRAP_MAGS/results/assembly2sp_count_sorted.tsv'


# ## Pearson Correlation Coefficients of Contig Relative Abundance Profiles 

# In[17]:

df_abund = pd.read_csv(f"{abund_dir}/relative_contig_abundance_profiles_crt_{mapping_tool}.tsv" , sep="\t")
df_A50prof = pd.read_csv(f"{abund_dir}/A50_concensus_profile_{mapping_tool}.tsv" , sep="\t")
df_A50cos = pd.read_csv(f"{abund_dir}/A50_cosine_distance_{mapping_tool}.tsv" , sep="\t")

# In[19]:

'''topmost abundant bins'''
bin_sub = []
count = 0
spcount_file = f'{outdir}/assembly2sp_count_sorted.tsv'
with open(spcount_file, 'r') as fp:
    next(fp)
    for line in fp:
        linep = line.strip('\n').split('\t')
        assembly = linep[0]
        if 'bin' in assembly and assembly != 'unbinned' and len(bin_sub) < topmost_bin_count:
            bin_sub.append(assembly)


'''src_dict'''
'''sample->read_count'''
'''src (sample read counts) = total read counts for each sample '''
src_path = f'{abund_dir}/sample_reads_count.txt'
src_df = pd.read_csv( src_path , sep="\t")
src_dict = dict(zip(src_df['sample'], src_df['total_reads']))
samples = list(src_df['sample'])


# In[20]:


def pca_pearson(df_abund, df_A50prof, samples, assembly):

    # subsetting & concatenating A50 profile with abundance 
    df_A50prof_tmp = df_A50prof[df_A50prof['assembly'] == assembly]
    df_A50prof_tmp.insert(1,'contig', ['concensus'])
    df_A50prof_tmp.insert(2,'crt', [0])

    cols = ['assembly', 'contig', 'crt'] + samples
    df_abund_tmp = df_abund[df_abund['assembly'] == assembly][cols]
    df = pd.concat([df_A50prof_tmp, df_abund_tmp], ignore_index=True, sort=False)
    
    assert len(df) > 0

    # setting contig column as index column
    df.set_index("contig", inplace = True)

    # swap rows and columns for pca
    dfT = df[samples].T

    # pearson correlation between contig RRCPM  profiles with respect to bin
    df_corr = dfT.corr(method='pearson')

    # https://towardsdatascience.com/pca-using-python-scikit-learn-e653f8989e60
    # Standardizing the features
    scaled_df = StandardScaler().fit_transform(df_corr)

    pca = PCA(n_components=10)
    principalComponents = pca.fit_transform(scaled_df)
    principalDf = pd.DataFrame(data = principalComponents, 
        columns = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10'])

    assert len(principalDf) == len(df)

    # Inserting the column at the beginning in the DataFrame
    finalDf = principalDf.copy()
    finalDf.insert(loc = 0, column = 'crt', value = list(df['crt']))
    finalDf.insert(loc = 0, column = 'contig', value = df.index)
    finalDf.insert(loc = 0, column = 'assembly', value = list(df['assembly']))

    # labels (+/-) crispr
    finalDf.loc[finalDf["crt"] > 0, "crt_sign"] = "+"
    finalDf.loc[finalDf["crt"] == 0, "crt_sign"] = "-"
    finalDf.loc[finalDf["crt"] == 'NA', "crt_sign"] = "-"
    
    return finalDf


# ## Plot PCA of Contig RRCPM Pearson Correlations 

# In[21]:


def pca_plot(finalDf, df_A50cos, outdir, assembly):

    # contigs without crispr (including outliers)
    CRTabsent = finalDf[finalDf['crt_sign'] == '-']
    # contigs without crispr (excluding outliers)
    CRTabsent_outliers_removed = CRTabsent.merge(df_A50cos[df_A50cos['outlier']==0][['assembly','contig']], on = ['assembly','contig'], how = 'inner')
    CRTabsent.reset_index(inplace=True)
    # N50 contigs without crispr
    CRTabsent_N50contigs = CRTabsent.merge(df_A50cos[df_A50cos['N50_contig_boolean']==1][['assembly','contig']], on = ['assembly','contig'], how="inner")

    # contigs with crispr (including outliers)
    CRTpresent = finalDf[finalDf['crt_sign'] == '+']
    # contigs with crispr (excluding outliers)
    CRTpresent_outliers_removed = CRTpresent.merge(df_A50cos[df_A50cos['outlier']==0][['assembly','contig']], on = ['assembly','contig'], how = 'inner')
    CRTpresent.reset_index(inplace=True)
    # N50 contigs with crispr
    CRTpresent_N50contigs = CRTpresent.merge(df_A50cos[df_A50cos['N50_contig_boolean']==1][['assembly','contig']], on = ['assembly','contig'], how="inner")

    # concensus (A50 profile)
    A50prof = finalDf[finalDf['contig']=='concensus']

    # N50 contigs
    N50contigs = finalDf.merge(df_A50cos[df_A50cos['N50_contig_boolean']==1][['assembly','contig']], on = ['assembly','contig'], how="inner")
    
    print("##########################################")
    print("Assembly: ", assembly)
    print("##########################################")
    print("contigs without crispr (including outliers)", len(CRTabsent))
    print("contigs without crispr (excluding outliers)", len(CRTabsent_outliers_removed))
    print("delta: ", (len(CRTabsent_outliers_removed)-len(CRTabsent))*100/len(CRTabsent), "prcnt change")
    print("##########################################")
    print("contigs with crispr (including outliers)", len(CRTpresent))
    print("contigs with crispr (excluding outliers)", len(CRTpresent_outliers_removed))
    print("delta: ", (len(CRTpresent_outliers_removed)-len(CRTpresent))*100/len(CRTpresent), "prcnt change")
    
    info_type = ['assembly',
                 'contigs_without_crispr_including_outliers',
                 'contigs_without_crispr_excluding_outliers',
                 'contigs_with_crispr_including_outliers',
                 'contigs_with_crispr_excluding_outliers']
    info_rslt = [assembly,
                       len(CRTabsent),
                       len(CRTabsent_outliers_removed),
                       len(CRTpresent),
                       len(CRTpresent_outliers_removed)]
    contig_counts_dict = dict(zip(info_type,info_rslt))

    plt.figure(figsize=(10,7.5))
    scatter = plt.scatter(CRTabsent["PC1"], CRTabsent["PC2"], c="lightsteelblue")
    #scatter = plt.scatter(CRTpresent_outliers_removed["PC1"], CRTpresent_outliers_removed["PC2"], c="orange")
    scatter = plt.scatter(CRTpresent["PC1"], CRTpresent["PC2"], c="orange")
    scatter = plt.scatter(A50prof['PC1'], A50prof['PC2'], marker="d", edgecolors= 'k', linewidths=1.5, s = 100, color='red')
    
    #scatter = plt.scatter(CRTabsent_N50contigs['PC1'], CRTabsent_N50contigs['PC2'], edgecolors= 'k', linewidths=1.5, color='orange')
    #scatter = plt.scatter(CRTpresent_N50contigs['PC1'], CRTpresent_N50contigs['PC2'], edgecolors= 'k', linewidths=1.5)
    #scatter = plt.scatter(N50contigs['PC1'], N50contigs['PC2'], c = "green")
    
    
    plt.title('Contig RRCPM Pearson Correlations for %s' % (assembly), fontsize=20, loc='left') 
    
    # plt.legend(["-crt, +outliers", "+crt, -outliers","consensus"], fontsize=15, 
    # title = "Type", title_fontsize = 15, 
    # bbox_to_anchor=(1.02, 1), loc='upper left')
    
    plt.legend(["-crt", "+crt","consensus"], fontsize=15, 
    title = "Type", title_fontsize = 15, 
    bbox_to_anchor=(1.02, 1), loc='upper left')

    plt.xlabel('PC1', fontsize=20)
    plt.ylabel('PC2', fontsize=20)
    
    try: 
        os.mkdir(f'{outdir}/python_figs') 
    except OSError as error: 
        print(error) 
        
    plt.savefig(f"{outdir}/python_figs/{assembly}_rrcpm_pearson_{mapping_tool}.pdf", bbox_inches='tight', format="pdf")  
    plt.savefig(f"{outdir}/python_figs/{assembly}_rrcpm_pearson_{mapping_tool}.png", bbox_inches='tight', format="png")  
    plt.show()

    ############################## output df
    # rrcpm_pearson_corr_pc1 = list(CRTabsent['PC1']) + list(CRTpresent_outliers_removed['PC1']) + list(A50prof['PC1'])
    # rrcpm_pearson_corr_pc2 = list(CRTabsent['PC2']) + list(CRTpresent_outliers_removed['PC2']) + list(A50prof['PC2'])
    # type = ['CRTabsent']*len(CRTabsent) + ['CRTpresent_outliers_removed']*len(CRTpresent_outliers_removed) + ['concensus']*len(A50prof)
    
    rrcpm_pearson_corr_pc1 = list(CRTabsent['PC1']) + list(CRTpresent['PC1']) + list(A50prof['PC1'])
    rrcpm_pearson_corr_pc2 = list(CRTabsent['PC2']) + list(CRTpresent['PC2']) + list(A50prof['PC2'])
    type = ['CRTabsent']*len(CRTabsent) + ['CRTpresent']*len(CRTpresent) + ['concensus']*len(A50prof)
    
    assembly_ls = [assembly]*len(type)

    assert len(assembly_ls) == len(rrcpm_pearson_corr_pc1) == len(rrcpm_pearson_corr_pc2) == len(type)

    df_out = pd.DataFrame(zip(assembly_ls,rrcpm_pearson_corr_pc1,rrcpm_pearson_corr_pc2,type),
     columns = ['assembly', 'rrcpm_pearson_corr_pc1', 'rrcpm_pearson_corr_pc2', 'type'])

    return df_out, contig_counts_dict


# In[22]:


file_out1 = f'{outdir}/assembly_rrcpm_pearson_pca_{mapping_tool}.tsv'
file_out2 = f'{outdir}/assembly_rrcpm_pearson_pca_contig_counts_{mapping_tool}.tsv'

count = 0
for assembly in bin_sub:
    finalDf =  pca_pearson(df_abund,df_A50prof, samples, assembly)
    result = pca_plot(finalDf,df_A50cos, outdir, assembly)
    df_out = result[0]
    contig_counts_dict = result[1]
    new_row = [str(contig_counts_dict.get(i)) for i in contig_counts_dict.keys()]
    new_row_string = "\t".join(new_row).strip("\t") + "\n"
    
    if count == 0:
        df_out.to_csv(file_out1, mode='w', header=True, sep='\t', index=False)
        contig_counts_header = "\t".join(list(contig_counts_dict.keys())).strip("\t") + "\n"
        with open(file_out2, 'w') as f2:
            f2.write(contig_counts_header)
            f2.write(new_row_string) 
        count += 1
    else :
        df_out.to_csv(file_out1, mode='a', header=False, sep='\t', index=False)
        with open(file_out2, 'a') as f2:
            f2.write(new_row_string) 




# %%
