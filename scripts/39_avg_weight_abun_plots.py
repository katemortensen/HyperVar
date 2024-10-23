#!/usr/bin/env python
# coding: utf-8

# 

# In[1]:


import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
# import seaborn as sns
import argparse
# import random
import math 
import os


# In[ ]:


# Initialize parser
parser = argparse.ArgumentParser()
 
# Adding optional argument
parser.add_argument("-accessions", "--accessions", help = "path to file containing accession numbers (ex: \'/usr/home/wkdir/accessions.txt\' )", required=True)
parser.add_argument("-topmost_bin_count", "--topmost_bin_count", help = "quantity of topmost spacer abundant bins for comparison", required=True)
parser.add_argument("-repeat_guide", "--repeat_guide", help = "path to file containing accession numbers (ex: \'/usr/home/wkdir/accessions.txt\' )", required=True)
parser.add_argument("-abund_dir", "--abund_dir", help = "path to directory where CONTIG_RELATIVE_ABUNDANCE output exists (ex: \'/usr/home/wkdir/CONTIG_RELATIVE_ABUNDANCE\')", required=True)
parser.add_argument("-cdhit_dir", "--cdhit_dir", help = "path to directory where CONTIG_RELATIVE_ABUNDANCE output exists (ex: \'/usr/home/wkdir/CD-HIT-EST\')", required=True)
parser.add_argument("-outdir", "--outdir", help = "path to results directory", required=True)
parser.add_argument("-mapping_tool", "--mapping_tool", help = "(salmon or minimap)", required = True )

args = parser.parse_args()

accessions = args.accessions
topmost_bin_count = int(args.topmost_bin_count)
repeat_guide = args.repeat_guide
abund_dir = args.abund_dir
cdhit_dir = args.cdhit_dir
outdir = args.outdir
mapping_tool = args.mapping_tool
spcount_file = f'{outdir}/assembly2sp_count_sorted.tsv'

print(f'Comparing subset of {topmost_bin_count} most spacer abundant bins.')
print("output directory: ", outdir )

try: 
    os.mkdir(outdir) 
except OSError as error: 
    print(error)  

# stop and comment out hardcoded params

# In[2]:

# # paths
# mapping_tool = 'minimap'
# mapping_tool = 'salmon'

# accessions = '/home/kmorten/ADMB/SRR_list.txt'
# topmost_bin_count = 5 
# repeat_guide = 'rNA'
# abund_dir = '/home/kmorten/ADMB/METAWRAP_MAGS/CONTIG_RELATIVE_ABUNDANCE'
# cdhit_dir = '/home/kmorten/ADMB/METAWRAP_MAGS/CD-HIT-EST'
# outdir = '/home/kmorten/ADMB/METAWRAP_MAGS/results'
# spcount_file = '/home/kmorten/ADMB/METAWRAP_MAGS/results/assembly2sp_count_sorted.tsv'

# accessions = '/home/kmorten/ADMB/SRR_list.txt'
# topmost_bin_count = 5 
# repeat_guide = 'rNA'
# abund_dir = '/home/kmorten/ADMB/DASTOOL_MAGS/CONTIG_RELATIVE_ABUNDANCE'
# cdhit_dir = '/home/kmorten/ADMB/DASTOOL_MAGS/CD-HIT-EST'
# outdir = '/home/kmorten/ADMB/DASTOOL_MAGS/results'
# spcount_file = '/home/kmorten/ADMB/DASTOOL_MAGS/results/assembly2sp_count_sorted.tsv'

# accessions = '/home/kmorten/ADMB/SRR_list.txt'
# topmost_bin_count = 5 
# repeat_guide = 'rNA'
# abund_dir = '/home/kmorten/ADMB/MAGSCOT_MAGS/CONTIG_RELATIVE_ABUNDANCE'
# cdhit_dir = '/home/kmorten/ADMB/MAGSCOT_MAGS/CD-HIT-EST'
# outdir = '/home/kmorten/ADMB/MAGSCOT_MAGS/results'
# spcount_file = '/home/kmorten/ADMB/MAGSCOT_MAGS/results/assembly2sp_count_sorted.tsv'


# In[3]:

df = pd.read_csv(f"{abund_dir}/relative_contig_abundance_profiles_crt_{mapping_tool}.tsv" , sep="\t")
# df_dist = pd.read_csv(f"{abund_dir}/A50_cosine_distance_{mapping_tool}.tsv" , sep="\t")
df_A50 = pd.read_csv(f"{abund_dir}/A50_concensus_profile_{mapping_tool}.tsv" , sep="\t")
df_sp = pd.read_csv(f"{cdhit_dir}/sp2assembly_{repeat_guide}.out", sep="\t")
df_abund =  pd.read_csv(f"{abund_dir}/contig_abund_{mapping_tool}.tsv" , sep="\t")

# %% 
'''sample accessions'''
acc = open(accessions, "r").read().split('\n')
acc = acc.remove('')
 
'''topmost abundant bins'''
bin_sub = []
count = 0
with open(spcount_file, 'r') as fp:
    next(fp)
    for line in fp:
        linep = line.strip('\n').split('\t')
        assembly = linep[0]
        if 'bin' in assembly and assembly != 'unbinned' and len(bin_sub) <= topmost_bin_count:
            bin_sub.append(assembly)


# In[4]:
'''Distribution of Contig Approximated Spacer Abundances (Coassembly v. In-common)'''

# distribution of co-assembly (only) and in-common (co-assembly & bins) 
# abundances (y = freq, x = abund)

# Duplicate Contigs (contigs with many spacers)
# A contig can have multiple spacers
# Spacers from a contig can appear in different clusters 
# df_sp[df_sp[['assembly','contig','clust_diversity']].duplicated()].head()

# contigs can have multiple spacers... some contigs will be counted more than once depending on cluster diversity
df_clustdiv = df_sp[['assembly','contig','clust_diversity']].drop_duplicates() # contig only counted 1x for a particular clust_diversity
df_spcontig_abund = df_clustdiv.merge(df_abund, left_on=['assembly','contig'], right_on=['assembly','contig'],how="left")

# Note: If coassembly and individual assemblies have not been mapped by raw reads and thus don't have an abundance profile...
assert len(df_spcontig_abund[df_spcontig_abund.isna().any(axis=1)]) == 0
# assert len(df_spcontig_abund.dropna()) == len(df_spcontig_abund) - len(df_spcontig_abund[df_spcontig_abund.isna().any(axis=1)])
# df_spcontig_abund = df_spcontig_abund.dropna()

# Note: there are redundant nodes in the df_sp because some contigs have multiple spacers which can experience different sp cluster diversities than other spacers in the same contig
assert len(df_spcontig_abund) == len(df_spcontig_abund.drop_duplicates())
df_spcontig_abund.head()


# In[5]:

coassembly_only = df_spcontig_abund[df_spcontig_abund['clust_diversity'] == 'C']
assert len(coassembly_only) == len(coassembly_only[['contig','arrcpm','awrrcpm']].drop_duplicates())
coassembly_only = np.log10(coassembly_only['arrcpm'])

individual_only = df_spcontig_abund[df_spcontig_abund['clust_diversity'] == 'I']
assert len(individual_only) == len(individual_only[['contig','arrcpm','awrrcpm']].drop_duplicates())
individual_only = np.log10(individual_only['arrcpm'])

incommon = df_spcontig_abund[df_spcontig_abund['clust_diversity'].isin(['CI','CIM'])]
incommon = incommon[['contig', 'arrcpm', 'awrrcpm']].drop_duplicates()
incommon = np.log10(incommon['arrcpm'])

lower_bound = math.log10(min(df_spcontig_abund['arrcpm']))
upper_bound = math.log10(max(df_spcontig_abund['arrcpm']))
bins = np.linspace(lower_bound, upper_bound, 100)

# plt.hist(coassembly_only, bins, alpha=0.5, label='sp contigs from exclv coassembly comprised clusters')
# plt.hist(incommon, bins, alpha=0.5, label='sp contigs from individial & co-assembly comprised clusters')
plt.hist(coassembly_only, bins, alpha=0.5, label= r'co-assembly $\oplus$' ) # co-assembly exclusively
plt.hist(individual_only, bins, alpha=0.5, label= r'individual $\oplus$' ) # individual exclusively
plt.hist(incommon, bins, alpha=0.5, label=r'individual $\cap$ co-assembly') # co-assembly/individual intersection

plt.legend(loc='upper right', title = "Cluster Composition Type")
plt.xlabel('Abundance as log10(arrcpm)')
plt.ylabel('Frequency')
#plt.title('Distribution of Spacer Contig Abundnaces')
plt.title('Distribution of Contig Abundance for Spacer Cluster Type')

try: 
    os.mkdir(f'{outdir}/python_figs') 
except OSError as error: 
    print(error)  

plt.savefig("%s/python_figs/dist_spcontig_abund_%s.pdf" % (outdir, mapping_tool), format="pdf")
plt.savefig("%s/python_figs/dist_spcontig_abund_%s.png" % (outdir, mapping_tool), format="png")  

plt.show()

############################## save df as tsv for rplot 
cluster_composition = ['coassembly_only']*len(coassembly_only) + ['incommon']*len(incommon) + ['individual']*len(individual_only)
log10_arrcpm = list(coassembly_only) + list(incommon) + list(individual_only)
assert len(cluster_composition) == len(log10_arrcpm)

df_out = pd.DataFrame(zip(cluster_composition, log10_arrcpm), columns = ['cluster_composition', 'log10_arrcpm'])
df_out.head()
df_out.to_csv(outdir + '/dist_spcontig_abund_%s.tsv' % (mapping_tool), sep="\t", index=False) 


# #### Distribution of Spacer Count in Contigs

# In[ ]:

# Duplicate Contigs (contigs with many spacers)
df_sp = pd.read_csv(f"{cdhit_dir}/sp2assembly_{repeat_guide}.out", sep="\t")

contig2spcount = dict(zip(list(set(df_sp['contig'])), [1]*len(list(set(df_sp['contig'])))))

for contig_name in list(set(df_sp['contig'])):
    df_tmp = df_sp[df_sp['contig'] == contig_name]
    spacers = len(list(set(df_tmp['sp_loci'])))
    contig2spcount.update({contig_name: spacers})

df_sp['raw_sp_count'] = [contig2spcount.get(contig_name) for contig_name in df_sp['contig']]
df_sp.head()


# In[ ]:

# distribution of spacers in contig

coassembly_only = df_sp[df_sp['clust_diversity'] == 'C']
coassembly_only = coassembly_only[['contig','clust_diversity','raw_sp_count']].drop_duplicates()
coassembly_only = coassembly_only['raw_sp_count']

individual_only = df_sp[df_sp['clust_diversity'] == 'I']
individual_only = individual_only[['contig','clust_diversity','raw_sp_count']].drop_duplicates()
individual_only = individual_only['raw_sp_count']

incommon = df_sp[df_sp['clust_diversity'].isin(['CI','CIM'])]
incommon = incommon[['contig', 'raw_sp_count']].drop_duplicates()
incommon = incommon['raw_sp_count']

lower_bound = min(df_sp['raw_sp_count'])
upper_bound = max(df_sp['raw_sp_count'])
bins = np.linspace(lower_bound, upper_bound, 100)

# plt.hist(coassembly_only, bins, alpha=0.5, label='sp contigs from exclv coassembly comprised clusters')
# plt.hist(incommon, bins, alpha=0.5, label='sp contigs from individial & co-assembly comprised clusters')
plt.hist(coassembly_only, bins, alpha=0.5, label=r'co-assembly $\oplus$')
plt.hist(individual_only, bins, alpha=0.5, label=r'individual $\oplus$')
plt.hist(incommon, bins, alpha=0.5, label=r'individual $\cap$ co-assembly')
plt.legend(loc='upper right')
plt.xlabel('Spacer Count in Contig')
plt.ylabel('Frequency')
#plt.title('Histogram of Spacer Contig Abundnaces')
plt.title('Distribution of Spacer Counts in Contig')

try: 
    os.mkdir(f'{outdir}/python_figs') 
except OSError as error: 
    print(error)  

plt.savefig(f"{outdir}/python_figs/dist_spcount_{mapping_tool}.pdf", format="pdf")  
plt.savefig(f"{outdir}/python_figs/dist_spcount_{mapping_tool}.png", format="png") 
plt.show()


# #### Distribution of Average Weighted Contig Abundances
# 

# In[ ]:
df_arrcpm = pd.read_csv(f"{abund_dir}/contig_abund_{mapping_tool}.tsv" , sep="\t")
df_arrcpm = df_arrcpm[df_arrcpm['N50_contig_boolean']==1]

coassembly = df_arrcpm[df_arrcpm['assembly'] == 'coassembly']['arrcpm'] 
pattern = '|'.join(['ERR', 'SRR'])
indiv = df_arrcpm[df_arrcpm['assembly'].str.contains(pattern, na=False)]['arrcpm'] 
mag = df_arrcpm[df_arrcpm['assembly'].str.contains('bin')]['arrcpm']
assert len(coassembly) + len(indiv) + len(mag) == len(df_arrcpm)

################### transform data
scale_factor = 10**8
coassembly = np.log10(coassembly.apply(lambda x: x*scale_factor))
indiv = np.log10(indiv.apply(lambda x: x*scale_factor))
mag = np.log10(mag.apply(lambda x: x*scale_factor))

lower_bound = min(np.log10(df_arrcpm['arrcpm'].apply(lambda x: x*scale_factor)))
upper_bound = max(np.log10(df_arrcpm['arrcpm'].apply(lambda x: x*scale_factor)))
##################

print("lower_bound:", lower_bound)
print("upper_bound:", upper_bound)
bins = np.linspace(lower_bound, upper_bound, 50)

plt.hist(coassembly, bins, alpha=0.5, label='Co-assembly')
plt.hist(indiv, bins, alpha=0.5, label='Individual')
plt.hist(mag, bins, alpha=0.5, label='MAG')
plt.legend(loc='upper right')
plt.xlabel('Abundance log10(arrcpm*10^8)')
plt.ylabel('Frequency')
plt.title('Distribution of N50 Contig Abundance')

try: 
    os.mkdir(f'{outdir}/python_figs') 
except OSError as error: 
    print(error)  

plt.savefig(f"{outdir}/python_figs/dist_N50contig_abund_{mapping_tool}.pdf", format="pdf")  
plt.savefig(f"{outdir}/python_figs/dist_N50contig_abund_{mapping_tool}.png", format="png")  
plt.show()

############################## save df as tsv for rplot 
assembly_type = ['coassembly']*len(coassembly) + ['individual']*len(indiv) + ['mag']*len(mag)
log10_arrcpm_10e8 = list(coassembly) + list(indiv) + list(mag)
assert len(assembly_type) == len(log10_arrcpm_10e8)

df_out = pd.DataFrame(zip(assembly_type, log10_arrcpm_10e8), columns = ['assembly_type', 'log10_arrcpm_10e8'])
df_out.to_csv(outdir + '/dist_N50contig_abund_%s.tsv' % (mapping_tool), sep="\t", index=False) 

# %%
