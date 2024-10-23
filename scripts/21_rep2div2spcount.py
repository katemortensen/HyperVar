#!/usr/bin/env python
# coding: utf-8

# %% 

# Author: Kate Mortensen
# Date: 5-9-2023
# Purpose: data file to create stacked barplot of spacer counts with respect to repeat guide and assembly type 

# In[2]:

import pandas as pd
import numpy as np
import argparse
import os


# In[3]:

# Initialize parser
parser = argparse.ArgumentParser()
 
# Adding optional argument
parser.add_argument("-cdhit_dir", "--cdhit_dir", help = "path to directory where CD-HIT-EST output exists (ex: \'/usr/home/wkdir/CD-HIT-EST\')", required=True)
parser.add_argument("-repeat_guides", "--repeat_guides", help = "path to repeat guides fasta (/wkdir/CD-HIT-EST/top25_representative_repeats.fa)", required=True)
parser.add_argument("-outdir", "--outdir", help = "path to save output directory", required=True)

# Read arguments from command line
args = parser.parse_args()
cdhit_dir = args.cdhit_dir
repeat_guides = args.repeat_guides
outdir = args.outdir

try: 
    os.mkdir(outdir) 
except OSError as error: 
    print(error) 

# %%

# cdhit_dir = '/home/kmorten/ADMB/DASTOOL_MAGS/CD-HIT-EST'
# repeat_guides = "/home/kmorten/ADMB/DASTOOL_MAGS/CD-HIT-EST/top25_representative_repeats.fa"
# outdir = "/home/kmorten/ADMB/DASTOOL_MAGS/results"

# try: 
#     os.mkdir(outdir) 
# except OSError as error: 
#     print(error) 

 # %% 

repr_reps = ['rNA', 'rALL']
with open(repeat_guides, 'r') as fp:
    lines = fp.readlines()
    for line in lines:
        if line.startswith(">"):
            repr_reps.append(line.strip(">").strip("\n").strip(" "))

# %% 
'''rep2div2spcount'''
'''repeat_guide > diversity_of_assembly_type(s) > spacer_count'''
'''Dictionary of Counts of Unique Spacers shared by different assembly types'''

rep2div2spcount = {}
for rep in repr_reps:
    path = f'{cdhit_dir}/cd-hit-est-{rep}_spacer.out.clstr'        
    with open(path, "r") as fp:
        lines = fp.readlines()
        
        '''clust2assembly'''
        '''cluster > assemblies'''
        '''Dictionary of Spacer Clusters and the Assemblies for the spacers in a cluster originate'''
        clust2assembly = {}
        count = -1
        sp_count = 0
        for line in lines:
            if line.startswith(">") :
                count += 1
                clust_name = "clust_%s" % (count)
                assemblies = []
            if line[0] != '>':
                assembly = line.split(">")[1].split("_")[0]
                assemblies.append(assembly)
                sp_count += 1
            clust2assembly.update({clust_name: assemblies})
            
        '''clust2diversity'''
        '''cluster > diversity_type'''
        '''Dictionary of Cluster Diversity Types'''
        clust2diversity = {}
        for cluster in clust2assembly.keys():
            assemblies = clust2assembly.get(cluster)
            ltype = []
            for assembly in assemblies:
                if 'unbinned' in assembly:
                    continue
                elif any(prefix in assembly for prefix in ['ERR', 'SRR']):
                    ltype.append('I')
                elif 'coassembly' in assembly:
                    ltype.append('C')
                elif 'bin' in assembly:
                    ltype.append('M')
            type = ''.join(sorted(set(ltype)))
            clust2diversity.update({cluster:type})
            
        '''diversity2spcount'''
        '''diversity_category > sp_count'''
        '''Dictionary of assembly type(s) and representive spacers shared (spacer clusters)'''
        diversity2spcount = {}
        C = 0 # coassembly (exclusivley)
        I = 0 # indiv assembly (exclusivley)
        M = 0 # mag (exclusivley)
        CI = 0 # coassembly & indiv (in common)
        CM = 0 # coassemlby & mag (in common)
        IM = 0 # indiv & mag (in common)
        CIM = 0 # all (in common)
        for cluster in clust2diversity.keys():
            if clust2diversity.get(cluster) == 'C':
                C += 1
            elif clust2diversity.get(cluster) == 'I':
                I += 1
            elif clust2diversity.get(cluster) == 'M':
                M += 1
            elif clust2diversity.get(cluster) == 'CI':
                CI += 1
            elif clust2diversity.get(cluster) == 'CM':
                CM += 1
            elif clust2diversity.get(cluster) == 'IM':
                IM += 1
            elif clust2diversity.get(cluster) == 'CIM':
                CIM += 1   
        diversity2spcount = {'C':C,'I':I,'M':M,'CI':CI,'CM':CM,'IM':IM,'CIM':CIM}
        
    rep2div2spcount[rep] = diversity2spcount

    
# %%
path_out = outdir + '/rep2div2spcount.tsv'
with open(path_out, "w") as fp_out:
    header = ['repeat_guide', 'assembly_type_diversity', 'spacer_count']
    fp_out.write("\t".join(header).strip("\t") + "\n")
    for rep in rep2div2spcount.keys():
        #for div in rep2div2spcount.get(rep):
        for div in ['C','I','CI']:
            if div == 'C':
                assembly_type_diversity = 'Co-assembly'
            elif div == 'I':
                assembly_type_diversity = 'Individual Assembly'
            elif div == 'CI':
                assembly_type_diversity = 'In Common'
            spcount = rep2div2spcount[rep].get(div)
            new_row = [rep, assembly_type_diversity, str(spcount)]
            fp_out.write("\t".join(new_row).strip("\t") + "\n")
            