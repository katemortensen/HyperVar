#!/usr/bin/env python
# coding: utf-8


# Author: Kate Mortensen
# Date: 01-25-2023
# Purpose: Create a dictionary of spacer clusters and their members.

# Ex: {Clust 1 : bin1, SRR1, SRR2, ...}
# Ven Diagram: Spacers from individual assemblies only, co-assembly only, and in common with both individual assemblies and co-assembly

# Purpose: After finding representative repeat-associated spacers, clustering these repeat-associated spacers, find all individual assemblies and bins that are associated with each representative spacer for clustered spacers. 

# In[2]:
import pandas as pd
import numpy as np
import argparse


# In[3]:

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

print("Input File: ", input)
print("Output File: ", output)
print("Repeat Guide: ", rep)

# %% 
# input= "/home/kmorten/ADMB/CD-HIT-EST/cd-hit-est-rNA_spacer.out.clstr"
# rep = input.strip('_spacer.out.clstr').split('-')[-1]
# output= f"/home/kmorten/ADMB/CD-HIT-EST/sp2assembly_{rep}"


# input = "/home/kmorten/ADMB/DASTOOL_MAGS/CD-HIT-EST/cd-hit-est-rNA_spacer.out.clstr"
# rep = input.strip('_spacer.out.clstr').split('-')[-1]
# output = f"/home/kmorten/ADMB/DASTOOL_MAGS/CD-HIT-EST/sp2assembly_{rep}"

# input = "/home/kmorten/ADMB/MAGSCOT_MAGS/CD-HIT-EST/cd-hit-est-rNA_spacer.out.clstr"
# rep = input.strip('_spacer.out.clstr').split('-')[-1]
# output = f"/home/kmorten/ADMB/MAGSCOT_MAGS/CD-HIT-EST/sp2assembly_{rep}"

# %% 
'''clust2sp'''
'''cluster > spacers'''
'''Dictionary of Spacer Clusters and their Spacers'''

clust2spacers = {}
count = -1
sp_count = 0
with open(input) as fp:
    lines = fp.readlines()
    for line in lines:
        if line.startswith(">") :
            count += 1
            clust_name = "clust_%s" % (count)
            spacers = []
        if line[0] != '>':
            spacer = line.split(">")[1]
            spacers.append(spacer)
            sp_count += 1
        clust2spacers.update({clust_name: spacers})

'''clust2assembly'''
'''cluster > assemblies'''
'''Dictionary of Spacer Clusters and the Assemblies for the spacers in a cluster originate'''

clust2assembly = {}
count = -1
sp_count = 0
with open(input) as fp:
    lines = fp.readlines()
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
    spacers = clust2assembly.get(cluster)
    ltype = []
    for spacer in spacers:
        if any(prefix in spacer for prefix in ['ERR', 'SRR']):
            ltype.append('I')
        if 'coassembly' in spacer:
            ltype.append('C')
        if 'bin' in spacer:
            ltype.append('M')
    type = ''.join(sorted(set(ltype)))
    clust2diversity.update({cluster:type})


# In[49]:
'''Convert cluster input file into a dataframe'''
count = -1
with open(input) as fp_in:
    lines = fp_in.readlines()
    file_out = output
    with open(file_out, "w") as fp_out:
        fp_out.write("cluster\tsp_len\tassembly\tcontig\tsp_loci\tscore\tclust_member_qty\tclust_diversity\n")
        for line in lines:
            if line.startswith(">") :
                count += 1
                cluster = "clust_%s" % (count)
            if line[0] != '>':
                member = line.split(" ")
                sp_len = member[0].split("\t")[1].strip("nt,")
                assembly = member[1].strip(">").split("_")[0]
                contig = member[1].strip(">%s_" % (assembly)).split(":")[0]
                sp_loci = member[1].split(">%s_%s" % (assembly,contig))[1]
                sp_loci = sp_loci.strip("...")
                score = member[-1].strip("\n")
                clust_member_qty = len(clust2spacers.get(cluster))
                clust_diversity = clust2diversity.get(cluster)
                
                # Writing data to a file
                new_row = f"{cluster}\t{sp_len}\t{assembly}\t{contig}\t{sp_loci}\t{score}\t{clust_member_qty}\t{clust_diversity}\n"
                fp_out.write(new_row)
    fp_out.close


