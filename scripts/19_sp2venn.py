#!/usr/bin/env python
# coding: utf-8

# %% 

# Author: Kate Mortensen
# Date: 04-29-2022
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
 
# Adding optional argument
parser.add_argument("-input", "--input", help = "path to cd-hit-est-r#-spacer.out.clstr")
parser.add_argument("-output", "--output", help = "path to save output including file name.")

# Read arguments from command line
args = parser.parse_args()
input = args.input
output = args.output
print("Input File: ", input)
print("Output File: ", output)



    
# %% 
# input= "/home/kmorten/ADMB/CD-HIT-EST/cd-hit-est-rNA_spacer.out.clstr"
# output= "/home/kmorten/ADMB/CD-HIT-EST/sp2venn_rNA.out"

# input = "/home/kmorten/ADMB/DASTOOL_MAGS/CD-HIT-EST/cd-hit-est-rNA_spacer.out.clstr"
# output = "/home/kmorten/ADMB/DASTOOL_MAGS/CD-HIT-EST/sp2venn_rNA.out"

# input = "/home/kmorten/ADMB/MAGSCOT_MAGS/CD-HIT-EST/cd-hit-est-rNA_spacer.out.clstr"
# output = "/home/kmorten/ADMB/MAGSCOT_MAGS/CD-HIT-EST/sp2venn_rNA.out"

    
# %%

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

# Venn Diagram - diversity of Spacer Clusters

unbinned_count = 0
bin_count = 0
indiv_count = 0
coassem_count = 0
bin_coassem_count = 0 # bin & coassembly
bin_indiv_count = 0 # bin & indiv assemlby 
coassem_indiv_count = 0 # indiv & coassembly
common_count = 0 # bin & indiv & coassembly 

for cluster in list(clust2assembly.keys()) :
    unbinned_temp = 0
    bin_temp = 0 # bin with spacer
    indiv_temp = 0 # individual assembly with spacer
    coassem_temp = 0 # coassembly with spacer
    assemblies = list(clust2assembly.get(cluster))
    for assembly in assemblies :
        if 'unbinned' in assembly: # for unbinned
            unbinned_temp += 1
        elif 'bin' in assembly : # for bin
            bin_temp += 1
        elif 'coassembly' in assembly: # for coassembly
            coassem_temp += 1
        elif any(prefix in assembly for prefix in ['ERR', 'SRR']) : # for individual
            indiv_temp += 1
    
    if bin_temp > 0 and indiv_temp > 0 and coassem_temp > 0: # spacer in common
        common_count += 1
    elif bin_temp > 0 and indiv_temp > 0 and coassem_temp==0  :
        bin_indiv_count += 1
    elif bin_temp > 0 and coassem_temp > 0 and indiv_temp==0 :
        bin_coassem_count += 1
    elif coassem_temp > 0 and indiv_temp > 0 and bin_temp==0 : 
        coassem_indiv_count += 1

    elif bin_temp > 0 and indiv_temp == 0 and coassem_temp ==0: # bin only has spacer
        bin_count += 1
    elif indiv_temp > 0 and bin_temp == 0 and coassem_temp == 0: # indiv only has spacer
        indiv_count += 1
    elif coassem_temp > 0 and indiv_temp == 0 and bin_temp == 0 : # coassembly only has spacer
        coassem_count += 1
    elif unbinned_temp > 1 and indiv_temp == 0 and bin_temp == 0 and coassem_temp ==0:
        unbinned_count += 1
assert unbinned_count + bin_count + indiv_count + coassem_count + bin_coassem_count + bin_indiv_count + coassem_indiv_count + common_count == len(list(clust2assembly.keys()))

print("The number of clustered spacers and the assemblies from with the spacers originate from within cluster.")
print(f"# of Clusters with Spacers from Indiv Assemblies (exclusively): {indiv_count}")
print(f"# of Clusters with Spacers from Coassembly exclusively: {coassem_count}")
print(f"# of Clusters with Spacers from MAGs (bins) exclusively: {bin_count}")
print(f"# of Clusters with Spacers from both Indiv Assemblies & Co-assembly, excluding MAGs : {coassem_indiv_count}")
print(f"# of Clusters with Spacers from both Indiv Assemblies & MAGs, excluding Co-assembly: {bin_indiv_count}")
print(f"# of Clusters with Spacers from both Coassembly & MAG Assemblies, excluding Iniv Assemblies: {bin_coassem_count}")
print(f"# of Clusters with Spacers from both MAGs & Indiv Assemblies & Coassembly: {common_count}")
print(f"Total # of Spacer Clusters: {len(list(clust2assembly.keys()))}")



# %%

'''clust2diversity'''
'''cluster > diversity_type'''
'''Dictionary of Cluster diversity Types'''

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

    
# %% 
'''sp2clust'''
'''spacer > cluster(s) for which spacer occurs > cluster_diversity_type'''
'''Dictionary of Spacers, the Clusters they belong to, and Cluster Diversity (Type)'''

sp2clust = {}
with open(input) as fp:
    lines = fp.readlines()
    cluster = None
    for line in lines:
        if line[0] == '>':
            cluster = "clust_%s" % (line.split(' ')[1].strip('\n'))
        if line[0] != '>':
            assembly = line.split('>')[1].split('_')[0]
            spacer = line.split(assembly)[1].split('...')[0].strip('_')
            diversity_type = clust2diversity.get(cluster)
            if sp2clust.get(spacer) == None:
                sp2clust[spacer] = {cluster:diversity_type}
            else:
                sp2clust[spacer].update({cluster:diversity_type})

# %%
'''sp2fidelity'''
'''spacer > cluster diversity type occurrences (for which the spacer occurs)'''   
   
sp2fidelity = {}
for spacer in sp2clust.keys():
    if sp2fidelity.get(spacer) == None:
        sp2fidelity[spacer] = {'I':0,'C':0,'M':0,'CI':0,'IM':0,'CM':0,'CIM':0}
    for cluster in sp2clust[spacer].keys():
        fidelity_type = sp2clust[spacer].get(cluster)
        new_value = sp2fidelity[spacer].get(fidelity_type) + 1
        sp2fidelity[spacer][fidelity_type] = new_value

# %%
# spacer accounting
spacers = []
with open(input) as fp:
    for line in lines:
        if line[0] != '>':
            assembly = line.split('>')[1].split('_')[0]
            spacer = line.split(assembly)[1].split('...')[0].strip('_')
            spacers.append(spacer)
total_spacers_redun = len(spacers)
total_spacers_nonredun = len(set(spacers))

assert len(clust2diversity.keys())<=total_spacers_nonredun<=total_spacers_redun 

'''fidelityDict'''
'''cluster_fidelity_type > quantity of spacers that experience respective fidelity type'''
fidelityDict = {'I':0,'C':0,'M':0,'CI':0,'IM':0,'CM':0,'CIM':0}
for spacer in sp2fidelity.keys():
    for type in sp2fidelity[spacer].keys():
        # cluster diversity type occurences experienced by spacer
        type_occurrence = sp2fidelity[spacer].get(type)
        new_value = fidelityDict.get(type) + type_occurrence
        fidelityDict[type] = new_value
        
print("The number of spacers that belong to a particular category of spacer diversity.") 
print(f"# of Spacers from Indiv Assemblies (exclusively) comprised clusters: {fidelityDict.get('I')}")
print(f"# of Spacers from Coassembly (exclusively) comprised clusters: {fidelityDict.get('C')}")
print(f"# of Spacers from MAG (exclusively) comprised clusters: {fidelityDict.get('M')}")
print(f"# of Spacers from both Indiv Assemblies & Coassemlby comprised clusters: {fidelityDict.get('CI')}")
print(f"# of Spacers from both Indiv Assemblies & MAG comprised clusters: {fidelityDict.get('IM')}")
print(f"# of Spacers from both Coassembly & MAGs Assemblies comprised clusters: {fidelityDict.get('CM')}")
print(f"# of Spacers from both Indiv Assemblies & Coassembly & MAGs comprised clusters: {fidelityDict.get('CIM')}")
print(f"Total # of Spacers (redundant): {total_spacers_redun}")
print(f"Total # of Spacers (non-redundant): {total_spacers_nonredun}")
print(f"Total # of Unique Representative Spacers (Spacer Clusters)): {len(clust2diversity.keys())}")

# %% 

with open(output, "w") as f:
    print("The number of clustered spacers and the assemblies from with the spacers originate from within cluster.", file=f)
    print(f"# of Clusters with Spacers from Indiv Assemblies (exclusively): {indiv_count}", file=f)
    print(f"# of Clusters with Spacers from Coassembly (exclusively): {coassem_count}", file=f)
    print(f"# of Clusters with Spacers from MAGs (bins) (exclusively): {bin_count}", file=f)
    print(f"# of Clusters with Spacers from both Indiv & Coassembly Assemblies: {coassem_indiv_count}", file=f)
    print(f"# of Clusters with Spacers from both Indiv & MAG Assemblies: {bin_indiv_count}", file=f)
    print(f"# of Clusters with Spacers from both Coassembly & MAG Assemblies: {bin_coassem_count}", file=f)
    print(f"# of Clusters with Spacers from both MAGs & Indiv Assemblies & Coassembly: {common_count}", file=f)
    print(f"Total # of Spacer Clusters: {len(list(clust2assembly.keys()))}", file=f)
    print("\n", file=f)
    print("The number of spacers that belong to a particular category of spacer diversity.", file=f) 
    print(f"# of Spacers from Indiv Assemblies (exclusively) comprised clusters: {fidelityDict.get('I')}", file=f)
    print(f"# of Spacers from Coassembly (exclusively) comprised clusters: {fidelityDict.get('C')}", file=f)
    print(f"# of Spacers from MAG (exclusively) comprised clusters: {fidelityDict.get('M')}", file=f)
    print(f"# of Spacers from both Indiv Assemblies & Coassemlby comprised clusters: {fidelityDict.get('CI')}", file=f)
    print(f"# of Spacers from both Indiv Assemblies & MAG comprised clusters: {fidelityDict.get('IM')}", file=f)
    print(f"# of Spacers from both Coassembly & MAGs Assemblies comprised clusters: {fidelityDict.get('CM')}", file=f)
    print(f"# of Spacers from both Indiv Assemblies & Coassembly & MAGs comprised clusters: {fidelityDict.get('CIM')}", file=f)
    print(f"Total # of Spacers (redundant): {total_spacers_redun}", file=f)
    print(f"Total # of Spacers (non-redundant): {total_spacers_nonredun}", file=f)
    print(f"Total # of Unique Representative Spacers (Spacer Clusters)): {len(clust2diversity.keys())}", file=f)


# %%
