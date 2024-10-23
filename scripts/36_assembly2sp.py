#!/usr/bin/env python
# coding: utf-8

# In[2]:
# Author: Kate Mortensen
# Last Modified: 5/3/2023
Purpose = "Report spacer counts for assemblies and assembly types."


# In[5]:

import pandas as pd
import numpy as np
import argparse
import os


# In[5]:
# Initialize parser
parser = argparse.ArgumentParser()
 
# Arguments
parser = argparse.ArgumentParser(description=Purpose)
parser.add_argument("-cdhit_dir", "--cdhit_dir", help = "path to directory where CD-HIT-EST output exists (ex: \'/usr/home/wkdir/CD-HIT-EST\')", required=True)
parser.add_argument("-abund_dir", "--abund_dir", help = "path to directory where CONTIG_RELATIVE_ABUNDANCE output exists (ex: \'/usr/home/wkdir/CONTIG_RELATIVE_ABUNDANCE\')", required=True)
parser.add_argument("-repeat_guides", "--repeat_guides", help = "path to file containing repeat guides used to find CRISPR artifacts in CRISPRone (ex: \'/usr/home/wkdir/CD-HIT-EST/top25_representative_repeats.fa\' )", required=True)
parser.add_argument("-outdir", "--outdir", help = "path to save output directory", required=True)
args = parser.parse_args()

cdhit_dir = args.cdhit_dir
abund_dir = args.abund_dir
nodeIDs_file = f'{abund_dir}/nodeIDs.txt'
repeat_guides_file = args.repeat_guides
outdir = args.outdir

try: 
    os.mkdir(outdir) 
except OSError as error: 
    print(error) 
    
print("CD-HIT-EST directory: ", cdhit_dir )
print('CONTIG_RELATIVE_ABUNDANCE directory:', abund_dir)
print("path to repeat_guides file: ", repeat_guides_file )
print("output directory: ", outdir )


# In[6]:

# cdhit_dir = '/home/kmorten/ADMB/MAGSCOT_MAGS/CD-HIT-EST'
# repeat_guides_file = '/home/kmorten/ADMB/MAGSCOT_MAGS/CD-HIT-EST/top25_representative_repeats.fa'
# accessions = '/home/kmorten/ADMB/MAGSCOT_MAGS/SRR_list.txt'
# outdir = '/home/kmorten/ADMB/MAGSCOT_MAGS/results'
# nodeIDs_file = '/home/kmorten/ADMB/MAGSCOT_MAGS/CONTIG_RELATIVE_ABUNDANCE/nodeIDs.txt'

# cdhit_dir = '/home/kmorten/ADMB/DASTOOL_MAGS/CD-HIT-EST'
# accessions = '/home/kmorten/ADMB/SRR_list.txt'
# repeat_guides = "/home/kmorten/ADMB/DASTOOL_MAGS/CD-HIT-EST/top25_representative_repeats.fa"
# outdir = "/home/kmorten/ADMB/DASTOOL_MAGS/results"
# nodeIDs_file = '/home/kmorten/ADMB/DASTOOL_MAGS/CONTIG_RELATIVE_ABUNDANCE/nodeIDs.txt'


# try: 
#     os.mkdir(outdir) 
# except OSError as error: 
#     print(error)  


# #### Spacer Count

# %%
repr_reps = ['rNA', 'rALL']
with open(repeat_guides_file, 'r') as fp:
    lines = fp.readlines()
    for line in lines:
        if line.startswith(">"):
            repr_reps.append(line.strip(">").strip("\n").strip(" "))

all_assemblies = []
with open(nodeIDs_file, 'r') as fp:
    next(fp)
    for line in fp:
        all_assemblies.append(line.split('\t')[0])
all_assemblies = set(all_assemblies)

# In[4]:
'''assembly2sp'''
'''assembly > assembly_type, repeat_guide > spacer_count'''

assembly2sp = {}
for assembly in all_assemblies:
    assembly_type = 'Individual'
    if any(prefix in assembly.split('.')[0] for prefix in ['bin', 'coassembly']):
        assembly_type = 'Co-assembled'
    assembly2sp[assembly] = {'assembly_type':assembly_type}
    for rep_guide in repr_reps:
        if assembly2sp[assembly].get('rep_guide') == None:
            assembly2sp[assembly].update({'rep_guide':{rep_guide:{'sp_count':0}}})
        else:
            assembly2sp[assembly]['rep_guide'].update({rep_guide:{'sp_count':0}})
   
for rep_guide in repr_reps:
    clstr_file = cdhit_dir + '/concat_' + rep_guide + '_spacer.seq' # (ex: /home/kmorten/oralMB/PRJEB45779/CD-HIT-EST/concat_rNA_spacer.seq)
    fp = open(clstr_file, 'r')
    for line in fp:
        if line[0] == '>':
            assembly = line.split('_')[0].strip('>')
            sp_count = assembly2sp[assembly]['rep_guide'][rep_guide].get('sp_count') + 1
            assembly2sp[assembly]['rep_guide'][rep_guide].update({'sp_count': sp_count})

total_spacers1 = 0
for assembly in assembly2sp.keys():
    for rep in assembly2sp[assembly]['rep_guide'].keys():
        total_spacers1 += assembly2sp[assembly]['rep_guide'][rep_guide].get('sp_count')
  
# %% 
'''assembly2spcount'''
'''assembly > assembly_type, repeat_guide > spacer_count'''
'''(same as assembly2sp dictionary, except parsing different cd-hit-est output file to gather spacer counts)'''
assembly2spcount = {}
for assembly in all_assemblies:
    if any(prefix in assembly.split('.')[0] for prefix in ['bin']):
        assembly_type = 'MAG'
    elif any(prefix in assembly.split('.')[0] for prefix in ['coassembly']):
        assembly_type = 'Co-assembled'
    if any(prefix in assembly.split('.')[0] for prefix in ['ERR', 'SRR']):
        assembly_type = 'Individual'
    assembly2spcount[assembly] = {'assembly_type':assembly_type}
    for rep_guide in repr_reps:
        assembly2spcount[assembly][rep_guide] = {'sp_count':0}
            
for rep_guide in repr_reps:
    path = f'{cdhit_dir}/cd-hit-est-{rep_guide}_spacer.out.clstr'        
    with open(path) as fp:
        lines = fp.readlines()
        count = 0
        for line in lines:
            if line.startswith(">") :
                if count == 1:
                    assemblies = list(set(assemblies))
                    for assembly in assemblies:
                        new_spcount = assembly2spcount[assembly][rep_guide].get('sp_count') + 1
                        assembly2spcount[assembly][rep_guide].update({'sp_count':new_spcount})
                assemblies = []
                count = 0
            if line[0] != '>':
                assembly = line.split(">")[1].split("_")[0]
                assemblies.append(assembly)
                count = 1
# %% 
'''check point'''
for assembly in assembly2sp.keys():
    for rep in assembly2sp[assembly]['rep_guide'].keys():
        spcount1 = assembly2sp[assembly]['rep_guide'][rep_guide].get('sp_count')
        spcount2 = assembly2spcount[assembly][rep_guide].get('sp_count')
        spcount1 == spcount2
        
        

# %%
'''assemblytype2spcount'''
'''assembly_type > repeat_guide > spacer_count'''

assemblytype2spcount = {}

for assembly_type in ['Individual', 'Co-assembled']:
    for rep_guide in repr_reps:
        if assemblytype2spcount.get(assembly_type) == None:
            assemblytype2spcount[assembly_type] = {'rep_guide':{rep_guide:{'sp_count':0}}}
        else:
            assemblytype2spcount[assembly_type]['rep_guide'][rep_guide] = {'sp_count':0}

for assembly in list(assembly2sp.keys()):
    assembly_type = str(assembly2sp[assembly].get('assembly_type'))
    for rep_guide in list(assembly2sp[assembly]['rep_guide'].keys()):
        sp_count = int(assembly2sp[assembly]['rep_guide'][rep_guide].get('sp_count'))
        new_sp_count = sp_count + assemblytype2spcount[assembly_type]['rep_guide'][rep_guide].get('sp_count')
        assemblytype2spcount[assembly_type]['rep_guide'][rep_guide]['sp_count'] = new_sp_count

total_spacers2 = 0
for assembly_type in assemblytype2spcount.keys():
    for rep in assemblytype2spcount[assembly_type]['rep_guide'].keys():
        total_spacers2 += assemblytype2spcount[assembly_type]['rep_guide'][rep_guide].get('sp_count')

assert total_spacers1 == total_spacers2

# %% 
path_out1 = outdir + '/assembly2sp_summary.tsv'
with open(path_out1, "w") as fp_out:
    header = ['assembly', 'assembly_type', 'repeat_guide', 'spacer_count']
    fp_out.write("\t".join(header).strip("\t") + "\n")
    for assembly in list(assembly2sp.keys()):
        assembly_type = str(assembly2sp[assembly].get('assembly_type'))
        for rep_guide in list(assembly2sp[assembly]['rep_guide'].keys()):
            sp_count = str(assembly2sp[assembly]['rep_guide'][rep_guide].get('sp_count'))
            new_row = [assembly, assembly_type, rep_guide, sp_count]
            fp_out.write("\t".join(new_row).strip("\t") + "\n")
              
# %%
path_out2 = outdir + '/assemblytype2sp_summary.tsv'
header_accounted = 0
with open(path_out2, "w") as fp_out:
    if header_accounted == 0:
        header = ['assembly_type', 'repeat_guide', 'spacer_count']
        fp_out.write("\t".join(header).strip("\t") + "\n")
    for assembly_type in list(assemblytype2spcount.keys()):
        for repeat_guide in list(assemblytype2spcount[assembly_type]['rep_guide'].keys()):
            spacer_count = str(assemblytype2spcount[assembly_type]['rep_guide'][repeat_guide].get('sp_count'))
            #print('assembly_type: ', assembly_type, 'repeat_guide: ', repeat_guide, "spacer_count: ", spacer_count)
            new_row = [assembly_type, repeat_guide, spacer_count]
            fp_out.write("\t".join(new_row).strip("\t") + "\n")      
            
# %% 
path_out3 = outdir + '/assembly2sp_count_sorted.tsv'
header = ['assembly', 'assembly_type', 'spacer_count']
df = pd.DataFrame(columns = header, dtype = object)
for assembly in list(assembly2sp.keys()):
    assembly_type = str(assembly2sp[assembly].get('assembly_type'))
    sp_count = int(assembly2sp[assembly]['rep_guide']['rNA'].get('sp_count'))
    new_row = [assembly, assembly_type, sp_count]
    df.loc[len(df.index)] = new_row
df1 = df.sort_values(by='spacer_count', ascending=False).reset_index(drop=True)
df1.to_csv(path_out3, sep = '\t', index=False)
            
