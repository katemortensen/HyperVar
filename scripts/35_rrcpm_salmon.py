#!/usr/bin/env python
# coding: utf-8

# In[21]:


# Author: Kate Mortensen
# Last Modified: 2/28/2023
Purpose = "Purpose: to find relative abundance of contigs from MAGs, Individual Assemblies, and Co-assembly."


# %%
import pandas as pd 
import numpy as np
import os
import argparse


# ## Relative Abundance Profile of Contigs

# %%
# Initialize parser
parser = argparse.ArgumentParser()
 
# Adding optional argument
parser.add_argument("-wkdir", "--wkdir", help = "path to working directory")
parser.add_argument("-subwkdir", "--subwkdir", help = "path to binner working directory")
parser.add_argument("-abund_dir", "--abund_dir", help = "path to CONTIG_RELATIVE_ABUNDANCE output dir")
parser.add_argument("-crtdir", "--crtdir", help = "path to CRISPRone output dir")
parser.add_argument("-out", "--out", help = "path to save output (including file name)")

# Read arguments from command line
args = parser.parse_args()
wkdir = args.wkdir
subwkdir = args.subwkdir
abund_dir = args.abund_dir 
CRISPRONE_CRT = args.crtdir
out = args.out

# print("working directory: ", wkdir )
# print("sub-working directory: ", subwkdir )
# print("CRISPRone results directory: ", CRISPRONE_CRT)

# %%

# wkdir = "/home/kmorten/ADMB/"
# # subwkdir = "/home/kmorten/ADMB/"
# subwkdir = f"/home/kmorten/ADMB/METAWRAP_MAGS"
# # subwkdir = f"/home/kmorten/ADMB/MAGSCOT_MAGS"
# # subwkdir = f"/home/kmorten/ADMB/DASTOOL_MAGS"
# CRISPRONE_CRT = f"{subwkdir}/CRISPRONE/rNA"
# abund_dir = f"{subwkdir}/CONTIG_RELATIVE_ABUNDANCE"
# out = abund_dir + '/relative_contig_abundance_profiles_crt_salmon.tsv'

# ### Step 1 - CRT contig counts

# %%
'''src_dict'''
'''sample->read_count'''
'''src (sample read counts) = total read counts for each sample '''
src_path = abund_dir + '/sample_reads_count.txt'
src_df = pd.read_csv( src_path , sep="\t")
src_dict = dict(zip(src_df['sample'], src_df['total_reads']))
samples = list(src_df['sample'])

# %%
'''assembly2crt_count'''
'''assembly->assembly_type,contig->crt_count'''
assembly2crt_count = {}
nodeIDs_path = abund_dir + '/nodeIDs.txt'
nodeIDs_file = open(nodeIDs_path, 'r')
for line in nodeIDs_file:
    linep =  line.strip("\n").split("\t")
    assembly = linep[0]
    contig = linep[1]
    assembly_type = linep[2]
    if assembly2crt_count.get(assembly) == None:
        assembly2crt_count[assembly] = {'assembly_type':assembly_type,'contig':{contig:{'crt_count':0}}}
    else :
        assembly2crt_count[assembly]['contig'][contig] = {'crt_count':0} # check this works
assemblies = os.listdir(CRISPRONE_CRT)
for assembly in assemblies:
    path_crt = CRISPRONE_CRT + '/%s/%s.crt' % (assembly,assembly)
    if assembly2crt_count.get(assembly) != None and os.stat(path_crt).st_size != 0:
        contig2crt_tmp = {}
        crt_file = open(path_crt, 'r')
        contig = 'place_holder'
        for line in crt_file :
            if line[0:3] == 'SEQ': # new contig (SEQ)
                contig = line.split(' ')[1]
            if line[0:6] == 'CRISPR': # (CRISPR)
                crt_count = assembly2crt_count[assembly]['contig'][contig].get('crt_count') + 1
                assembly2crt_count[assembly]['contig'][contig]['crt_count'] = crt_count


# ### Step 2 - salmon output to calculate RCPM contig profiles for each sample 
# 
# #### RCPM = (# mapped reads)(10^6)/(total reads in sample)(contig length)

# In[25]:
# Dictionary creation for Relative Abundance Profile creation

'''node2bin'''
'''contig->assembly'''
# node info
nodeIDs_path = abund_dir + '/nodeIDs.txt'
df_nodes = pd.read_csv(nodeIDs_path , sep="\t")
# lookup assembly bin from contig dictionary
df_nodes_tmp = df_nodes[df_nodes['assembly_type']=='mag']
node2bin = dict(zip(df_nodes_tmp['contig'],df_nodes_tmp['assembly']))
assert len(df_nodes_tmp['contig']) == len(list(set(df_nodes_tmp['contig'])))


# In[26]:
# Relative abundance profiles of contigs across samples 

# start df out
cols = ["assembly_type", "assembly", "contig", "crt"] + samples

# quant runs 
#quant_runs = ['coassembly', 'reassembled_bins_concat'] + samples
quant_runs = ['coassembly', 'mag'] + samples

# path to output file
#path_out = abund_dir + '/relative_contig_abundance_profiles_crt_salmon.tsv'
path_out = out


with open(path_out, "w") as fp_out:
    # write header
    header = "\t".join(cols).strip("\t") + "\n"
    fp_out.write(header)

    for q in quant_runs :

        assembly_type = q
        # if q[0] == 'E': 
        if any(prefix in q[0:3] for prefix in ['ERR', 'SRR']):
            assembly_type = 'individual'
        # if q == 'reassembled_bins_concat':
        #     assembly_type = 'reassembled'

        sample2node2rcpm = {}
        type_contigs = []
        type_contigs_accounted = 0 # write contigs to list only 1x

        for s in samples:
            #print("Quant Run: ", q, "Sample: ", s)
            if q == 'mag':
                path_in = subwkdir + '/SALMON_bins_concat/alignment_files/%s.quant/quant.sf' % (s)
            else :
                path_in = wkdir + '/SALMON_%s/alignment_files/%s.quant/quant.sf' % (q,s)
        
            if os.path.exists(path_in):
                
                quantsf = open(path_in, 'r')
                
                # tot reads in sample
                tot_sample_reads = float(src_dict.get(s))
                
                header_accounted = 0 # skip file header
                samplenode_accounted = 0 # sample2node2rcmp initiate new entry
                
                for line in quantsf:
                    if header_accounted != 0:
                        # parse line
                        linep = line.strip('\n').split('\t')

                        # lookup real contig name 
                        #contig = ""
                        contig = linep[0]
                        
                        # accounting of assembly_type contigs
                        if type_contigs_accounted == 0:
                            type_contigs.append("%s-%s" % (assembly_type,contig))
                        
                        # number of reads mapped to contig from sample
                        qty_mapped_reads = float(linep[4])
                        # scale factor 
                        M = 10**6 
                        # length of contig
                        contig_lengths = float(linep[1] )
                        # RCPM
                        rcpm_tmp = (qty_mapped_reads * M) / (tot_sample_reads * contig_lengths)
                        # dictionary
                        if samplenode_accounted != 0:
                            sample2node2rcpm[s].update({contig:rcpm_tmp})
                        else:
                            # sample accounted for 
                            sample2node2rcpm[s] = {contig:rcpm_tmp}
                            samplenode_accounted += 1
                    else :
                        # header accounted for
                        header_accounted += 1

                # contigs accounted for 
                type_contigs_accounted += 1

                # check for duplicates
                assert len(type_contigs) == len(set(type_contigs))
        
                
        if os.path.exists(path_in):
                    
            # assembly & crt info 
            for type_contig in type_contigs:
                contig = type_contig.split("-")[1]
                # assembly
                assembly = q
                if q == 'mag':
                    assembly = node2bin.get(contig)
                # crt count
                crt_count = str(assembly2crt_count[assembly]['contig'][contig].get('crt_count'))
                
                assert assembly_type == assembly2crt_count[assembly].get('assembly_type')
                new_row = [assembly_type,assembly, contig, crt_count]
                rcpm_vec = [str(sample2node2rcpm[s][contig]) for s in samples]
                new_row = new_row + rcpm_vec
                new_row_string = "\t".join(new_row).strip("\t") + "\n"
                fp_out.write(new_row_string)

    del sample2node2rcpm





# %%
