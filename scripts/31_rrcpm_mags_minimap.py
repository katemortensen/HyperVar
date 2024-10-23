#!/usr/bin/env python
# coding: utf-8

# In[21]:


# Author: Kate Mortensen
# Last Modified: 8/15/2023
Purpose = "Purpose: to find relative abundance of contigs from MAGs from minimap sam file."


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
parser.add_argument("-aln_dir", "--aln_dir", help = "path to directory with sam files")
parser.add_argument("-wkdir", "--wkdir", help = "path to working directory")
parser.add_argument("-subwkdir", "--subwkdir", help = "path to binner working directory")
parser.add_argument("-abund_dir", "--abund_dir", help = "path to contig relative abundance directory")
parser.add_argument("-crtdir", "--crtdir", help = "path to CRISPRone output dir")
parser.add_argument("-out", "--out", help = "path to save output (including file name)")

# Read arguments from command line
args = parser.parse_args()
# wkdir = args.wkdir
subwkdir = args.subwkdir
aln_dir = args.aln_dir
abund_dir = args.abund_dir
CRISPRONE_CRT = args.crtdir
out = args.out


# %% 

# subwkdir = "/home/kmorten/ADMB/METAWRAP_MAGS/"
# subwkdir = "/home/kmorten/ADMB/DASTOOL_MAGS/"
# abund_dir = f"{subwkdir}/CONTIG_RELATIVE_ABUNDANCE"
# aln_dir = f"{subwkdir}/MINIMAP_bins_concat/"
# CRISPRONE_CRT = f"{subwkdir}/CRISPRONE/rNA"
# out = f'{abund_dir}/relative_contig_abundance_profiles_crt_minimap.tsv'

# %% 
# sam_files = []
# for x in os.listdir(aln_dir):
#     if x.endswith(".sam"):
#         sam_files.append(x)
        
# sam_file = aln_dir + '/' + sam_files[-1]

# %% 

def contig_rrcpm_by_sample(sam_file):

    '''Parse SAM file'''
    '''contig = reference = rname'''
    '''readID = query = qname'''

    SQdict = {}
    qname2count = {}
    rname2qname2counts = {}

    with open(sam_file, "r") as fp: 

        lines = fp.readlines()
        switch = 0
        count = 0
        for line in lines:
            
            '''SQdict'''
            '''SN > LN'''
            '''Reference Segment ID (SN) > Segment Length (LN) '''
            '''Maps segments IDs from MAGs to segment length (data from sam file)'''

            if line.startswith('@SQ'):
                linep = line.strip('\n').split('\t')
                sq = linep[1].split('SN:')[1]
                ln = int(linep[2].split('LN:')[1])
                SQdict[sq] = ln
    
            if line.startswith('@PG'):
                switch = 1
            if switch == 1 and line[0] != '@':
                linep = line.strip('\n').split('\t')
                qname = linep[0] # query template name
                rname = linep[2] # reference seq name 
                
                if rname != '*' and qname != '*':
                    count += 1
                    
                    '''qname2count'''
                    '''qname > count'''
                    '''Dictonary of reads to mapping occurences'''
                    
                    if qname2count.get(qname) == None:
                        qname2count[qname] = 1
                    else:
                        new_count = qname2count.get(qname) + 1
                        qname2count[qname] = new_count
                                        
                    '''rname2qname2counts'''
                    '''rname > qname > counts'''
                    '''Dictonary of reference segments (mag contigs), read IDs of mapped reads, and total instances read was mapped to contig, adn redistributed counts'''

                    if rname2qname2counts.get(rname) == None:
                        rname2qname2counts[rname] = {qname: 1}
                    elif rname2qname2counts[rname].get(qname) == None:
                        rname2qname2counts[rname][qname] = 1
                    else:
                        new_count = 1 + rname2qname2counts[rname].get(qname)
                        rname2qname2counts[rname][qname] = new_count

        sum1 = sum(qname2count.values())
        sum2 = sum(sum(rname2qname2counts[rname].values()) for rname in rname2qname2counts.keys())
        assert count == sum1 == sum2

    '''Redistributed Read Count (RRC)'''
    '''Contig RRC = (contig readID count)/(total readID count across all contigs in all mags)'''
    '''contig2rrc'''
    '''contig > rrc'''

    contig2rrc = {}
    for rname in rname2qname2counts.keys():
        contig_rrc = 0 
        for qname in rname2qname2counts[rname].keys():
            readID_counts = rname2qname2counts[rname].get(qname)
            total_readID_counts = qname2count.get(qname)
            readID_rrc = readID_counts/total_readID_counts
            contig_rrc += readID_rrc
            contig2rrc[rname] = contig_rrc
        
    '''Redistributed Relative Counts per Million (RRCPM)'''
    '''Contig RRCPM = (# mapped reads_redistributed to contig)(10^6)/(total reads)(context length)'''
    '''contig2rrcpm'''
    '''contig > rrcpm'''

    contig2rrcpm = {}
    total_reads = len(qname2count.keys())
    scale_factor = 10**6
    for rname in contig2rrc.keys():
        context_len = SQdict.get(rname)
        rrc = contig2rrc.get(rname)
        rrcpm = (rrc * scale_factor)/(total_reads * context_len)
        contig2rrcpm[rname] = rrcpm
    
    return contig2rrcpm



'''assembly2crt_count'''
'''assembly->assembly_type,contig->crt_count'''

assembly2crt_count = {}
nodeIDs_path = abund_dir + '/nodeIDs.txt'
nodeIDs_file = open(nodeIDs_path, 'r')
with open(nodeIDs_path) as f:
    next(f)
    for line in f:
        linep =  line.strip("\n").split("\t")
        assembly = linep[0]
        contig = linep[1]
        assembly_type = linep[2]
        if assembly2crt_count.get(assembly) == None:
            assembly2crt_count[assembly] = {'assembly_type':assembly_type,'contig':{contig:{'crt_count':0}}}
        else :
            assembly2crt_count[assembly]['contig'][contig] = {'crt_count':0} 

assemblies = os.listdir(CRISPRONE_CRT)
for assembly in assemblies:
    path_crt = CRISPRONE_CRT + '/%s/%s.crt' % (assembly,assembly)
    if os.stat(path_crt).st_size != 0:
        contig2crt_tmp = {}
        crt_file = open(path_crt, 'r')
        contig = 'place_holder'
        for line in crt_file :
            if line[0:3] == 'SEQ': # new contig (SEQ)
                contig = line.split(' ')[1]
            if line[0:6] == 'CRISPR': # (CRISPR)
                crt_count = assembly2crt_count[assembly]['contig'][contig].get('crt_count') + 1
                assembly2crt_count[assembly]['contig'][contig]['crt_count'] = crt_count


'''Relative Abundance Profile (RAP)'''
'''Contig Relative Abundance Profile (aka Contig RRCPM by Sample)'''
'''Ex: Contig RAP = [RRCPM_SRR1, RRCPM_SRR2, ... RRCPM_SRRN]'''

sam_files = []
for x in os.listdir(aln_dir):
    if x.endswith(".sam"):
        sam_files.append(x)

sample2contig2rrcpm = {}
for sam in sam_files:
    sample = sam.split('_')[0]
    path_in = aln_dir + '/' + sam
    contig2rrcpm = contig_rrcpm_by_sample(path_in)
    for contig in contig2rrcpm.keys():
        if sample2contig2rrcpm.get(sample) == None:
            sample2contig2rrcpm[sample] = {contig:contig2rrcpm.get(contig)}
        else:
            sample2contig2rrcpm[sample][contig]=contig2rrcpm.get(contig)

samples = sample2contig2rrcpm.keys()


# %% 
# Relative abundance profiles of contigs across samples 

# start df out
cols = ["assembly_type", "assembly", "contig", "crt"] + list(samples)

# path to output file
# path_out = abund_dir + "/relative_contig_abundance_profiles_crt_minimap.tsv"
path_out  = out

print('Writting output to: ', path_out)
with open(path_out, "w") as fp_out:
    # write header
    header = "\t".join(cols).strip("\t") + "\n"
    fp_out.write(header)
    for assembly in assembly2crt_count.keys():
        assembly_type = assembly2crt_count[assembly].get('assembly_type')
        if assembly_type == 'mag':
            for contig in assembly2crt_count[assembly].get('contig'):
                crt_count = assembly2crt_count[assembly]['contig'][contig].get('crt_count')
                contig_rrcpms = []
                for sample in sample2contig2rrcpm.keys():
                    if sample2contig2rrcpm[sample].get(contig) != None:
                        rrcpm = sample2contig2rrcpm[sample].get(contig)
                    else:
                        rrcpm = 0
                    contig_rrcpms.append(rrcpm)

                new_row = [assembly_type,assembly, contig, str(crt_count)]
                rcpm_vec = [str(rrcpm) for rrcpm in contig_rrcpms]
                new_row = new_row + rcpm_vec
                new_row_string = "\t".join(new_row).strip("\t") + "\n"
                fp_out.write(new_row_string)
            







# %%
