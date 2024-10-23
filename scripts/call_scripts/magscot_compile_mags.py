#!/usr/bin/env python3

# Author: Kate Mortensen
# Date: 04-13-2022
Purpose = "Purpose: create magscot bins (to extract specified contig from specified assembly and save to file)"

# %%
from Bio import SeqIO
import pandas as pd
import argparse
from sys import argv
import os.path

# %%
# arguments
parser = argparse.ArgumentParser(description=Purpose)
parser.add_argument('-f', type = str, help = "enter path to assembly fasta file from which you would like to extract the contig", required = True)
parser.add_argument('-i', type = str, help = "enter path to magscot output file (ex: final_assembly.contigs_to_bins.tsv)", required = True)
parser.add_argument('-o', type = str, help = "enter output path", required = True)
args = parser.parse_args()

coassembly = args.f
input_path = args.i
out_path = args.o

# %%
# coassembly = '/home/kmorten/ADMB/COASSEMBLY/final_assembly.fasta'
# input_path = '/home/kmorten/ADMB/MAGSCOT_MAGS/MAGScoT.refined.contig_to_bin.out'
# out_path = '/home/kmorten/ADMB/MAGSCOT_MAGS/magscot_bins'


# %%
# make output directory
isExist = os.path.exists(out_path)
if not isExist:
   os.makedirs(out_path)
   
# %%
record_dict = SeqIO.index(coassembly, "fasta")

# %%
contigs2bins = {}

with open(input_path) as f:
    next(f)
    for line in f:
        linep = line.strip('\n').split("\t")
        bin = linep[0]
        bin_path = out_path + '/' + bin + '.fa'
        contig = linep[1]
        seq = record_dict[contig].seq
        if contigs2bins.get(bin) == None:
            contigs2bins[bin] = {'contigs':{contig:seq}}
        else:
            contigs2bins[bin]['contigs'].update({contig:seq})
            
            
# %%

for bin in contigs2bins.keys():
    fp = out_path + '/'  + bin + '.fa'
    if os.path.exists(fp):
        os.remove(fp)
    with open(fp, "a") as fp_out:
        for contig in contigs2bins[bin].get('contigs'): 
            seq = str(contigs2bins[bin]['contigs'].get(contig))
            SeqIO.write(record_dict[contig], fp_out, 'fasta')
            #fp_out.write('>' + contig + '\n')
            #fp_out.write(seq + '\n')

# Resources:
# https://biopython.org/wiki/SeqRecord

