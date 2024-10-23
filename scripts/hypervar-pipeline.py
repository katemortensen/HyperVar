#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Update: 10/21/2024

# %%
from __future__ import print_function
import sys
import os
from pathlib import Path
import subprocess
import time
from time import gmtime, strftime
import argparse
import statistics
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import csv
from inspect import currentframe, getframeinfo
import multiprocessing
import shutil
from multiprocessing import Pool

# %%
# arguments
Purpose="Explore hypervariable regions in MAG curration."

parser = argparse.ArgumentParser(description=Purpose)
parser.add_argument("-accessions", "--accessions", help = "path to directory of accessions", required=True)
parser.add_argument("-threads", "--threads", help = "threads", required=False)
parser.add_argument("-metawrap", "--metawrap", help = "execute metawrap binning", required=False)
parser.add_argument("-dastool", "--dastool", help = "execute dastool binning", required=False)
parser.add_argument("-magscot", "--magscot", help = "execute magscot binning", required=False)
parser.add_argument("-mags", "--mags", help = "text file with paths to mag directories (if already prepared)", required=False)
parser.add_argument("-salmon", "--salmon", help = "lite-weight mapping estimation via salmon quasi-mapping", required=True)
parser.add_argument("-minimap", "--minimap", help = "fast, lite-weight mapping via minimap", required=True)
parser.add_argument("-out", "--out", help = "path to output directory", required=True)

args = parser.parse_args()

accessions = args.accessions
threads = int(args.threads)
wkdir = args.out

metawrap = args.metawrap
dastool = args.dastool 
magscot = args.magscot 
mags = args.mags
salmon = args.salmon
minimap = args.minimap 

# %% 
        
'''Pre-sets'''
# scripts = os.getcwd()
scripts = os.path.dirname(os.path.abspath(__file__))
bin = '%s/bin' % (scripts.strip('scripts'))
print(f'bin path: {bin}')
python = 'python3'


'''Software'''
crisprone = f'{bin}/CRISPRone-main/crisprone-local-20240911_KM.py'
summarize = f'{bin}/CRISPRpk/metaCRT/summarize-crispr-new.py'
crisprgraph = f'{bin}/CRISPRpk/pipelines/crispr_ann_ort.py'

'''Params'''
mapping_tools = []
if salmon == True:
    mapping_tools.append('salmon')
if minimap == True:
    mapping_tools.append('minimap')

subwkdirs = []
mag_dirs = []
if metawrap == True:
    subwkdir = f'{wkdir}/METAWRAP_MAGS'
    subwkdirs.append(subwkdir)
if dastool == True:
    subwkdir = f'{wkdir}/DASTOOL_MAGS'
    subwkdirs.append(subwkdir)
if magscot == True:
    subwkdir = f'{wkdir}/MAGSCOT_MAGS'
    subwkdirs.append(subwkdir)
if mags != False:
    # MAGSCOT   /path/to/magscot_bins
    # BinnerX   /path/to/binnerx_bins
    for line in open(args.mags, 'r'):
        linep = line.strip('\n').split('\t')
        binning_set_name = linep[0]
        path_to_mags = linep[1]
        subwkdirs.append(f'{wkdir}/{binning_set_name}')
        mag_dirs.append(path_to_mags)

individuals = f'{wkdir}/INDIVIDUALS'
clean_reads = f'{wkdir}/clean_reads'
coassembly_dir = f'{wkdir}/COASSEMBLY'
mag_dirs = [f'{subwkdirs[0]}/BIN_REASSEMBLY/reassembled_bins',f'{subwkdirs[1]}/bins',f'{subwkdirs[2]}/bins']
taxa_dirs = [f'{path}/BIN_CLASSIFICATION' for path in subwkdirs]
annot_dirs = [f'{path}/BIN_ANNOTATION' for path in subwkdirs]
crisprone_dirs = [f'{path}/CRISPRONE' for path in subwkdirs]
cdhit_dirs = [f'{path}/CD-HIT-EST' for path in subwkdirs]
minimap_mag_dirs = [f'{path}/MINIMAP_bins_concat/' for path in subwkdirs]
salmon_dirs = [f'{path}/SALMON_bins_concat/' for path in subwkdirs]
abund_dirs =  [f'{path}/CONTIG_RELATIVE_ABUNDANCE' for path in subwkdirs]
result_dirs = [f'{path}/results' for path in subwkdirs]
crispr_graphs = [f'{path}/CRISPR_GRAPHS' for path in subwkdirs]

samples = open(accessions, "r").read().strip('\n').split('\n')
minimap_indiv_dirs = [f'{wkdir}/MINIMAP_{sample}' for sample in samples]  
minimap_coassembly_dir = [f'{wkdir}/MINIMAP_coassembly']
aln_dirs = minimap_coassembly_dir + minimap_indiv_dirs
suffixes = [i.split('_')[-1] for i in aln_dirs]


# %%
'''Clean Reads'''
'''TODO create python script and option to remove human contamination (or alternative contamination)'''
cmd = f'{scripts}/01_clean_reads.sh {accessions} {individuals}'
subprocess.check_output(cmd, shell=True)

# %%
'''Concatenate Reads for Coassembly'''
cmd = f'''{scripts}/02_concat_reads.sh \
    {wkdir} \
    {accessions} \
    {threads}'''
subprocess.check_output(cmd, shell=True)

# %% 
'''Coassembly and Individual Metagenomic Assemblies'''
cmd = f'{scripts}/03_assembly_reads.sh {wkdir} {threads}'
subprocess.check_output(cmd, shell=True)

# %% 
'''MetaWRAP Binning '''
cmd = f'''{scripts}/04_metawrap.sh \
    {accessions} \
    {wkdir}/METAWRAP_MAGS \
    {threads}'''
subprocess.check_output(cmd, shell=True)

# %%
'''DASTool Binning'''
if dastool==True:
    cmd = f'{scripts}/05_dastool.sh'
    subprocess.check_output(cmd, shell=True)

# %%
'''MAGScoT Binning'''
if magscot==True:
    cmd = f'{scripts}/06_magscot.sh'
    subprocess.check_output(cmd, shell=True)

# %%
'''Bin Taxonomy '''
for i in range(len(subwkdirs)):
    cmd = f'''{scripts}/07_bin_clasification.sh \
        {mag_dirs[i]} \
        {taxa_dirs[i]} \
        {threads}'''
    subprocess.check_output(cmd, shell=True)

# %%
'''Bin Annotation'''
for i in range(len(subwkdirs)):
    cmd = f'''{scripts}/08_bin_annot.sh \
        {mag_dirs[i]} \
        {annot_dirs[i]} \
        {threads}'''
    subprocess.check_output(cmd, shell=True)

# %%
'''CRISPRone (no repeat guide) for coassembly'''
cmd = f'''{scripts}/09_crisprone_coassembly.sh \
    {python} \
    {crisprone}
    {wkdir} \
    {wkdir}/CRISPRONE '''
subprocess.check_output(cmd, shell=True)

'''CRISPRone (no repeat guide) for individual assemblies'''
cmd = f'''{scripts}/10_crisprone_indiv.sh \
    {python} \
    {crisprone} \
    {wkdir} \
    {wkdir}/CRISPRONE \
    {individuals}'''
subprocess.check_output(cmd, shell=True)

for i in range(len(subwkdirs)):
    
    '''CRISPRone (no repeat guide) for MAGs'''
    cmd = f'''{scripts}/11_crisprone_bins.sh \
        {python} \
        {crisprone} \
        {wkdir}/CRISPRONE
        {subwkdirs[i]} \
        {crisprone_dirs[i]} \
        {mag_dirs[i]} \
        {threads}'''
    subprocess.check_output(cmd, shell=True)
    
    '''CRISPRone Summarize'''
    cmd = f'''{scripts}/13_crisprone_summarize.sh \
        {python} \
        {summarize} \
        {subwkdirs[i]} \
        {crisprone_dirs[i]} \
        {threads}'''
    subprocess.check_output(cmd, shell=True)


# %% 
for i in range(len(subwkdirs)):
    
    '''Cluster Assembly Repeats with CD-HIT-EST'''
    # outputs: 
    # <subwkdir>/CD-HIT-EST/cd-hit-est.out
    # <subwkdir>/CD-HIT-EST/cd-hit-est.out.clstr
    cmd = f'''{scripts}/14_cdhit_rep_clustr.sh \
        {cdhit_dirs[i]} \
        {crisprone_dirs[i]} \
        {bin}/cdhit-4.8.1/cd-hit-est'''
    subprocess.check_output(cmd, shell=True)
    
    '''Representative Repeats'''
    # outputs:
    # <subwkdir>/CD-HIT-EST/top25_representative_repeats.fa
    topmost_repeat_guides = 25
    cmd = f'''{python} {scripts}/15_reprsnt_repeats.py \
        --clstr {cdhit_dirs[i]}/cd-hit-est.out.clstr \
        --repcon {cdhit_dirs[i]}/concat_repeatcon.seq \
        --outdir {cdhit_dirs[i]} \
        --topmost_rep_guides {topmost_repeat_guides}'''
    subprocess.check_output(cmd, shell=True)


# %% 
for i in range(len(subwkdirs)):
    
    '''CRISPRone Guided'''
    cmd = f'''{scripts}/16_crisprone_guided.sh \
        {python} \
        {scripts}/call_scripts/extract_contig.py \
        {crisprone} \
        {mag_dirs[i]} \
        {cdhit_dirs[i]}/top25_representative_repeats.fa \
        {cdhit_dirs[i]} \
        {crisprone_dirs[i]} \
        {wkdir} \
        {coassembly_dir}/final_assembly.fasta'''
    subprocess.check_output(cmd, shell=True)
    
    '''CRISPRone Summary Script for Spacer Discovery'''
    '''Note: not working because of summarize script'''
    # output: <subwkdir>/CRISPRONE/<repeat_guide>/<assembly>/<assembly>-spacer.seq
    cmd = f'''{scripts}/17_spdisc.sh \
        {python} \
        {summarize} \
        {crisprone_dirs[i]} '''
    subprocess.check_output(cmd, shell=True)

    '''CD-HIT-EST for Spacers'''
    # output: 
    # <subwkdir>/CD-HIT-EST/concat_<repeat_guide>_spacer.seq
    # <subwkdir>/CD-HIT-EST/cd-hit-est-<repeat_guide>_spacer.out
    # <subwkdir>/CD-HIT-EST/cd-hit-est-<repeat_guide>_spacer.out.clstr
    cmd = f'''{scripts}/18_cdhit_sp_clustr.sh \
        {bin}/cdhit-4.8.1/cd-hit-est \
        {wkdir} \
        {cdhit_dirs[i]} \
        {crisprone_dirs[i]} '''
    subprocess.check_output(cmd, shell=True)

    '''Spacers Venn Diagram by Assembly Type'''
    # output :
    # sp2venn_rNA.summary
    repeat_guide = 'rNA'
    cmd = f'''{python} {scripts}/19_sp2venn.py \
        --input {cdhit_dirs[i]}/cd-hit-est-{repeat_guide}_spacer.out.clstr \
        --output {cdhit_dirs[i]}/sp2venn_{repeat_guide}.summary'''
    subprocess.check_output(cmd, shell=True)
    
    '''Spacer Cluster info'''
    # output: 
    # sp2assembly_rNA.out
    repeat_guide = 'rNA'
    cmd = f'''{python} {scripts}/20_sp2assembly.py \
        --input {cdhit_dirs[i]}/cd-hit-est-{repeat_guide}_spacer.out.clstr \
        --output {cdhit_dirs[i]}/sp2assembly_{repeat_guide}.out'''
    subprocess.check_output(cmd, shell=True)
    
    '''Co-assembly, Individual Assembly, & In Common Unique Spacer Counts'''
    # output: 
    # rep2div2spcount.tsv
    cmd = f'''{python} {scripts}/21_rep2div2spcount.py \
        --cdhit_dir {cdhit_dirs[i]} \
        --repeat_guides {cdhit_dirs[i]}/top25_representative_repeats.fa \
        --outdir {result_dirs[i]}'''
    subprocess.check_output(cmd, shell=True)

    '''Co-assembly, Individual Assembly, & In Common Unique Spacer Counts - Rplots'''
    # outputs: 
    # rep2div2spcount.png/pdf
    cmd = f'''Rscript {scripts}/plot_scripts/rep2div2spcount_plot.R \
        --input {result_dirs[i]}/rep2div2spcount.tsv \
        --outdir {result_dirs[i]}/R_figs'''
    subprocess.check_output(cmd, shell=True)      


# %% 
'''Node IDs & Collect Sample Read Counts - Co-assembly & Individuals'''

'''Collect Sample Read Counts - Co-assembly & Individuals'''
# outputs:
# ../CONTIG_RELATIVE_ABUNDANCE/sample_reads_count.txt
cmd = f'{scripts}/25_tot_sample_reads.sh \
    {accessions} \
    {clean_reads} \
    {wkdir}/CONTIG_RELATIVE_ABUNDANCE'
subprocess.check_output(cmd, shell=True)

'''Node IDs - Co-assembly & Individuals'''
# outputs:
# ../CONTIG_RELATIVE_ABUNDANCE/nodeIDs.txt
cmd = f'''{scripts}/26_nodeIDs_nonmags.sh \
    {wkdir} \
    {wkdir} \
    {coassembly_dir} \
    {individuals} \
    {accessions}'''
subprocess.check_output(cmd, shell=True)


# %%
'''Node IDs & Collect Sample Read Counts- MAGs'''
for i in range(len(subwkdirs)):
    
    '''Collect Sample Read Counts'''
    # outputs:
    # ../CONTIG_RELATIVE_ABUNDANCE/sample_reads_count.txt
    cmd = f'{scripts}/25_tot_sample_reads.sh {abund_dirs[i]}'
    subprocess.check_output(cmd, shell=True)

    '''Node IDs'''
    # outputs:
    # ../CONTIG_RELATIVE_ABUNDANCE/nodeIDs.txt
    cmd = f'''{scripts}/27_nodeIDs.sh \
        {wkdir} \
        {subwkdirs[i]} \
        {coassembly_dir} \
        {individuals} \
        {mag_dirs[i]} \
        {accessions}'''
    subprocess.check_output(cmd, shell=True)
       
    '''Node ID Conversion'''
    # outputs:
    # ../CONTIG_RELATIVE_ABUNDANCE/nodeID_conversion.txt
    cmd = f'''{scripts}/28_nodeID_conversion.sh \
        {subwkdirs[i]} \
        {mag_dirs[i]} '''
    subprocess.check_output(cmd, shell=True)

    '''Concat Bins'''
    # output: bins_concat.fasta 
    cmd = f'''{python} {scripts}/29_concat_bins.py \
        --mag_dir {mag_dirs[i]} \
        --outdir {subwkdirs[i]}/bins_concat '''
    subprocess.check_output(cmd, shell=True)


# %% 
'''WARNING: Mapping scripts are task and memory intensive. \
Each loop requires ~50GB. Avoid parallelizing and multithreading \
(unless your server can handle it).'''


'''Mapping - Sample Reads to Co-assembly'''
'''(via Salmon and/or Minimap)'''

# samples = open(accessions, "r").read().split('\n')
# samples.remove('')
assert len(samples) > 0 

'''Minimap & Samtools'''
# outputs:
# /subwkdir/<accession>_aln.sam 
if minimap == True:  
    
    try: 
        os.mkdir(f'''{wkdir}/MINIMAP_coassembly''') 
    except OSError as error: 
        print(error)  
    
    for sample in samples :
        cmd = f'''{bin}/minimap2/minimap2 -ax sr \
            {coassembly_dir}/final_assembly.fasta \
            {clean_reads}/{sample}/final_pure_reads_1.fastq \
            {clean_reads}/{sample}/final_pure_reads_2.fastq \
            | {bin}/samtools-1.14/samtools \
            sort -o {wkdir}/'MINIMAP_coassembly'/{sample}_aln.sam & '''
        subprocess.check_output(cmd, shell=True) 
        
        '''Coverage'''
        cmd = f'''{bin}/samtools-1.14/samtools depth \
            {wkdir}/'MINIMAP_coassembly'/{sample}_aln.sam \
            > {wkdir}/'MINIMAP_coassembly'/{sample}_aln.coverage'''
        subprocess.check_output(cmd, shell=True) 

''''Mapping (Estimation) with SALMON'''
'''SALMON - light weight reads mapping estimation'''
'''Map clean reads form all samples to concatenated mags'''
# outputs: SALMON_bins_concat
if salmon == True:
    cmd = f'''{scripts}/32_salmon_coassembly.sh \
        {bin}/salmon-latest_linux_x86_64/bin/ \
        {wkdir} \
        {accessions} \
        {clean_reads} \
        {coassembly_dir} \
        {wkdir} \
        {threads}'''
    subprocess.check_output(cmd, shell=True)

# %% 
'''Mapping - Sample Reads to Individual Assemblies'''
'''(via Salmon and/or Minimap)'''

samples = open(accessions, "r").read().split('\n')
samples.remove('')
assert len(samples) > 0 

for i in samples:
    '''Minimap & Samtools'''
    # outputs:
    # /subwkdir/<accession>_aln.sam 
    if minimap == True:  
        
        try: 
            os.mkdir(f'''{wkdir}/MINIMAP_{i}''') 
        except OSError as error: 
            print(error)  
        
        for sample in samples :
            cmd = f'''{bin}/minimap2/minimap2 -ax sr \
                {individuals}/{i}/ASSEMBLY/final_assembly.fasta \
                {clean_reads}/{sample}/final_pure_reads_1.fastq \
                {clean_reads}/{sample}/final_pure_reads_2.fastq \
                | {bin}/samtools-1.14/samtools \
                sort -o {wkdir}/MINIMAP_{i}/{sample}_aln.sam & '''
            subprocess.check_output(cmd, shell=True) 
            
            '''Coverage'''
            cmd = f'''{bin}/samtools-1.14/samtools depth \
                {wkdir}/MINIMAP_{i}/{sample}_aln.sam \
                > {wkdir}/MINIMAP_{i}/{sample}_aln.coverage'''
            subprocess.check_output(cmd, shell=True) 

    ''''Mapping (Estimation) with SALMON'''
    '''SALMON - light weight reads mapping estimation'''
    '''Map clean reads form all samples to concatenated mags'''
    # outputs: SALMON_bins_concat
    if salmon == True:
        cmd = f'''{scripts}/34_salmon_bins_concat.sh \
            {bin}/salmon-latest_linux_x86_64/bin/ \
            {wkdir} \
            {accessions} \
            {clean_reads} \
            {individuals} \
            {wkdir}/SALMON_{i} \
            {threads}'''
        subprocess.check_output(cmd, shell=True)

# %%    
'''Mapping - Sample Reads to MAGS'''
'''(via Salmon and/or Minimap)'''
for i in range(len(subwkdirs)):
    
    samples = open(accessions, "r").read().split('\n')
    samples.remove('')
    
    '''Minimap & Samtools'''
    # outputs:
    # /subwkdir/<accession>_aln.sam 
    if minimap == True:  
        
        try: 
            os.mkdir(minimap_mag_dirs[i]) 
        except OSError as error: 
            print(error)  
        
        for sample in samples :
            cmd = f'''{bin}/minimap2/minimap2 -ax sr \
                {subwkdirs[i]}/bins_concat/bins_concat.fasta \
                {clean_reads}/{sample}/final_pure_reads_1.fastq \
                {clean_reads}/{sample}/final_pure_reads_2.fastq \
                | {bin}/samtools-1.14/samtools \
                sort -o {minimap_mag_dirs[i]}/{sample}_aln.sam & '''
            subprocess.check_output(cmd, shell=True) 
            
            '''Coverage'''
            cmd = f'''{bin}/samtools-1.14/samtools depth \
                {minimap_mag_dirs[i]}/{sample}_aln.sam \
                > {minimap_mag_dirs[i]}/{sample}_aln.coverage'''
            subprocess.check_output(cmd, shell=True) 

    ''''Mapping (Estimation) with SALMON'''
    '''SALMON - light weight reads mapping estimation'''
    '''Map clean reads form all samples to concatenated mags'''
    # outputs: SALMON_bins_concat
    if salmon == True:
        cmd = f'''{scripts}/34_salmon_bins_concat.sh \
            {bin}/salmon-latest_linux_x86_64/bin/ \
            {subwkdirs[i]} \
            {accessions} \
            {clean_reads} \
            {subwkdirs[i]}/bins_concat/bins_concat.fasta \
            {subwkdirs[i]}/SALMON_bins_concat \
            {threads}'''
        subprocess.check_output(cmd, shell=True)

# %% 
'''WARNING: Contig Relative Abundance scripts are task and memory intensive. Each loop requires ~50GB. Avoid parallelizing and multithreading (unless your server can handle it). \

How the relative abundance section is runs:
1. Calculate relative abundance (RRCPM aka redistributed read counts per million) for individual metagenomic assemblies and co-assemlby with respect to mapping method (salmon or minimap) \
2. Calculate relative abundance (RRCPM aka redistributed read counts per million) for MAGs from all MAG methods (MetaWRAP, DAS Tool, MAGScoT) with to mapping method \
3. A contig relative abundace profile tsv file is produce for individual, co-assemlby, and MAGs with respect to MAG method and mapping tool.'''

# %% 
'''Contig Relative Abundance - Co-assembly & Individual Assemblies'''

samples = open(accessions, "r").read().split('\n')
samples.remove('')
    
''''Abundance with Minimap'''
if minimap == True:  

        '''Contig Relative Abundance - Minimap'''            
        '''RRCPM Profile with Minimapped reads'''
        # inputs:
        # ../<abund_dir>/nodeIDs.txt
        # outputs: 
        # ../<abund_dir>/relative_contig_abundance_profiles_crt_<assembly>_minimap.tsv

        minimap_indiv_dirs = [f'{wkdir}/MINIMAP_{sample}' for sample in samples]  
        minimap_coassembly_dir = [f'{wkdir}/MINIMAP_coassembly']
        aln_dirs = minimap_coassembly_dir + minimap_indiv_dirs
        suffixes = [i.split('_')[-1] for i in aln_dirs]
        
        '''NOTE: This step is memory intensive.'''
        for i in range(len(aln_dirs)):
            cmd = f'''{python} {scripts}/30_rrcpm_nonmags_minimap.py \
                --subwkdir {wkdir} \
                --aln_dir {aln_dirs[i]} \
                --abund_dir {wkdir}/CONTIG_RELATIVE_ABUNDANCE \
                --crtdir {wkdir}/CRISPRONE/rNA \
                --out {wkdir}/CONTIG_RELATIVE_ABUNDANCE/relative_contig_abundance_profiles_crt_{suffixes[i]}_minimap.tsv \
                --ref {suffixes[i]}''' 
            subprocess.check_output(cmd, shell=True) 
            
''''Abundance with SALMON'''
if salmon == True :

    '''Contig Relative Abundance - SALMON'''
    '''RRCPM Profile with Minimapped reads'''
    # output: 
    # ../<abund_dir>/relative_contig_abundance_profiles_crt_nonmags_salmon.tsv

    cmd = f'''{python} {scripts}/35_rrcpm_salmon.py \
        --wkdir {wkdir} \
        --subwkdir {wkdir} \
        --abund_dir {wkdir}/CONTIG_RELATIVE_ABUNDANCE \
        --crtdir {wkdir}/CRISPRONE/rNA \
        --out {wkdir}/CONTIG_RELATIVE_ABUNDANCE/relative_contig_abundance_profiles_crt_nonmags_salmon.tsv '''
    subprocess.check_output(cmd, shell=True)



# %% 
'''Contig Relative Abundance - MAGs'''

''''Abundance with Minimap'''
if minimap == True:  

    '''Contig Relative Abundance - Minimap'''            
    '''RRCPM Profile with Minimapped reads'''
    # inputs:
    # ../<abund_dir>/nodeIDs.txt
    # outputs: 
    # ../<abund_dir>/relative_contig_abundance_profiles_crt_minimap.tsv
    for i in range(len(subwkdirs)):
        cmd = f'''{python} {scripts}/31_rrcpm_mags_minimap.py \
            --subwkdir {subwkdirs[i]} \
            --aln_dir {minimap_mag_dirs[i]} \
            --abund_dir {abund_dirs[i]} \
            --crtdir {crisprone_dirs[i]}/rNA \
            --out {abund_dirs[i]}/relative_contig_abundance_profiles_crt_mags_minimap.tsv'''
        subprocess.check_output(cmd, shell=True) 


# %% 
if minimap == True:  

# concatenate abundance profiles (coassembly, individuals, mags)
    aln_dirs = minimap_coassembly_dir + minimap_indiv_dirs
    suffixes = [i.split('_')[-1] for i in aln_dirs]
    for i in range(len(subwkdirs)):
        file1 = f'{abund_dirs[i]}/relative_contig_abundance_profiles_crt_mags_minimap.tsv'
        file2 = f'{abund_dirs[i]}/relative_contig_abundance_profiles_crt_minimap.tsv'
        
        if os.path.exists(file1):
            shutil.copyfile(file1, file2)
            with open(file1, 'r') as f1:
                lines1 = len(f1.readlines()) 
            with open(file2, 'r') as f2:
                lines2 = len(f2.readlines())
            assert lines1 == lines2
        
        lines3 = 0
        for j in suffixes:
            file3 = f'{wkdir}/CONTIG_RELATIVE_ABUNDANCE/relative_contig_abundance_profiles_crt_{j}_minimap.tsv'
            if os.path.exists(file3):
                with open(file3, 'r') as f3:
                    next(f3)
                    with open(file2, 'a') as f2:
                        f2.write(f3.read())          
                with open(file3, 'r') as f3:
                    next(f3)
                    lines3 += len(f3.readlines())
        with open(file1, 'r') as f1:
            next(f1)
            lines1 = len(f1.readlines())           
        with open(file2, 'r') as f2:
            next(f2)
            lines2 = len(f2.readlines())
        assert lines2 == lines1 + lines3
    
    
               
# %%   
''''Abundance with SALMON'''
if salmon == True :

    '''Contig Relative Abundance - SALMON'''
    '''RRCPM Profile with Minimapped reads'''
    # output: 
    # ../<abund_dir>/relative_contig_abundance_profiles_crt_salmon.tsv
    
    # NOTE: script also includes individual & coassemlby rrcpm profiles with mags (unlike minimap script).
    # no need to concatenate indivs & coassembly to mag rrcpm output file.
    for i in range(len(subwkdirs)):
        cmd = f'''{python} {scripts}/35_rrcpm_salmon.py \
            --wkdir {wkdir} \
            --subwkdir {subwkdirs[i]} \
            --abund_dir {abund_dirs[i]} \
            --crtdir {crisprone_dirs[i]}/rNA \
            --out {abund_dirs[i]}/relative_contig_abundance_profiles_crt_salmon.tsv'''
        subprocess.check_output(cmd, shell=True)

# %% 
'''Note: End of Contig Relative Abundance section.
The section below deals with Spacer counts, etc. and is placed here due to dependence on output files for following sections.'''

# %% 
for i in range(len(subwkdirs)):
    
    '''Spacer Counts for Assembly & Assembly Type'''
    # outputs:
    # ../results/assembly2sp_count_sorted.tsv
    # ../results/assembly2spcount_summary.tsv
    
    cmd = f'''{python} {scripts}/36_assembly2sp.py \
        --cdhit_dir {cdhit_dirs[i]} \
        --abund_dir {abund_dirs[i]} \
        --repeat_guides {cdhit_dirs[i]}/top25_representative_repeats.fa \
        --outdir {result_dirs[i]}'''
    subprocess.check_output(cmd, shell=True)
    
    '''Spacer Counts by Assembly and Assembly Type - Rplots'''
    # outputs:
    # ../results/R_figs/sp_count_rNA.pdf/png
    # ../results/R_figs/sp_count_rALL.pdf/png

    topmost_bin_count = 5
    fig_height = 5
    fig_width = 7
    cmd = f'''Rscript {scripts}/plot_scripts/sp_count_plot.R \
        --input {result_dirs[i]}/assembly2sp_summary.tsv \
        --accessions {accessions} \
        --outdir {result_dirs[i]}/R_figs \
        --basename 'sp_count' \
        --topmost_bin_count {topmost_bin_count} '''
    subprocess.check_output(cmd, shell=True)    

# %% 
'''Spacer Counts for Assembly & Assembly Type - All MAG Tools'''
# outputs:
# wkdir/results/assembly2sp_summary_all.tsv
# wkdir/results/assembly2sp_sorted_all.tsv

input_files = ','.join(f'{i}/results/assembly2sp_summary.tsv' for i in subwkdirs)

cmd = f'''{python} {scripts}/37_assembly2sp_all.py \
    --input_files {input_files} \
    --outdir {wkdir}/results \
    --output_file assembly2sp_summary_all.tsv'''
subprocess.check_output(cmd, shell=True)

input_files = ','.join(f'{i}/results/assembly2sp_count_sorted.tsv' for i in subwkdirs)

cmd = f'''{python} {scripts}/37_assembly2sp_all.py \
    --input_files {input_files} \
    --outdir {wkdir}/results \
    --output_file assembly2sp_count_sorted_all.tsv'''
subprocess.check_output(cmd, shell=True)


# %% 
'''Spacer Counts by Assembly and Assembly Type - All MAG Tools - Rplots'''

# outputs:
# ../results/R_figs/sp_count_rNA_all.pdf/png
# ../results/R_figs/sp_count_rALL_all.pdf/png

topmost_bin_count = 5
fig_height = 5
fig_width = 7
cmd = f'''Rscript {scripts}/plot_scripts/sp_count_plot_all.R \
    --input {wkdir}/results/assembly2sp_summary_all.tsv \
    --accessions {accessions} \
    --outdir {wkdir}/results/R_figs \
    --basename 'sp_count' \
    --topmost_bin_count {topmost_bin_count} '''
subprocess.check_output(cmd, shell=True) 


# %% 
'''End of Spacer Count section. 

The section below deals with Assembly Consensus Profiles (awrrcpm of N50 contigs)
and Contig Abundance Profiles (rrcpm).

'''

# %% 
'''Average Weighted Abundance Profile (awrrcpm)'''
# inputs:
# ../<abund_dir>/relative_contig_abundance_profiles_crt_%s.tsv
# ../<abund_dir>/sample_reads_count.txt
# ../<abund_dir>/nodeIDs.txt
# ../<abund_dir>/weighted_relative_contig_abundance_profiles_crt_*.tsv
# outputs:
# ../<abund_dir>/weighted_relative_contig_abundance_profiles_crt_*.tsv
# ../<abund_dir>/A50_concensus_profile_*.tsv
# ../<abund_dir>/contig_abund_*.tsv
# ../<abund_dir>/A50_cosine_distance_*.tsv

for i in range(len(subwkdirs)):
    for j in mapping_tools:
        cmd = f'''{python} {scripts}/38_avg_weight_abun.py  \
            --abund_dir {abund_dirs[i]} \
            --rrcpm_profile {abund_dirs[i]}/relative_contig_abundance_profiles_crt_{j}.tsv  \
            --mapping_tool {j}'''
        subprocess.check_output(cmd, shell=True)
        
# %%
'''Abundance Profile (rrcpm) - Python Plots'''
# inputs:
# ../<abund_dir>/relative_contig_abundance_profiles_crt_%s.tsv
# ../<abund_dir>/A50_concensus_profile_*.tsv
# ../<cdhit_dir>/sp2assembly_<repeat_guide>.out
# ../<abund_dir>/contig_abund_*.tsv
# output : 
# ../results/dist_N50contig_abund_<mapping_tool>.tsv
# ../results/dist_spcontig_abund_<mapping_tool>.tsv
# ../results/python_figs/dist_N50contig_abund_<mapping_tool>.png/pdf
# ../results/python_figs/dist_spcontig_abund_<mapping_tool>.png/pdf
# ../results/python_figs/dist_spcount_<mapping_tool>.png/pdf

topmost_bin_count = 5
repeat_guide = 'rNA'
for i in range(len(subwkdirs)):
    for j in mapping_tools:
        cmd = f'''{python} {scripts}/39_avg_weight_abun_plots.py \
            --accessions {accessions} \
            --topmost_bin_count {topmost_bin_count} \
            --repeat_guide {repeat_guide} \
            --abund_dir {abund_dirs[i]} \
            --cdhit_dir {cdhit_dirs[i]} \
            --outdir {result_dirs[i]} \
            --mapping_tool {j}'''
        subprocess.check_output(cmd, shell=True)

# %%
'''Distribution of Spacer Contig Abundance (rrcpm ) - Rplots'''
# Note: rrcpm values averaged across sample (not contig)
# output : 
# dist_spcontig_abund.pdf/png

for i in range(len(subwkdirs)):
    for j in mapping_tools:    
        cmd = f'''Rscript {scripts}/plot_scripts/dist_spcontig_abund_plot.R\
            --input {result_dirs[i]}/dist_spcontig_abund_{j}.tsv \
            --outdir {result_dirs[i]}/R_figs \
            --basename 'dist_spcontig_abund_{j}' '''
        subprocess.check_output(cmd, shell=True)

# %% 
'''Distribution of N50 Contig Abundance (rrcpm) - Rplots'''
# Note: rrcpm values averaged across sample (not contig)
# output : 
# ../results/R_figs/dist_N50contig_abund.pdf/png

for i in range(len(subwkdirs)):
    for j in mapping_tools:
        cmd = f'''Rscript {scripts}/plot_scripts/dist_N50contig_abund_plot.R\
            --input {result_dirs[i]}/dist_N50contig_abund_{j}.tsv \
            --outdir {result_dirs[i]}/R_figs \
            --basename 'dist_N50contig_abund_{j}' '''
        subprocess.check_output(cmd, shell=True)

# %% 
'''Contig Relative Abundance Profiling (rrcpm)'''
# inputs:
# ../<abund_dir>/relative_contig_abundance_profiles_crt_%s.tsv
# ../<abund_dir>/A50_concensus_profile_*.tsv
# ../<abund_dir>/A50_cosine_distance_*.tsv
# outputs:
# .../results/assembly_rrcpm_pearson_pca_<mapping_tool>.tsv
# .../results/python_figs/<assembly>_rrcpm_pearson_<mapping_tool>.png/pdf
# .../results/assembly_rrcpm_pearson_pca_contig_counts_<mapping_tool>.tsv

topmost_bin_count = 5
for i in range(len(subwkdirs)):
    for j in mapping_tools:
        cmd = f'''{python} {scripts}/40_profilingRRCPM.py \
            --topmost_bin_count {topmost_bin_count} \
            --abund_dir {abund_dirs[i]} \
            --outdir {result_dirs[i]} \
            --mapping_tool {j}'''
        subprocess.check_output(cmd, shell=True)
 
# %% 
'''Contig CRT Outliers Accounting (rrcpm) - All MAG Tools'''
# outputs:
# .../results/assembly_rrcpm_pearson_pca_contig_counts_all.tsv

count = 0
output_file = f'{wkdir}/results/assembly_rrcpm_pearson_pca_contig_counts_all.tsv'
for i in range(len(subwkdirs)):
    mag_tool = subwkdirs[i].split('/')[-1].split('_')[0]
    for j in mapping_tools:
        input_file = f'{subwkdirs[i]}/results/assembly_rrcpm_pearson_pca_contig_counts_{j}.tsv'
        with open(input_file,'r') as fp_in:
            if count == 0:
                with open(output_file, 'w') as fp_out:
                    first_line = fp_in.readline()
                    header = f'mapping_tool\tmag_tool\t{first_line}'
                    fp_out.write(header)
                count += 1
            else:
                next(fp_in)
            with open(output_file, 'a') as fp_out:
                for line in fp_in:
                    new_line = f'{j}\t{mag_tool}\t{line}'
                    fp_out.write(new_line)
                

# %% 
'''Contig Relative Abundance Pearson Correlation - Rplots'''
# outputs:
# ../results/R_figs/assembly_rrcpm_pearson_pca_<mapping_tool>.pdf/png

for i in range(len(subwkdirs)):
    for j in mapping_tools:     
        cmd = f'''Rscript {scripts}/plot_scripts/profilingRRCPM_plot.R \
            --input {result_dirs[i]}/assembly_rrcpm_pearson_pca_{j}.tsv \
            --outdir {result_dirs[i]}/R_figs \
            --basename 'assembly_rrcpm_pearson_pca_{j}' '''
        subprocess.check_output(cmd, shell=True) 

# %% 
'''
End of Abundance Profiling Section (rrcpm, wrrcpm, etc.). 

The section below deals with Spacer Clustering.

'''
# %% 
'''CD-HIT-EST for Spacers - Individuals Only'''
# output: 
# <wkdir>/CD-HIT-EST/concat_<repeat_guide>_spacer.seq
# <wkdir>/CD-HIT-EST/cd-hit-est-<repeat_guide>_spacer.out
# <wkdir>/CD-HIT-EST/cd-hit-est-<repeat_guide>_spacer.out.clstr
cmd = f'''{scripts}/41_cdhit_sp_clustr_indiv.sh \
    {bin}/cdhit-4.8.1/cd-hit-est \
    {wkdir} \
    {wkdir}/CD-HIT-EST \
    {wkdir}/CRISPRONE \
    {accessions} \
    rNA '''
subprocess.check_output(cmd, shell=True)

# %%
'''Pseudo Spacer Presence in Individuals'''
# output: 
# <wkdir>/results/indiv2sp_presence.out

repeat_guide = 'rNA'
cmd = f'''{python} {scripts}/42_indiv2spclust_presence.py \
    --input {wkdir}/CD-HIT-EST/cd-hit-est-{repeat_guide}_spacer.out.clstr \
    --output {wkdir}/results/indiv2spclust_presence.tsv'''
subprocess.check_output(cmd, shell=True)
    

# %% 
'''Unidirectional Heirarchical Clustering of Non-redundant Spacers (Spacer Clusters) across Time (Samples)'''
# output: 
# <wkdir>/results/dendogram_nrspacers_time.png/pdf

cmd = f'''{python} {scripts}/43_dendogram_nrspacers.py \
    --input {wkdir}/results/indiv2spclust_presence_matrix.tsv \
    --n 1000 \
    --outdir {wkdir}/results/python_figs'''
subprocess.check_output(cmd, shell=True)

# %%
'Read Counts & Spacer Counts'
# inputs: 
# <wkdir>/SRR_list.txt
# ../<contig_abund>/sample_reads_count.txt
# <wkdir>/results/assembly2sp_count_sorted_all.tsv
# outputs:
# <wkdir>/python_figs/spacer_counts_by_sample.pdf
# <wkdir>/python_figs/read_counts_by_sample.pdf

cmd = f'''{python} {scripts}/44_spacers_by_sample.py \
    --accessions {accessions} \
    --spcount_file {wkdir}/results/assembly2sp_count_sorted_all.tsv \
    --read_count_file {wkdir}/CONTIG_RELATIVE_ABUNDANCE/sample_reads_count.txt \
    --outdir {wkdir}/results '''
subprocess.check_output(cmd, shell=True)

# %% 
'''The following sections use t-SNE, UMAP, and k-means to cluster 
spacer-containing contigs from a MAG by rrcpm values with respect to sample.
Points on output graphs are contigs and are colored 
by a red gradient for average rrcpm abundance across samples (red=high, white=low)'''

# %% 
'Spacer contig abundance profiles by MAG'
# inputs:
# ../<abund_dir>/relative_contig_abundance_profiles_crt_%s.tsv
# ../results/assembly2sp_count_sorted.tsv
# ../<cdhit_dirs>/concat_rNA_spacer.seq
# outputs:
# ../results/spcontig_matrix_<bin>_<mapping_tool>.tsv

topmost_bin_count = 5
for i in range(len(subwkdirs)):
    for j in mapping_tools:    
        cmd = f'''{python} {scripts}/45_spcontig_matrix.py \
            --accessions {accessions} \
            --spacer_seq_file {cdhit_dirs[i]}/concat_rNA_spacer.seq \
            --spacer_counts_file {result_dirs[i]}/assembly2sp_count_sorted.tsv \
            --abund_file {abund_dirs[i]}/relative_contig_abundance_profiles_crt_{j}.tsv \
            --topmost_bin_count {topmost_bin_count} \
            --basename spcontig_matrix_{j} \
            --outdir {subwkdirs[i]}/results'''
        subprocess.check_output(cmd, shell=True)


        
# %% 
'UMAP & K-means Clustering of spacer-containing contigs by bin'
# inputs:
# ../<abund_dir>/contig_abund_<mapping_tool>.tsv
# ../results/spcontig_matrix_<mapping_tool>_<bin>.tsv
# ../results/assembly2sp_count_sorted.tsv
# ../<cdhit_dirs>/concat_rNA_spacer.seq
# outputs:
# ../results/python_figs/UMAP_kmeans_spcontig_<mapping_tool>_<bin>.pdf/png
    
topmost_bin_count = 5
for i in range(len(subwkdirs)):
    for j in mapping_tools:    
        cmd = f'''{python} {scripts}/46_UMAP_kmeans.py \
        --spacer_seq_file {cdhit_dirs[i]}/concat_rNA_spacer.seq \
        --spacer_counts_file {result_dirs[i]}/assembly2sp_count_sorted.tsv \
        --topmost_bin_count {topmost_bin_count} \
        --avg_abund_file {abund_dirs[i]}/contig_abund_{j}.tsv \
        --abund_file {abund_dirs[i]}/relative_contig_abundance_profiles_crt_{j}.tsv \
        --mapping_tool {j} \
        --outdir {subwkdirs[i]}/results'''
        subprocess.check_output(cmd, shell=True)

# %% 
'tSNE & K-means Clustering of spacer-containing contigs by bin'
# inputs:
# ../<abund_dir>/contig_abund_<mapping_tool>.tsv
# ../results/spcontig_matrix_<mapping_tool>_<bin>.tsv
# ../results/assembly2sp_count_sorted.tsv
# ../<cdhit_dirs>/concat_rNA_spacer.seq
# outputs:
# ../results/python_figs/UMAP_kmeans_spcontig_<mapping_tool>_<bin>.pdf/png
    
topmost_bin_count = 5
for i in range(len(subwkdirs)):
    for j in mapping_tools:    
        cmd = f'''{python} {scripts}/47_tSNE_kmeans.py \
        --spacer_seq_file {cdhit_dirs[i]}/concat_rNA_spacer.seq \
        --spacer_counts_file {result_dirs[i]}/assembly2sp_count_sorted.tsv \
        --topmost_bin_count {topmost_bin_count} \
        --avg_abund_file {abund_dirs[i]}/contig_abund_{j}.tsv \
        --abund_file {abund_dirs[i]}/relative_contig_abundance_profiles_crt_{j}.tsv \
        --mapping_tool {j} \
        --outdir {subwkdirs[i]}/results'''
        subprocess.check_output(cmd, shell=True)


