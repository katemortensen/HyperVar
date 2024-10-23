#!/bin/bash

# Author: Kate Mortensen
# Last Modified: 10/21/2024
# Purpose: execute initial binning of contigs in co-assembly
# Purpose: execute metaWRAP's bin refinement module 
# Purpose: execute MetaWRAP's bin reassembly module to re-assemble refined bins into MAGs

# setup variables
# . 00_vars.sh

# conda env
# module load anaconda
# source activate metawrap-env

ACCESSIONS=$1
SUBWKDIR=$2
THREADS=$3

mkdir -p $INITIAL_BINNING

metawrap binning \
    -o $INITIAL_BINNING \
    -t $THREADS \
    -a $COASSEMBLY/final_assembly.fasta \
    --metabat2 --maxbin2 --concoct $SUBWKDIR/*fastq 

wait 

mkdir -p $BIN_REFINEMENT

metawrap bin_refinement \
    -o $BIN_REFINEMENT \
    -A $INITIAL_BINNING/metabat2_bins/ \
    -B $INITIAL_BINNING/maxbin2_bins/ \
    -C $INITIAL_BINNING/concoct_bins/ \
    -t $THREADS -c 50 -x 5 # 50% completion/5% contamination

mkdir -p $BIN_REASSEMBLY

S=$(cat $ACCESSIONS)
for SRR in ${S[@]}; do
    ln -s $CLEAN/$SRR/*_1.fastq $SUBWKDIR/${SRR}_1.fastq
    ln -s $CLEAN/$SRR/*_2.fastq $SUBWKDIR/${SRR}_2.fastq
done

metawrap reassemble_bins \
    -o $BIN_REASSEMBLY \
    -1 $SUBWKDIR/*_1.fastq \
    -2 $SUBWKDIR/*_2.fastq \
    -t $THREADS -m $MEM \
    -c 50 -x 5 \ # 50% completion/5% contamination
    -b $BIN_REFINEMENT/metawrap_*_bins

wait

# B=$(ls /home/kmorten/ADMB/BIN_REASSEMBLY/reassembled_bins)
# for BIN in ${B[@]}; do
#     echo "${BIN%.fa}" >> $SUBWKDIR/bins.txt
# done

S=$(cat $ACCESSIONS)
for SRR in ${S[@]}; do
    unlink $SUBWKDIR/${SRR}_1.fastq
    unlink $SUBWKDIR/${SRR}_2.fastq
done





