#!/bin/bash

# Author: Kate Mortensen
# Last Modified: 7/28/2023
# Purpose: discover CRISPR artifacts in assemblies using CRISPRone

python=$1
CRISPRONE=$2
WKDIR=$3
CRISPR=$4
INDIVIDUALS=$5

mkdir -p $CRISPR 

# CRISPRone for individual assemblies

S=$(ls $INDIVIDUALS)
for SRR in ${S[@]}; do

    FNA="$INDIVIDUALS/$SRR/ASSEMBLY/final_assembly.fasta"
    ln -s $FNA ${FNA%.fasta}.fa

    $python $CRISPRONE \
        --outbase $CRISPR/rNA \
        --outprefix $SRR \
        --fna ${FNA%.fasta}.fa \
	--no_tracrRNA True

    unlink ${FNA%.fasta}.fa
done

# find all assemblies that contain CRISPR artifacts

NOTEMPTY="$CRISPR/rNA_assemblies_with_crispr.txt"
rm -f $NOTEMPTY
A=$(ls -d $CRISPR/rNA/*)
for ASSEMBLY in ${A[@]}; do

    if [ -s $ASSEMBLY/*.crt ]; then
	echo "${ASSEMBLY##*/}" >> $NOTEMPTY
    fi

done

wait

# list all nodes with crispr from individual assemblies and write to file 

printf "grabbing crispr containing contigs form individual assemblies \n"
S=$(cat $WKDIR/SRR_list.txt)

for SRR in ${S[@]}; do

    grep "NODE" $CRISPR/rNA/$SRR/$SRR.crt > $CRISPR/rNA/$SRR/${SRR}_crt_nodes.txt

    sed 's/SEQ: //' $CRISPR/rNA/$SRR/${SRR}_crt_nodes.txt | \
    awk -F "[ ]+" '{print $1}' \
    > $CRISPR/rNA/$SRR/tmp \
    && mv $CRISPR/rNA/$SRR/tmp $CRISPR/rNA/$SRR/${SRR}_crt_nodes.txt

done


