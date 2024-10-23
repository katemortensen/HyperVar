#!/bin/bash

# Author: Kate Mortensen
# Last Modified: 7/28/2023
# Purpose: discover CRISPR artifacts in assemblies using CRISPRone for co-assembly


python=$1
CRISPRONE=$2
WKDIR=$3
CRISPR=$4

mkdir -p $CRISPR 

# CRISPRone for co-assembly

FNA="$COASSEMBLY/final_assembly.fasta"
ln -s $FNA ${FNA%.fasta}.fa
$python $CRISPRONE \
	--outbase $CRISPR/rNA \
	--outprefix coassembly \
	--fna ${FNA%.fasta}.fa \
	--no_tracrRNA True
wait

unlink ${FNA%.fasta}.fa

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

# list all nodes with crispr from coassembly and write to file

printf "grabbing crispr containing contigs from coassembly \n"
BIN="coassembly"
grep "NODE" $CRISPR/rNA/$BIN/${BIN}.crt > $CRISPR/rNA/$BIN/${BIN}_crt_nodes.txt
sed 's/SEQ: //' $CRISPR/rNA/$BIN/${BIN}_crt_nodes.txt | \
    awk -F "[ ]+" '{print $1}' \
    > $CRISPR/rNA/$BIN/tmp \
    && mv $CRISPR/rNA/$BIN/tmp $CRISPR/rNA/$BIN/${BIN}_crt_nodes.txt

