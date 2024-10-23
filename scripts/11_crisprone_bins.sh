#!/bin/bash

# Author: Kate Mortensen
# Last Modified: 7/28/2023
# Purpose: discover CRISPR artifacts in assemblies using CRISPRone


python=$1
CRISPRONE=$2
COASSEMBLY_INDIV_CRISPR=$3
SUBWKDIR=$4
CRISPR=$5
BINS=$6
THREADS=$7


mkdir -p $CRISPR
mkdir -p $CRISPR/rNA
mkdir -p $SUBWKDIR/commands
CMD=$SUBWKDIR/commands/crisprone_cmds.txt
rm -f $CMD

# copy over CRISPRone for Coassembly and Individual Assemblies

cp -r $COASSEMBLY_INDIV_CRISPR/rNA/* $CRISPR/rNA/.

# CRISPRone for MAGs

B=$(ls $BINS)
for BIN in ${B[@]}; do

    FNA="$BINS/$BIN"

#    $python $CRISPRONE \
#        --outbase $CRISPR/rNA \
#        --outprefix ${BIN%.fa} \
#        --fna $FNA \
#        --no_tracrRNA True

    echo "$python $CRISPRONE --outbase $CRISPR/rNA --outprefix ${BIN%.fa} --fna $FNA --no_tracrRNA True" >> $CMD

done

parallel -j $THREADS < $CMD 

# find all assemblies that contain CRISPR artifacts

NOTEMPTY="$CRISPR/rNA_assemblies_with_crispr.txt"
rm -f $NOTEMPTY
A=$(ls -d $CRISPR/rNA/*)
for ASSEMBLY in ${A[@]}; do

    if [ -s $ASSEMBLY/*.crt ]; then
	echo "${ASSEMBLY##*/}" >> $NOTEMPTY
    fi

done

# list all nodes with crispr from mags and write to file 

B=$(ls $BINS/)
for BIN in ${B[@]}; do
    BIN=${BIN%.fa}
    grep "NODE" $CRISPR/rNA/$BIN/$BIN.crt \
        > $CRISPR/rNA/$BIN/${BIN}_crt_nodes.txt

    sed 's/SEQ: //' $CRISPR/rNA/$BIN/${BIN}_crt_nodes.txt | \
    awk -F "[ ]+" '{print $1}' \
    > $CRISPR/rNA/$BIN/tmp \
    && mv $CRISPR/rNA/$BIN/tmp \
    $CRISPR/rNA/$BIN/${BIN}_crt_nodes.txt
done




