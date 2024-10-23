#!/bin/bash

# Author: Kate Mortensen
# Last Modified: 9/4/2023
# Purpose: create co-assembly


WKDIR=$1
THREADS=$2
#cd $WKDIR

mkdir -p $COASSEMBLY

R1_COUNTS="$(cat $WKDIR/READS_1.fastq|wc -l)/4|bc"
R2_COUNTS="$(cat $WKDIR/READS_2.fastq|wc -l)/4|bc"

if [ $R1_COUNTS == $R2_COUNTS ]; then
    metawrap assembly \
	-1 $WKDIR/READS_1.fastq \
	-2 $WKDIR/READS_2.fastq \
	-m $MEM -t $THREADS \
	--metaspades -o $COASSEMBLY 
    wait 
else
    echo "read files have unequal lengths" 
fi


S=$(ls $CLEAN)
for SRR in ${S[@]}; do
    mkdir -p $INDIV/$SRR

    R1_COUNTS="$(cat $CLEAN/$SRR/*_1.fastq|wc -l)/4|bc"
    R2_COUNTS="$(cat $CLEAN/$SRR/*_2.fastq|wc -l)/4|bc"

    if [ $R1_COUNTS == $R2_COUNTS ]; then
        metawrap assembly \
        -1 $CLEAN/$SRR/*_1.fastq \
        -2 $CLEAN/$SRR/*_2.fastq \
        -m $MEM -t $THREADS \
        --metaspades -o $INDIV/$SRR/ASSEMBLY
    else
	echo "read files for indiv accessions have unequal lengths"
    fi 
done

