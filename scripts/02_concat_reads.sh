#!/bin/bash

# Author: Kate Mortensen
# Last Modified: 10/21/2023
# Purpose: concatenated cleaned reads into forward and reverse fastq files

# setup variables
# . 00_vars.sh

WKDIR=$1
ACCESSIONS=$2

touch $WKDIR/READS_1.fastq
touch $WKDIR/READS_2.fastq

S=$(cat $ACCESSIONS)

for SRR in ${S[@]}; do

    R1="$CLEAN/$SRR/final_pure_reads_1.fastq"
    R2="$CLEAN/$SRR/final_pure_reads_2.fastq"

    R1_COUNTS="$(cat $R1|wc -l)/4|bc"
    R2_COUNTS="$(cat $R2|wc -l)/4|bc"

    if [ "$R1_COUNTS" == "$R2_COUNTS" ]; then

        cat $R1 >> $WKDIR/READS_1.fastq
        wait
        cat $R2 >> $WKDIR/READS_2.fastq
        wait

    else
        echo "$SRR unequeal read counts or already logged" 
    fi
done

