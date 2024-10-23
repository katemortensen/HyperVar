#!/bin/bash

# Author: Kate Mortensen
# Last Modified: 5/15/23
# Purpose: clean reads

# setup variables
# . 00_vars.sh 

# conda env
# module load anaconda
# source activate metawrap-env

ACCESSION=$1
RAW=$2


mkdir -p $CLEAN

# bash function for cleaning reads
clean_reads(){

    mkdir -p $CLEAN/$SRR
    
    R1="$RAW/$SRR/${SRR}_1.fastq"
    R2="$RAW/$SRR/${SRR}_2.fastq"
    gunzip ${R1}.gz
    gunzip ${R2}.gz

    metawrap read_qc --skip-bmtagger -1 \
    $R1 -2 $R2 -t $THREADS -o $CLEAN/$SRR 

}

S=$(cat $ACCESSIONS)
for SRR in ${S[@]}; do
    # clean_reads $SRR $RAW $CLEAN &
    clean_reads $SRR $RAW $CLEAN 
done


