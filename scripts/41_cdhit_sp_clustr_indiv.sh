#!/bin/bash

# Author: Kate Mortensen
# Last Modified: 8/29/2023
# Purpose: collect spacers discovered by CRISPRone and cluster using CD-HIT-EST

cd_hit_est=$1
WKDIR=$2
CDHIT=$3
CRISPR=$4
ACCESSIONS=$5
RTAG=$6

# cd_hit_est=/home/kmorten/Hypervariable-region-aware-co-assembly-of-metagenomes/hypervar-pipeline/bin/cdhit-4.8.1/cd-hit-est
# WKDIR=/home/kmorten/ADMB
# CDHIT=$WKDIR/CD-HIT-EST
# CRISPR=$WKDIR/CRISPRONE
# ACCESSIONS=$WKDIR/SRR_list.txt
# RTAG=rNA

mkdir -p $CDHIT
rm -f $CDHIT/concat_${RTAG}_spacer.seq

echo "concatenate all spacer files together"

rm -f $CDHIT/concat_${RTAG}_spacer.seq
A=$(cat $ACCESSIONS)
for ASSEMBLY in ${A[@]}; do
    TARGET_ASSEMBLY="$CRISPR/$RTAG/$ASSEMBLY"
    cat $TARGET_ASSEMBLY/$ASSEMBLY-spacer.seq >> $CDHIT/concat_${RTAG}_spacer.seq
done


echo "runnning cd-hit-est on concated repeat-associated spacers to get clusters of repeats"

$cd_hit_est -i $CDHIT/concat_${RTAG}_spacer.seq \
-sc 1 -sf -1 -d 200 -o $CDHIT/cd-hit-est-${RTAG}_spacer.out 
    




