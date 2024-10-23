#!/bin/bash

# Author: Kate Mortensen
# Last Modified: 8/29/2023
# Purpose: collect spacers discovered by CRISPRone and cluster using CD-HIT-EST

cd_hit_est=$1
WKDIR=$2
CDHIT=$3
CRISPR=$4

mkdir -p $CDHIT
chmod -R 775 $WKDIR

R=$(ls -d $CRISPR/*/)

echo "concatenate all spacer files together"

for REP in ${R[@]}; do
     RTAG=${REP##$CRISPR/}
     RTAG=${RTAG%/}
     A=$(ls $CRISPR/$RTAG)
     rm -f $CDHIT/concat_${RTAG}_spacer.seq
     for ASSEMBLY in ${A[@]}; do
	    TARGET_ASSEMBLY="$CRISPR/$RTAG/$ASSEMBLY"
	    cat $TARGET_ASSEMBLY/$ASSEMBLY-spacer.seq >> $CDHIT/concat_${RTAG}_spacer.seq
     done
done

echo "runnning cd-hit-est on concated repeat-associated spacers to get clusters of repeats"
for REP in ${R[@]}; do
    RTAG=${REP##$CRISPR/}
    RTAG=${RTAG%/}
    
    $cd_hit_est -i $CDHIT/concat_${RTAG}_spacer.seq \
	-sc 1 -sf -1 -d 200 -o $CDHIT/cd-hit-est-${RTAG}_spacer.out 
    
done




