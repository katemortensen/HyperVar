#!/bin/bash

# Author: Kate Mortensen
# Last Modified: 8/29/2023
# Purpose: cluster assembly repeats using CD-HIT-EST 


CDHIT=$1
CRISPR=$2
cd_hit_est=$3

mkdir -p $CDHIT

# concatenate all repeat consensus file together

rm -f "$CDHIT/concat_repeatcon.seq"
    
A=$(ls -d $CRISPR/rNA/*)
for ASSEMBLY in ${A[@]}; do
    cat $ASSEMBLY/repeatcon.seq >> $CDHIT/concat_repeatcon.seq
done

# run cd-hit-est on concated repeat consensus to get clusters of repeats

$cd_hit_est \
    -i $CDHIT/concat_repeatcon.seq \
    -sc 1 -sf -1 -d 200 -o $CDHIT/cd-hit-est.out


