#!/bin/bash

# Author: Kate Mortensen
# Last Modified: 08/29/2023
# Purpose: Spacer discovery using CRISPRone software

python=$1
SUMMARIZE=$2
CRISPR=$3

R=$(ls -d $CRISPR/*/) # include rNA dir

for REP in ${R[@]}; do
     RTAG=${REP##$CRISPR/}
     RTAG=${RTAG%/}
     A=$(ls $CRISPR/$RTAG)

     # run summary script on each bin to find spacers associated with guide repeat
     for ASSEMBLY in ${A[@]}; do

	 TARGET_DIR="$CRISPR/$RTAG/$ASSEMBLY"

	 rm -f $TARGET_DIR/$ASSEMBLY-spacer.seq

	 echo "REPEAT: ${RTAG}"
	 echo "ASSEMBLY: ${ASSEMBLY}"

	 python $SUMMARIZE \
		-f $TARGET_DIR/$ASSEMBLY.crt \
		-spacer $TARGET_DIR/$ASSEMBLY-spacer.seq \
		-maxmis 2 
     done
done 

# change spacer fasta headers to match bin name
for REP in ${R[@]}; do
	RTAG=${REP##$CRISPR/}
    RTAG=${RTAG%/}
    A=$(ls $CRISPR/$RTAG)
    
    for ASSEMBLY in ${A[@]}; do
		TARGET_DIR="$CRISPR/$RTAG/$ASSEMBLY"

		echo "REPEAT: ${RTAG}"
		echo "ASSEMBLY: ${ASSEMBLY}"
		sed -i "s/>/>${ASSEMBLY}_/" $TARGET_DIR/$ASSEMBLY-spacer.seq

    done
done

wait



