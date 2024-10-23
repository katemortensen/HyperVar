#!/bin/bash

# Author: Kate Mortensen
# Last Modified: 8/29/2023
# Purpose: run crisprone-summarize script on each assembly and find repeats

python=$1
SUMMARIZE=$2
SUBWKDIR=$3
CRISPR=$4
THREADS=$5

# loop

#A=$(ls -d $CRISPR/rNA/*)
#for ASSEMBLY in ${A[@]}; do
#    $python $SUMMARIZE \
#	-f $ASSEMBLY/*.crt \
#	-repeat $ASSEMBLY/repeat.seq \
#	-consensus $ASSEMBLY/repeatcon.seq 
#done

# in parallel 

mkdir -p $SUBWKDIR/commands
CMD=$SUBWKDIR/commands/crisprone_summary_cmds.txt
rm -f $CMD

A=$(ls -d $CRISPR/rNA/*)
for ASSEMBLY in ${A[@]}; do
    echo "$python $SUMMARIZE -f $ASSEMBLY/*.crt -repeat $ASSEMBLY/repeat.seq -consensus $ASSEMBLY/repeatcon.seq" >> $CMD
 
done

parallel -j $THREADS < $CMD

