#!/usr/bin/bash

# Author: Kate Mortensen
# Last Modified: 8-2-2023
# Purpose: create nodeID conversion key

printf "rename contigs in concatenated mags (not an issue for coassembly) \n"
printf "keep a record of bin names, old contig names, and new contig names \n"


SUBWKDIR=$1
BINS=$2

CONTIG_RELATIVE_ABUNDANCE="$SUBWKDIR/CONTIG_RELATIVE_ABUNDANCE"

NODEIDS=$CONTIG_RELATIVE_ABUNDANCE/nodeIDs_bin_conversion.txt
if [ -s $NODEIDS ];then
    rm -f $NODEIDS
fi
touch $NODEIDS
printf '%s\t%s\t%s\n' "bin" "old_contig_name" "new_contig_name" >> $NODEIDS

printf "creating key for contig names from reassembled MAGs \n"
B=$(ls $BINS)
COUNTER=0
for BIN in ${B[@]}; do
        N=$(grep ">" $BINS/$BIN | sed 's/>//')
    for NODE in ${N[@]}; do
	let COUNTER++
	NEW_NAME="NODE_${COUNTER}_length${NODE##*length}"
        printf '%s\t%s\t%s\n' "${BIN%.fa}" "$NODE" "$NEW_NAME" >> $NODEIDS
    done
done

