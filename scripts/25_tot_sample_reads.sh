#!/usr/bin/bash

# Author: Kate Mortensen
# Last Modified: 8/14/2023
# Purpose: count reads from accessions 


ACCESSIONS=$1
CLEAN=$2
CONTIG_RELATIVE_ABUNDANCE=$3

echo "$CONTIG_RELATIVE_ABUNDANCE"

CLEAN_READS_COUNT=$CONTIG_RELATIVE_ABUNDANCE/sample_reads_count.txt
echo "$CLEAN_READS_COUNT"

mkdir -p $CONTIG_RELATIVE_ABUNDANCE

if [ -s $CLEAN_READS_COUNT ]; then
    rm -f $CLEAN_READS_COUNT
fi
touch $CLEAN_READS_COUNT
printf '%s\t%s\n' "sample" "total_reads" >> $CLEAN_READS_COUNT

echo "counting reads from individual accessions"
S=$(cat $ACCESSIONS)
for SRR in ${S[@]}; do
    FORWARD=$CLEAN/$SRR/final_pure_reads_1.fastq
    REVERSE=$CLEAN/$SRR/final_pure_reads_2.fastq
    echo "$FORWARD"
    echo "$REVERSE"

    FWD_COUNT=$(grep -c "@" $FORWARD) 
    REV_COUNT=$(grep -c "@" $REVERSE)
    wait 
    TOT_COUNT=$(($FWD_COUNT + $REV_COUNT))
    printf '%s\t%s\n' "$SRR" "$TOT_COUNT" >> $CLEAN_READS_COUNT
done


