#!/usr/bin/bash

# Author: Kate Mortensen
# Last Modified: 8-2-2023
# Purpose: find all clean node IDs from reassembled bins and write to file - to be used in relative contig abundance profile estimation


WKDIR=$1
SUBWKDIR=$2
COASSEMBLY=$3
INDIV=$4
BINS=$5
ACCESSIONS=$6 

CONTIG_RELATIVE_ABUNDANCE="$SUBWKDIR/CONTIG_RELATIVE_ABUNDANCE"

printf "create nodeIDs.txt \n"

NODEIDS=$CONTIG_RELATIVE_ABUNDANCE/nodeIDs.txt
if [ -s $NODEIDS ];then
    rm -f $NODEIDS
fi
touch $NODEIDS
printf '%s\t%s\t%s\n' "assembly" "contig" "assembly_type" >> $NODEIDS


##################################
#    COASSEMBLY 
##################################

printf "collecting node IDs from metagenomic coassembly \n"
N=$(grep ">"  $COASSEMBLY/final_assembly.fasta | sed 's/>//')
for NODE in ${N[@]}; do
    printf '%s\t%s\t%s\n' "coassembly" "$NODE" "coassembly" >> $NODEIDS
done

##################################
#    INDIVDUALS 
##################################

printf "collecting node IDs for individual assemblies \n"
S=$(cat $ACCESSIONS)
for SRR in ${S[@]}; do
    N=$(grep ">" $INDIV/$SRR/ASSEMBLY/final_assembly.fasta | sed 's/>//')
    for NODE in ${N[@]}; do
        printf '%s\t%s\t%s\n' "$SRR" "$NODE" "individual" >> $NODEIDS
    done
done

##################################
#    BINS (MAGs) 
##################################

printf "collecting node IDs from reassembled bins \n"
B=$(ls $BINS)
for BIN in ${B[@]}; do
    N=$(grep ">" $BINS/$BIN | sed 's/>//')
    for NODE in ${N[@]}; do 
        printf '%s\t%s\t%s\n' "${BIN%.fa}" "$NODE" "mag" >> $NODEIDS
    done
done


