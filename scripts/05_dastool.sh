#!/bin/bash

# Author: Kate Mortensen
# Last Modified: 10/21/2024
# Purpose: DAS Tool binning of co-assembly contigs 

# setup variables
# . ./00_vars.sh

mkdir -p $DASTOOL_MAGS

# prepare contigs2bins input file from initial binning done through metawrap

PATHDIRS=$(ls -d $INITIAL_BINNING/*_bins)
for PATHDIR in ${PATHDIRS[@]}; do
    DIR=${PATHDIR##*/}
    METHOD=${DIR%_}
    FILE_OUT=$DASTOOL_MAGS/${METHOD}_contigs2bins.tsv
    rm -f $FILE_OUT
    BINS=$(ls $INITIAL_BINNING/$METHOD/)
    for BIN in ${BINS[@]}; do
	FILE_IN=$INITIAL_BINNING/$METHOD/$BIN
	lines=($(grep ">" $FILE_IN))
	for line in ${lines[@]}; do 
	    echo -e "${line##*>}\t${METHOD%_bins}${BIN%.fa}" >> $FILE_OUT
	done
    done
done    

# contigs2bins (exclude unbinned)
PATHDIRS=$(ls -d $INITIAL_BINNING/*_bins)
for PATHDIR in ${PATHDIRS[@]}; do
    DIR=${PATHDIR##*/}
    METHOD=${DIR%_}
    FILE_OUT=$DASTOOL_MAGS/${METHOD}_contigs2bins.tsv
    rm -f $FILE_OUT
    BINS=$(ls $INITIAL_BINNING/$METHOD/)
    for BIN in ${BINS[@]}; do
	if [ -z "$(echo $BIN | grep -o "unbin")" ]; then
            FILE_IN=$INITIAL_BINNING/$METHOD/$BIN
            lines=($(grep ">" $FILE_IN))
            for line in ${lines[@]}; do
                echo -e "${line##*>}\t${METHOD%_bins}${BIN%.fa}" >> $FILE_OUT
            done
	else
	    echo "skip $BIN ("bin" of unbinned contigs)"
	fi
    done
done



# DAS Tool

cd $DASTOOL_MAGS

contig2bin_inputs=$(find $DASTOOL_MAGS -path '*contigs2bins.tsv' | tr '\n' ',')
OUT_PREFIX="dastool_out"

DAS_Tool -i ${contig2bin_inputs%,} \
    -c $COASSEMBLY/final_assembly.fasta \
    -o $OUT_PREFIX \
    --write_bin_evals \
    --write_bins \
    --write_unbinned \
    -t $THREADS 

cd $CURRENT/scripts

# fix bin naming scheme for standardization and compatibility with sp2assembly.py script

mv $DASTOOL_MAGS/dastool_out_DASTool_bins $DASTOOL_MAGS/bins

B=$(ls $DASTOOL_MAGS/bins)
for BIN in ${B[@]}; do
    NEW_BIN=$(echo "$BIN" | sed -r 's/[_]+/./g')
    OLD_PATH=$DASTOOL_MAGS/bins/$BIN
    NEW_PATH=$DASTOOL_MAGS/bins/$NEW_BIN
    echo "old: $OLD_PATH"
    echo "new: $NEW_PATH"
    mv $OLD_PATH $NEW_PATH
done


