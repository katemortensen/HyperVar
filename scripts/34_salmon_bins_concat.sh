#!/bin/bash

# Author: Kate Mortensen
# Last Modified: 7/19/2023
# Purpose: execute SALMON for quantification

# setup variables
# . ./00_vars.sh

SALMON=$1
SUBWKDIR=$2
ACCESSIONS=$3
CLEAN=$4
BINS_CONCAT=$5
OUTDIR=$6
THREADS=$7

export PATH="$SALMON:$PATH"

printf "software paths: \n"
which salmon

printf "symlink cleaned reads temporarily for salmon mapping \n"

pushd $SUBWKDIR
S=$(cat $ACCESSIONS)
for SRR in ${S[@]}; do
    ln -s $CLEAN/$SRR/*_1.fastq $SUBWKDIR/${SRR}_1.fastq
    ln -s $CLEAN/$SRR/*_2.fastq $SUBWKDIR/${SRR}_2.fastq
done

printf "salmon quantification of concatenated mags \n"

salmon index \
    -p $THREADS \
    -t $BINS_CONCAT \
    -i $OUTDIR/assembly_index

S=$(cat $ACCESSIONS)
for SRR in ${S[@]}; do
    salmon quant \
	-i $OUTDIR/assembly_index \
	--libType IU \
	-1 ${SRR}_1.fastq \
	-2 ${SRR}_2.fastq \
	-o $OUTDIR/alignment_files/$SRR.quant \
	--meta \
	-p $THREADS
done

printf "unlink cleaned fastas \n"

S=$(cat $ACCESSIONS)
for SRR in ${S[@]}; do
    unlink $SUBWKDIR/${SRR}_1.fastq
    unlink $SUBWKDIR/${SRR}_2.fastq
done

popd



