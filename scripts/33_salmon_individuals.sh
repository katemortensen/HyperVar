#!/bin/bash

# Author: Kate Mortensen
# Last Modified: 8/8/2023
# Purpose: execute SALMON for quantification (coassembly)

SALMON=$1
WKDIR=$2
ACCESSIONS=$3
CLEAN=$4
INDIVIDUALS=$5
OUTDIR=$6
THREADS=$7

printf "software paths: \n"
which salmon

printf "symlink cleaned reads temporarily for salmon mapping \n"

cd $WKDIR
S=$(cat $ACCESSIONS)
for SRR in ${S[@]}; do
    ln -s $CLEAN/$SRR/*_1.fastq $WKDIR/${SRR}_1.fastq
    ln -s $CLEAN/$SRR/*_2.fastq $WKDIR/${SRR}_2.fastq
done


printf "salmon quantification of individual assemblies \n"

for SRR in ${S[@]}; do
    salmon index \
	-p $THREADS \
	-t $INDIVIDUALS/$SRR/ASSEMBLY/final_assembly.fasta \
	-i $OUTDIR/SALMON_${SRR}/assembly_index
done

S=$(cat $ACCESSIONS)
for REF in ${S[@]}; do
    for SRR in ${S[@]}; do
	salmon quant \
	    -i $OUTDIR/SALMON_${REF}/assembly_index \
	    --libType IU \
	    -1 ${SRR}_1.fastq \
	    -2 ${SRR}_2.fastq \
	    -o $OUTDIR/SALMON_${REF}/alignment_files/$SRR.quant \
	    --meta \
	    -p $THREADS
done

printf "unlink cleaned fastas \n"

S=$(cat $ACCESSIONS)
for SRR in ${S[@]}; do
    unlink $WKDIR/${SRR}_1.fastq
    unlink $WKDIR/${SRR}_2.fastq
done


popd
