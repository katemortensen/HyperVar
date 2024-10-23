#!/bin/bash

# Author: Kate Mortensen
# Last Modified: 10/21/2024
# Purpose: MAGScoT binning of co-assembly contigs 

# setup variables
# . ./00_vars.sh

mkdir -p $MAGSCOT_MAGS
       	
# SINGLE COPY MARKER INDENTIFICAITON (GTDBtk rel 207)

BASENAME='final_assembly'
gzip $COASSEMBLY/final_assembly.fasta
cd $MAGSCOT_MAGS

### ORF detection with prodigal

#zcat $COASSEMBLY/$BASENAME.fasta.gz | prodigal -p meta -a $BASENAME.prodigal.faa -d $BASENAME.prodigal.ffn -o tmpfile

### ALTERNATIVE: Fast parallel ORF detection with prodigal
mkdir -p tmp_workfolder
zcat $COASSEMBLY/$BASENAME.fasta.gz | parallel -j $THREADS --block 999k --recstart '>' --pipe prodigal -p meta -a tmp_workfolder/$BASENAME.{#}.faa -d tmp_workfolder/$BASENAME.{#}.ffn -o tmpfile
cat tmp_workfolder/$BASENAME.*.faa > $BASENAME.prodigal.faa
cat tmp_workfolder/$BASENAME.*.ffn > $BASENAME.prodigal.ffn
rm -r tmp_workfolder tmpfile

### annotation of protein sequences using HMMer and GTDBtk r207 marker genes
hmmsearch -o $BASENAME.hmm.tigr.out --tblout $BASENAME.hmm.tigr.hit.out --noali --notextw --cut_nc --cpu $THREADS $MAGScoT/hmm/gtdbtk_rel207_tigrfam.hmm $BASENAME.prodigal.faa
hmmsearch -o $BASENAME.hmm.pfam.out --tblout $BASENAME.hmm.pfam.hit.out --noali --notextw --cut_nc --cpu $THREADS $MAGScoT/hmm/gtdbtk_rel207_Pfam-A.hmm $BASENAME.prodigal.faa

cat $BASENAME.hmm.tigr.hit.out | grep -v "^#" | awk '{print $1"\t"$3"\t"$5}' > $BASENAME.tigr
cat $BASENAME.hmm.pfam.hit.out | grep -v "^#" | awk '{print $1"\t"$4"\t"$5}' > $BASENAME.pfam
cat $BASENAME.pfam $BASENAME.tigr > $BASENAME.hmm


# contigs2bins input file from initial binning

PATHDIRS=$(ls -d $INITIAL_BINNING/*_bins)
FILE_OUT=$MAGSCOT_MAGS/$BASENAME.contigs_to_bin.tsv
rm -f $FILE_OUT
#echo -e "Bin ID\tContig ID\tBinning Tool" >> $FILE_OUT
for PATHDIR in ${PATHDIRS[@]}; do
    DIR=${PATHDIR##*/}
    METHOD=${DIR%_}
    PBINS=$(ls $INITIAL_BINNING/$METHOD/bin*.fa)
    for PBIN in ${PBINS[@]}; do
	FBIN=${PBIN##*/}
	BIN=${FBIN%.fa}
	FILE_IN=$PBIN
	lines=($(grep "NODE" $FILE_IN))
	for line in ${lines[@]}; do 
	    echo -e "${BASENAME}_${METHOD%_bins}_$BIN\t${line##*>}\t${METHOD%_bins}" >> $FILE_OUT
	done
    done
done    

# MAGScoT SCORING AND BIN REFINEMENT  
 
Rscript $MAGScoT/MAGScoT.R -i $FILE_OUT --hmm $BASENAME.hmm
 
# compile MAGs 

python3 $CURRENT/scripts/call_scripts/magscot_compile_mags.py \
    -f $COASSEMBLY/$BASENAME.fasta \
    -i $MAGSCOT_MAGS/MAGScoT.refined.contig_to_bin.out \
    -o $MAGSCOT_MAGS/bins 


# fix magscot bin names for standardization and compatability with sp2assembly.py script

B=$(ls $MAGSCOT_MAGS/bins)
for BIN in ${B[@]}; do
    NEW_BIN=$(echo "$BIN" | sed -r 's/[_]+/./g')
    echo "old name: $BIN"
    echo "new name $NEW_BIN:"
    mv $MAGSCOT_MAGS/bins/$BIN $MAGSCOT_MAGS/bins/$NEW_BIN
done




cd $CURRENT/scripts


