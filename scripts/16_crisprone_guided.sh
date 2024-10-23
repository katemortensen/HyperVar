#!/bin/bash

# Author: Kate Mortensen
# Last Modified: 10/21/2023
# Purpose: CRISPRone using representative repeats from clustered repeats as guides for CRISPR discovery 

python=$1
call_script=$2
CRISPRONE=$3
MAG_DIR=$4
REPRESENTATIVE_REPEATS=$5
CDHIT=$6
CRISPR=$7
WKDIR=$8
COASSEMBLY=$9

# python="python3"
# call_script='/home/kmorten/Hypervariable-region-aware-co-assembly-of-metagenomes/hypervar-pipeline/scripts/call_scripts/extract_contig.py'
# # CRISPRONE='/home/kmorten/Hypervariable-region-aware-co-assembly-of-metagenomes/hypervar-pipeline/bin/CRISPRone/crisprone-local-2022.py'
# CRISPRONE='/var/www/omics.informatics.indiana.edu/CRISPRone/crisprone-local-2022.py'
# MAG_DIR='/home/kmorten/ADMB//DASTOOL_MAGS/bins'
# REPRESENTATIVE_REPEATS='/home/kmorten/ADMB//DASTOOL_MAGS/CD-HIT-EST/top25_representative_repeats.fa'
# CDHIT='/home/kmorten/ADMB//DASTOOL_MAGS/CD-HIT-EST'
# CRISPR='/home/kmorten/ADMB//DASTOOL_MAGS/CRISPRONE'
# WKDIR='/home/kmorten/ADMB/'
# COASSEMBLY='/home/kmorten/ADMB/COASSEMBLY/final_assembly.fasta'

#####################################################################
# MAGs - repeat guides r0 to r25 (1 at a time )
#####################################################################

B=$(ls $MAG_DIR) # to loop over bins
R="$(grep ">" $REPRESENTATIVE_REPEATS)" # to loop over repeats

echo "Running CRISPRone (with repeat guides) for MAGs"
for REP in ${R[@]}; do
    RTAG=${REP##*>}
    #echo "Representative Repeat: $RTAG"
    # extracts repeat fasta from file
    # use python script to extract repeat contigs 
    if [ ! -f $CHDIT/$RTAG.fa ]; then
        $python $call_script \
            -n $RTAG -f $REPRESENTATIVE_REPEATS -o $CDHIT
    fi
    NOTEMPTY="$CRISPR/${RTAG}_assemblies_with_crispr.txt"
    if [ -s $NOTEMPTY ]; then 
	    rm -f $NOTEMPTY
    fi

    mkdir -p $CRISPR/$RTAG

    for BIN in ${B[@]}; do 
        BIN=${BIN%.fa}
        # copy over precomputed crisprone output (but exclude certain files)
        if [ ! -d $CRISPR/$RTAG/$BIN ]; then 
            mkdir -p $CRISPR/$RTAG/$BIN
            cp $CRISPR/rNA/$BIN/*.summary-filter $CRISPR/$RTAG/$BIN
            cp $CRISPR/rNA/$BIN/*.domtblout $CRISPR/$RTAG/$BIN
            cp $CRISPR/rNA/$BIN/*.out $CRISPR/$RTAG/$BIN
        fi

        [ -e $CRISPR/$RTAG/$BIN/*.crt ] && rm -f $CRISPR/$RTAG/$BIN/*.crt 
        [ -e $CRISPR/$RTAG/$BIN/*.crt.ori ] && rm -f $CRISPR/$RTAG/$BIN/*.crt.ori
        [ -e $CRISPR/$RTAG/$BIN/*.des ] && rm -f $CRISPR/$RTAG/$BIN/*.des 
        [ -e $CRISPR/$RTAG/$BIN/*-sm.faa ] && rm -f $CRISPR/$RTAG/$BIN/*-sm.faa 
        [ -e $CRISPR/$RTAG/$BIN/*-sm.fna ] && rm -f $CRISPR/$RTAG/$BIN/*-sm.fna
        [ -e $CRISPR/$RTAG/$BIN/*-sm.gff ] && rm -f $CRISPR/$RTAG/$BIN/*-sm.gff
        [ -e $CRISPR/$RTAG/$BIN/*-sm.gff.sus ] && rm -f $CRISPR/$RTAG/$BIN/*-sm.gff.sus
        [ -e $CRISPR/$RTAG/$BIN/*-spacer.seq ] && rm -f $CRISPR/$RTAG/$BIN/*-spacer.seq

        FNA="$MAG_DIR/$BIN.fa"

        wait

        echo "BIN: $BIN"
        echo "REPEAT: $RTAG"

        # run CRISPRone with -repeat flag for guidance by repeats
        $python $CRISPRONE \
            --maxmis 2 \
            --outbase $CRISPR/$RTAG \
            --outprefix $BIN \
            --fna $FNA \
            --no_tracrRNA True \
            --repeat $CDHIT/$RTAG.fa | tee
    
        wait 

        # find all bins with crispr artifacts present
        if [ -s $CRISPR/$RTAG/$BIN/*.crt ]; then
            echo "$BIN" >> $NOTEMPTY
        fi
    done # bins inner loop
done # representative repeat outer loop

wait 

#####################################################################
# MAGs - repeat guides rALL (r0-r25 all together )
#####################################################################

# CRISPRONE for all representative repeats
RTAG="rALL"
#echo "Representative Repeat: rALL (all top 25 representative repeats [r0:r25])"

NOTEMPTY="$CRISPR/${RTAG}_assemblies_with_crispr.txt"
if [ -s $NOTEMPTY ]; then 
    rm -f $NOTEMPTY
fi
mkdir -p $CRISPR/$RTAG

B=$(ls $MAG_DIR)
for BIN in ${B[@]}; do
    BIN=${BIN%.fa}

    # copy over precomputed crisprone output (but exclude certain files) 
    if [ ! -d $CRISPR/$RTAG/$BIN ]; then 
        mkdir -p $CRISPR/$RTAG/$BIN
        cp $CRISPR/rNA/$BIN/*.summary-filter $CRISPR/$RTAG/$BIN
        cp $CRISPR/rNA/$BIN/*.domtblout $CRISPR/$RTAG/$BIN
        cp $CRISPR/rNA/$BIN/*.out $CRISPR/$RTAG/$BIN
    fi

    [ -e $CRISPR/$RTAG/$BIN/*.crt ] && rm -f $CRISPR/$RTAG/$BIN/*.crt 
    [ -e $CRISPR/$RTAG/$BIN/*.crt.ori ] && rm -f $CRISPR/$RTAG/$BIN/*.crt.ori
    [ -e $CRISPR/$RTAG/$BIN/*.des ] && rm -f $CRISPR/$RTAG/$BIN/*.des 
    [ -e $CRISPR/$RTAG/$BIN/*-sm.faa ] && rm -f $CRISPR/$RTAG/$BIN/*-sm.faa 
    [ -e $CRISPR/$RTAG/$BIN/*-sm.fna ] && rm -f $CRISPR/$RTAG/$BIN/*-sm.fna
    [ -e $CRISPR/$RTAG/$BIN/*-sm.gff ] && rm -f $CRISPR/$RTAG/$BIN/*-sm.gff
    [ -e $CRISPR/$RTAG/$BIN/*-sm.gff.sus ] && rm -f $CRISPR/$RTAG/$BIN/*-sm.gff.sus
    [ -e $CRISPR/$RTAG/$BIN/*-spacer.seq ] && rm -f $CRISPR/$RTAG/$BIN/*-spacer.seq

    FNA="$MAG_DIR/$BIN.fa"

    wait

    echo "BIN: $BIN"
    echo "REPEAT: $RTAG"

    # run CRISPRone with -repeat flag for guidance by repeats
    $python $CRISPRONE \
        --maxmis 2 \
        --outbase $CRISPR/$RTAG \
        --outprefix ${BIN%.gff} \
        --fna $FNA \
        --no_tracrRNA True \
        --repeat $REPRESENTATIVE_REPEATS
    wait
   
   # find all bins with crispr artifacts present
   if [ -s $CRISPR/$RTAG/$BIN/*.crt ]; then
        echo "$BIN" >> $NOTEMPTY
   fi
done 

#####################################################################
# INDIVIDUAL ASSEMBLIES  - repeat guides r0 to r25 (1 at a time )
#####################################################################

R=$(grep ">" $REPRESENTATIVE_REPEATS) # to loop over repeats
S=$(ls $WKDIR/INDIVIDUALS)

echo "Running CRISPRone (with repeat guides) for Individual Assemblies"
for REP in ${R[@]}; do # for each represenative repeat 
    RTAG=${REP##*>}
    # extracts repeat fasta from file
    # use python script to extract repeat contigs
    # repeat fastas already extracted from previously run script
    if [ ! -f $CHDIT/$RTAG.fa ]; then
        $python $call_script \
            -n $RTAG -f $REPRESENTATIVE_REPEATS -o $CDHIT
    fi
    #grep -wEA1 --no-group-separator $RTAG $REPRESENTATIVE_REPEATS > $CHDIT/$RTAG.fa
    NOTEMPTY="$CRISPR/${RTAG}_assemblies_with_crispr.txt"
    mkdir -p $CRISPR/$RTAG

    for SRR in ${S[@]}; do
	# copy over precomputed crisprone output (but exclude certain files) 
	if [ ! -d $CRISPR/$RTAG/$SRR ]; then
	    # copy over precomputed crisprone output (but exclude certain files) 
	    mkdir -p $CRISPR/$RTAG/$SRR 
	    cp $CRISPR/rNA/$SRR/*.summary-filter $CRISPR/$RTAG/$SRR 
	    cp $CRISPR/rNA/$SRR/*.domtblout $CRISPR/$RTAG/$SRR 
	    cp $CRISPR/rNA/$SRR/*.out $CRISPR/$RTAG/$SRR 
    fi 

    [ -e $CRISPR/$RTAG/$SRR/*.crt ] && rm -f $CRISPR/$RTAG/$SRR/*.crt 
    [ -e $CRISPR/$RTAG/$SRR/*.crt.ori ] && rm -f $CRISPR/$RTAG/$SRR/*.crt.ori
    [ -e $CRISPR/$RTAG/$SRR/*.des ] && rm -f $CRISPR/$RTAG/$SRR/*.des 
    [ -e $CRISPR/$RTAG/$SRR/*-sm.faa ] && rm -f $CRISPR/$RTAG/$SRR/*-sm.faa 
    [ -e $CRISPR/$RTAG/$SRR/*-sm.fna ] && rm -f $CRISPR/$RTAG/$SRR/*-sm.fna
    [ -e $CRISPR/$RTAG/$SRR/*-sm.gff ] && rm -f $CRISPR/$RTAG/$SRR/*-sm.gff
    [ -e $CRISPR/$RTAG/$SRR/*-sm.gff.sus ] && rm -f $CRISPR/$RTAG/$SRR/*-sm.gff.sus
    [ -e $CRISPR/$RTAG/$SRR/*-spacer.seq ] && rm -f $CRISPR/$RTAG/$SRR/*-spacer.seq

    FNA="$WKDIR/INDIVIDUALS/$SRR/ASSEMBLY/final_assembly.fasta"

    # run CRISPRone with -repeat flag for guidance by repeats
    $python $CRISPRONE \
        --maxmis 2 \
        --outbase $CRISPR/$RTAG \
        --outprefix $SRR  \
        --fna $FNA \
        --no_tracrRNA True \
        --repeat $CDHIT/$RTAG.fa 

    wait
    
	# find all indivs with crispr artifacts present
	if [ -s $CRISPR/$RTAG/$SRR/*.crt ]; then
	    echo "$SRR" >> $NOTEMPTY
	fi
    done # indivs inner loop
done # representative repeat outer loop


#####################################################################
# INDIVIDUAL ASSEMBLIES - repeat guides rALL (r0-r25 all together )
#####################################################################

# CRISPRONE for all representative repeats
RTAG="rALL"
#echo "Representative Repeat: rALL (all top 25 representative repeats [r0:r25])"

NOTEMPTY="$CRISPR/${RTAG}_assemblies_with_crispr.txt"
mkdir -p $CRISPR/$RTAG

S=$(ls $WKDIR/INDIVIDUALS)
#S=$(ls $INDIV)
for SRR in ${S[@]}; do
    
    # copy over precomputed crisprone output (but exclude certain files) 
   if [ ! -d $CRISPR/$RTAG/$SRR ]; then
        mkdir -p $CRISPR/$RTAG/$SRR
        cp $CRISPR/rNA/$SRR/*.summary-filter $CRISPR/$RTAG/$SRR
        cp $CRISPR/rNA/$SRR/*.domtblout $CRISPR/$RTAG/$SRR
        cp $CRISPR/rNA/$SRR/*.out $CRISPR/$RTAG/$SRR
    fi

    [ -e $CRISPR/$RTAG/$SRR/*.crt ] && rm -f $CRISPR/$RTAG/$SRR/*.crt 
    [ -e $CRISPR/$RTAG/$SRR/*.crt.ori ] && rm -f $CRISPR/$RTAG/$SRR/*.crt.ori
    [ -e $CRISPR/$RTAG/$SRR/*.des ] && rm -f $CRISPR/$RTAG/$SRR/*.des 
    [ -e $CRISPR/$RTAG/$SRR/*-sm.faa ] && rm -f $CRISPR/$RTAG/$SRR/*-sm.faa 
    [ -e $CRISPR/$RTAG/$SRR/*-sm.fna ] && rm -f $CRISPR/$RTAG/$SRR/*-sm.fna
    [ -e $CRISPR/$RTAG/$SRR/*-sm.gff ] && rm -f $CRISPR/$RTAG/$SRR/*-sm.gff
    [ -e $CRISPR/$RTAG/$SRR/*-sm.gff.sus ] && rm -f $CRISPR/$RTAG/$SRR/*-sm.gff.sus
    [ -e $CRISPR/$RTAG/$SRR/*-spacer.seq ] && rm -f $CRISPR/$RTAG/$SRR/*-spacer.seq

    FNA="$WKDIR/INDIVIDUALS/$SRR/ASSEMBLY/final_assembly.fasta"

    # run CRISPRone with -repeat flag for guidance by repeats
    $python $CRISPRONE \
        --maxmis 2 \
        --outbase $CRISPR/$RTAG \
        --outprefix $SRR \
        --fna $FNA \
        --no_tracrRNA True \
        --repeat $REPRESENTATIVE_REPEATS
    wait
   
   # find all bins with crispr artifacts present
   if [ -s $CRISPR/$RTAG/$SRR/*.crt ]; then
        echo "$SRR" >> $NOTEMPTY
   fi
done 

#####################################################################
# CO-ASSEMLBLY  - repeat guides r0 to r25 (1 at a time )
#####################################################################

R=$(grep ">" $REPRESENTATIVE_REPEATS) # to loop over repeats

echo "Running CRISPRone (with repeat guides) for Co-assembly"
for REP in ${R[@]}; do # for each represenative repeat 
    RTAG=${REP##*>}
    # extracts repeat fasta from file
    # use python script to extract repeat contigs
    # repeat fastas already extracted from previously run script
    if [ ! -f $CHDIT/$RTAG.fa ]; then
        $python $call_script \
            -n $RTAG -f $REPRESENTATIVE_REPEATS -o $CDHIT
    fi
    #grep -wEA1 --no-group-separator $RTAG $REPRESENTATIVE_REPEATS > $CHDIT/$RTAG.fa
    NOTEMPTY="$CRISPR/${RTAG}_assemblies_with_crispr.txt"
    mkdir -p $CRISPR/$RTAG
    SRR="coassembly"
    if [ ! -d $CRISPR/$RTAG/$SRR ]; then
        # copy over precomputed crisprone output (but exclude certain files) 
        mkdir -p $CRISPR/$RTAG/$SRR
        cp $CRISPR/rNA/$SRR/*.summary-filter $CRISPR/$RTAG/$SRR
        cp $CRISPR/rNA/$SRR/*.domtblout $CRISPR/$RTAG/$SRR
        cp $CRISPR/rNA/$SRR/*.out $CRISPR/$RTAG/$SRR
    fi

    [ -e $CRISPR/$RTAG/$SRR/*.crt ] && rm -f $CRISPR/$RTAG/$SRR/*.crt 
    [ -e $CRISPR/$RTAG/$SRR/*.crt.ori ] && rm -f $CRISPR/$RTAG/$SRR/*.crt.ori
    [ -e $CRISPR/$RTAG/$SRR/*.des ] && rm -f $CRISPR/$RTAG/$SRR/*.des 
    [ -e $CRISPR/$RTAG/$SRR/*-sm.faa ] && rm -f $CRISPR/$RTAG/$SRR/*-sm.faa 
    [ -e $CRISPR/$RTAG/$SRR/*-sm.fna ] && rm -f $CRISPR/$RTAG/$SRR/*-sm.fna
    [ -e $CRISPR/$RTAG/$SRR/*-sm.gff ] && rm -f $CRISPR/$RTAG/$SRR/*-sm.gff
    [ -e $CRISPR/$RTAG/$SRR/*-sm.gff.sus ] && rm -f $CRISPR/$RTAG/$SRR/*-sm.gff.sus
    [ -e $CRISPR/$RTAG/$SRR/*-spacer.seq ] && rm -f $CRISPR/$RTAG/$SRR/*-spacer.seq

    FNA="$COASSEMBLY"

    # run CRISPRone with -repeat flag for guidance by repeats
    $python $CRISPRONE \
        --maxmis 2 \
        --outbase $CRISPR/$RTAG \
        --outprefix $SRR \
        --fna $FNA \
        --no_tracrRNA True \
        --repeat $CDHIT/$RTAG.fa 

    wait

    # find all indivs with crispr artifacts present
    if [ -s $CRISPR/$RTAG/$SRR/*.crt ]; then
        echo "$SRR" >> $NOTEMPTY
    fi
done # representative repeat outer loop

#####################################################################
# CO-ASSEMBLY - repeat guides rALL (r0-r25 all together )
#####################################################################

# CRISPRONE for all representative repeats
RTAG="rALL"
#echo "Representative Repeat: rALL (all top 25 representative repeats [r0:r25])"

NOTEMPTY="$CRISPR/${RTAG}_assemblies_with_crispr.txt"
mkdir -p $CRISPR/$RTAG

SRR="coassembly"
    
# copy over precomputed crisprone output (but exclude certain files) 
if [ ! -d $CRISPR/$RTAG/$SRR ]; then
   mkdir -p $CRISPR/$RTAG/$SRR
   cp $CRISPR/rNA/$SRR/*.summary-filter $CRISPR/$RTAG/$SRR
   cp $CRISPR/rNA/$SRR/*.domtblout $CRISPR/$RTAG/$SRR
   cp $CRISPR/rNA/$SRR/*.out $CRISPR/$RTAG/$SRR
fi

[ -e $CRISPR/$RTAG/$SRR/*.crt ] && rm -f $CRISPR/$RTAG/$SRR/*.crt 
[ -e $CRISPR/$RTAG/$SRR/*.crt.ori ] && rm -f $CRISPR/$RTAG/$SRR/*.crt.ori
[ -e $CRISPR/$RTAG/$SRR/*.des ] && rm -f $CRISPR/$RTAG/$SRR/*.des 
[ -e $CRISPR/$RTAG/$SRR/*-sm.faa ] && rm -f $CRISPR/$RTAG/$SRR/*-sm.faa 
[ -e $CRISPR/$RTAG/$SRR/*-sm.fna ] && rm -f $CRISPR/$RTAG/$SRR/*-sm.fna
[ -e $CRISPR/$RTAG/$SRR/*-sm.gff ] && rm -f $CRISPR/$RTAG/$SRR/*-sm.gff
[ -e $CRISPR/$RTAG/$SRR/*-sm.gff.sus ] && rm -f $CRISPR/$RTAG/$SRR/*-sm.gff.sus
[ -e $CRISPR/$RTAG/$SRR/*-spacer.seq ] && rm -f $CRISPR/$RTAG/$SRR/*-spacer.seq

FNA="$COASSEMBLY"

# run CRISPRone with -repeat flag for guidance by repeats
$python $CRISPRONE \
    --maxmis 2 \
    --outbase $CRISPR/$RTAG \
    --outprefix $SRR \
    --fna $FNA \
    --no_tracrRNA True \
    --repeat $REPRESENTATIVE_REPEATS
wait
 
# find all bins with crispr artifacts present
if [ -s $CRISPR/$RTAG/$SRR/*.crt ]; then
    echo "$SRR" >> $NOTEMPTY
fi


# References:
# https://unix.stackexchange.com/questions/253499/extracting-subset-from-fasta-file

