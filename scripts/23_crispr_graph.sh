#!/bin/bash

# Author: Kate Mortensen
# Last Modified: 07/15/2023
# Purpose: Create CRISPR graphs

# setup variables
#. 00_vars.sh
#source ./00_vars.sh

python=$1
CRISPRGRAPH=$2
WKDIR=$3
SUBWKDIR=$4
CRISPR_GRAPHS=$5
BIN_REASSEMBLY=$6
INDIV=$7
CDHIT=$8
HYPERVAR_PIPELINE=$9
COASSEMBLY_DIR=$10


echo "create CRISPR_GRAPH directory and meta file and color file"
mkdir -p $CRISPR_GRAPHS
if [ -s $CRISPR_GRAPHS/meta.txt ];then
    rm -f $CRISPR_GRAPHS/meta.txt
fi
touch $CRISPR_GRAPHS/meta.txt
cp $WKDIR/bin/color.txt $CRISPR_GRAPHS/.

echo "add indiv assembly paths to meta file"
COUNTER=0
mapfile -t C < $HYPERVAR_PIPELINE/color.txt
M=$(ls $BIN_REASSEMBLY/reassembled_bins/*.fa)
I=$(cat $WKDIR/SRR_list.txt)
for acc in ${I[@]}; do
    color=${C[$COUNTER]: -3}
    acc_path=$INDIV/$acc/ASSEMBLY/final_assembly.fasta
    echo -e "${acc}\t${color}\t${acc_path}" >> $CRISPR_GRAPHS/meta.txt
    COUNTER=$COUNTER+1
done

wait

echo "add mag paths to meta file"
M=$(ls $SUBWKDIR/bins.txt)
for acc in ${M[@]}; do
    acc=${acc%.fa}
    color=${C[$COUNTER]: -3}
    acc_path=$BIN_REASSEMBLY/reassembled_bins/${acc}.fa
    echo -e "${acc}\t${color}\t${acc_path}" >> $CRISPR_GRAPHS/meta.txt
done

wait 

echo "add coassembly path to meta file"
echo -e "${COASSEMBLY_DIR}/final_assembly.fasta" >> $CRISPR_GRAPHS/meta.txt

wait 

echo "copy repeat files from CDHIT to CRISPR_GRAPHS directory"
cp $CDHIT/r*.fa $CRISPR_GRAPHS/.

cd $CRISPR_GRAPHS

echo "run crispr_ann_ort.py"
R=$(ls $CRISPR_GRAPHS/*.fa)
for REP in ${R[@]}; do 
    RTAG=${REP##$CRISPR_GRAPHS/}
    RTAG=${RTAG%.fa}
    echo "$RTAG"
    $python $CRISPRGRAPH \
	--input $CRISPR_GRAPHS/meta.txt \
	--color $CRISPR_GRAPHS/color.txt \
	--repeat $REP \
	--prog CRISPRAlign \
	--out $RTAG
done

echo "compensate for the labeling issue/bug from crispr_ann_ort.py"
R=$(ls $CRISPR_GRAPHS/*.fa)
for REP in ${R[@]}; do
    RTAG=${REP##$CRISPR_GRAPHS/}
    RTAG=${RTAG%.fa}
    mv $CRISPR_GRAPHS/$RTAG/spacer_graph/repeat_nr.fa \
	    $CRISPR_GRAPHS/$RTAG/spacer_graph/repeat_nr
    mv $CRISPR_GRAPHS/$RTAG/spacer_graph/repeat_nr.fa.clstr  \
	    $CRISPR_GRAPHS/$RTAG/spacer_graph/repeat_nr.clstr
    mv $CRISPR_GRAPHS/$RTAG/spacer_graph/repeat_nr100.fa \
	    $CRISPR_GRAPHS/$RTAG/spacer_graph/repeat_nr100
    mv $CRISPR_GRAPHS/$RTAG/spacer_graph/spacer_nr.fa \
	    $CRISPR_GRAPHS/$RTAG/spacer_graph/spacer_nr
    mv $CRISPR_GRAPHS/$RTAG/spacer_graph/spacer_nr.fa.clstr  \
	    $CRISPR_GRAPHS/$RTAG/spacer_graph/spacer_nr.clstr
    mv $CRISPR_GRAPHS/$RTAG/spacer_graph/spacer_nr100.fa \
	    $CRISPR_GRAPHS/$RTAG/spacer_graph/spacer_nr100
    mv $CRISPR_GRAPHS/$RTAG/spacer_graph/spacer_nr100.fa \
	    $CRISPR_GRAPHS/$RTAG/spacer_graph/spacer_nr100
    mv $CRISPR_GRAPHS/$RTAG/spacer_graph/spacer_nr100.fa.clstr \
	    $CRISPR_GRAPHS/$RTAG/spacer_graph/spacer_nr100.clstr
done


echo "create pdf from dot file"
R=$(ls $CRISPR_GRAPHS/*.fa)
for REP in ${R[@]}; do
    RTAG=${REP##$CRISPR_GRAPHS/}
    RTAG=${RTAG%.fa}

    dot -Tpdf -o $CRISPR_GRAPHS/$RTAG/spacer_graph/${RTAG}.pdf \
	$CRISPR_GRAPHS/$RTAG/spacer_graph/spacer_graph_compressed.dot
done


cd $HYPERVAR_PIPELINE/scripts
