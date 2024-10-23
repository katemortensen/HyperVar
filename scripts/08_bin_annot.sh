#!/bin/bash

# Author: Kate Mortensen
# Last Modified: 10/21/2024
# Purpose: execute MetaWRAP's functional annotation module 

# setup variables
# . 00_vars.sh

# conda env
# module load anaconda
# source activate metawrap-env

BIN_DIR=$1
BIN_ANNOTATION=$2
THREADS=$3

mkdir -p $BIN_ANNOTATION    

metawrap annotate_bins \
    -b $BIN_DIR
    -o $BIN_ANNOTATION \
    -t $THREADS \



