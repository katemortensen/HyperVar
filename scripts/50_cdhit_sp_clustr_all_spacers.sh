#!/bin/bash

# Author: Kate Mortensen
# Last Modified: 1/11/2024
# Purpose: cluster all spacers (no repeat guide) using CD-HIT-EST


cd_hit_est=$1
SPACERS_FILE=$2
FILE_OUT=$3

echo "runnning cd-hit-est on concated spacers for all assemblies and methodolies (no repeat guide)"

touch $FILE_OUT
$cd_hit_est -i $SPACERS_FILE \
-sc 1 -sf -1 -d 200 -o $FILE_OUT 
    



