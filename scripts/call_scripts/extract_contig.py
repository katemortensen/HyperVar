#!/usr/bin/env python3

# Author: Kate Mortensen
# Date: 01-31-2022
Purpose = "Purpose: to extract specified contig from specified assembly and save to file"

import Bio.SeqIO as IO
import pandas as pd
import argparse
from sys import argv
import os.path

# arguments
parser = argparse.ArgumentParser(description=Purpose)
parser.add_argument('-f', type = str, help = "enter path to assembly fasta file from which you would like to extract the contig", required = True)
parser.add_argument('-n', type = str, help = "enter contig header of contig to extract", required = True)
parser.add_argument('-o', type = str, help = "enter output path", required = True)
args = parser.parse_args()

# input fasta file 
record_dict = IO.to_dict(IO.parse(args.f, "fasta"))
names = record_dict.keys()

# longest contig in assembly
#longest = list(names)[0]

# write down as a fasta file
sample_fasta = [record_dict[args.n]]
#output_file_path = os.path.dirname(args.n) + "/" + "%s.fa" % (args.n)
output_file_path = args.o + "/" + "%s.fa" % (args.n)
IO.write(sample_fasta, output_file_path , "fasta")



# Resources:
# https://ravinpoudel.github.io/GenomeQuest/README.html
# https://stackoverflow.com/questions/63106413/find-length-of-a-contig-in-one-fasta-using-the-header-of-another-fasta-as-query
