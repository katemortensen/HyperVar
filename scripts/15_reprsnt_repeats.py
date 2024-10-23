#!/usr/bin/env python
# coding: utf-8

# Author: Kate Mortensen
# Last Modified: 8/29/2023
# Purpose: To collect representative repeats for the top 25 most populous clusters of repeats 

# %% 
import pandas as pd
import numpy as np
import re
import pathlib
import argparse

# %% 
parser = argparse.ArgumentParser()
parser.add_argument("--clstr", metavar="clstr-file", help="/path/to/CD-HIT-EST/cd-hit-est.out.clstr", required=True)
parser.add_argument("--repcon", metavar="repeat-consensus", help="/path/to/CD-HIT-EST/concat_repeatcon.seq", required=True)
parser.add_argument("--outdir", metavar="output-dir", help="CD-HIT-EST directory", required=True)
parser.add_argument("--topmost_rep_guides", metavar="top-reps", help="CD-HIT-EST directory", required=True)

args = parser.parse_args()

# directory with CD-HIT-EST outputs, etc.
outdir = args.outdir+"/"
clstr = args.clstr
repcon = args.repcon
topmost_reps = int(args.topmost_rep_guides)

# %% 
# test vars
# outdir = "/home/kmorten/ADMB/DASTOOL_MAGS/CD-HIT-EST/" 
# clstr = "/home/kmorten/ADMB/DASTOOL_MAGS/CD-HIT-EST/cd-hit-est.out.clstr"
# repcon = "/home/kmorten/ADMB/DASTOOL_MAGS/CD-HIT-EST/concat_repeatcon.seq"
# topmost_reps = 25

# %% 
'''node2seq'''
'''nodeID > seq'''

node2seq = {}
with open(repcon) as f:
    for line in f:
        if line[0] == '>':
            node = line.strip('\n').strip('>')
        if line[0] != '>':
            seq = line.strip('\n')
            if node2seq.get(node) != None:
                assert node2seq.get(node) == seq
            else:
                node2seq[node] = seq
    

# %% 
'''rep2seq'''
'''rep > nodeID > seq'''

rep2seq = {}
rep_repeats = []

counter = 0
with open(clstr) as f:
    switch = 0
    for line in f:
        if line[0] == '>' and counter < topmost_reps:
            clust = line.strip('\n').strip(">Cluster ")
            rep = '>r%s' % (clust)
            counter += 1 
            switch = 0
        if line[0] != '>' and switch == 0:
            star = line.split(">")[1].strip('\n')[-1]
            if star == '*':
                nodeID = line.split(">")[1].strip('\n').split("...")[0]
                rep2seq[rep] = {nodeID:seq}
                rep_repeats.append(nodeID)
                switch += 1
assert len(rep2seq.keys()) == len(rep_repeats) == topmost_reps

# %% 
'''node2seqlist'''
'''node > list of sequences'''
'''purpose is to make sure there aren't multiple sequences for the same nodeID'''

node2seqlist = {}

with open(clstr) as f:
    for line in f:
        if line[0] != '>':
            nodeID = line.split(">")[1].strip('\n').split("...")[0]
            if node2seqlist.get(nodeID) == None:
                node2seqlist.update({nodeID:[]})
                 
left_out = []
with open(repcon) as f:
    for line in f:
        if line[0] == '>':
            node = line.strip('\n').strip('>')
        if line[0] != '>':
            seq = line.strip('\n')           
        if node2seqlist.get(node) == None:
            left_out.append(node)
        elif node2seqlist.get(node) != None:
            temp = node2seqlist.get(node)
            node2seqlist.update({node:temp.append(seq)})

for node in node2seqlist.keys():
    if node2seqlist.get(node) != None:
        temp = len(node2seqlist.get(node))
        assert temp == 1            


# %% 
# Create a dictionary of each cluster and their node location
arr = np.genfromtxt(clstr, delimiter="\n", dtype="str")

cluster=[]
repeats=[]
counter=0
clust_repeats = [] 
for i in arr :
    if ">" in i[0] :
        cluster.append(i.strip(">Cluster "))
        if counter > 0 :
            repeats.append(clust_repeats)
        clust_repeats = []
        counter += 1   
    else :
        clust_repeats.append(i.split(" ")[1].strip("..."))
repeats.append(clust_repeats)
assert len(repeats) == len(cluster)


# %% 
# topmost representative repeats
clust2rep = dict(zip(cluster[0:topmost_reps], rep_repeats[0:topmost_reps]))

# Find the sequence for each node location and write to a file
arr = np.genfromtxt(repcon, delimiter="\n", dtype="str")

reps = []
seqs = []
for key_cluster in list(clust2rep.keys()) :
    value_node = clust2rep.get(key_cluster)
    for i in range(len(arr)):
        if value_node in arr[i]:
            if value_node not in reps:
                reps.append(value_node)
                seqs.append(arr[i+1])
                counter += 1
            else:
                assert seqs[-1] == arr[i+1]
        
            
assert len(reps)==len(seqs)==len(clust2rep)

# %% 

# write nodes and seqs to a file
path_out = outdir + "top%s_representative_repeats.fa" % (topmost_reps)
file = open(path_out, 'w') 
#write to file
for i in range(len(list(clust2rep.keys()))) :
   #  header = ">Cluster %s Representative Repeat" %i
     header = ">r%s" %i
     file.write(header+"\n")
     file
     file.write(seqs[i]+"\n")
file.close()




# %%
