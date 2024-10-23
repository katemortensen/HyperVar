#!/usr/bin/env python
# coding: utf-8

# In[2]:


# Author: Kate Mortensen
# Date: 1/17/2023
Purpose = "Purpose: to find the average weigthed abundance profile of contigs in MAGs with respect of sample."


# In[1]:


import pandas as pd 
import numpy as np
from scipy.spatial import distance
from sklearn.metrics import pairwise_distances
from sklearn.metrics.pairwise import pairwise_kernels
import argparse


# In[7]:

# Initialize parser
parser = argparse.ArgumentParser()
 
# Adding optional argument
parser.add_argument("-abund_dir", "--abund_dir", help = "path to directory where CONTIG_RELATIVE_ABUNDANCE output exists (ex: \'/usr/home/wkdir/subwkdir/abund_dir/CONTIG_RELATIVE_ABUNDANCE\')", required=True)
parser.add_argument("-rrcpm_profile", "--rrcpm_profile", help = "path to relative contig abundance profile (salmon or minimap)", required = True )
parser.add_argument("-mapping_tool", "--mapping_tool", help = "(salmon or minimap)", required = True )

args = parser.parse_args()
abund_dir = args.abund_dir
mapping_tool = args.mapping_tool
rrcpm_profile = args.rrcpm_profile 
#print('CONTIG_RELATIVE_ABUNDANCE directory:', abund_dir)
# stop and comment out hardcoded params


# %% 
# mapping_tool = 'minimap'
# mapping_tool = 'salmon'

# abund_dir = '/home/kmorten/ADMB/METAWRAP_MAGS/CONTIG_RELATIVE_ABUNDANCE'
# rrcpm_profile = '/home/kmorten/ADMB/METAWRAP_MAGS/CONTIG_RELATIVE_ABUNDANCE/relative_contig_abundance_profiles_crt_%s.tsv' % (mapping_tool)

# abund_dir = '/home/kmorten/ADMB/DASTOOL_MAGS/CONTIG_RELATIVE_ABUNDANCE'
# rrcpm_profile = '/home/kmorten/ADMB/DASTOOL_MAGS/CONTIG_RELATIVE_ABUNDANCE/relative_contig_abundance_profiles_crt_%s.tsv' % (mapping_tool)


# %%

# ####  1) sort contigs by length within each assembly 
# ####  2) N50 - calculate theoretical N50 length
# ####  3) N50 contigs - the smallest array of contigs whose cumulative length accounts for at least 50% of the total cumulative length of all contigs
# ####  4) A50 profile

# In[8]:

# Dictionary creation for Average Abundance Profile

'''samples'''
'''src (sample read counts) = total read counts for each sample '''
'''NOTE: Outputs verified for section below. (KM 12/20/2022)'''
src_path = abund_dir + '/sample_reads_count.txt'
src_df = pd.read_csv( src_path , sep="\t")
src_dict = dict(zip(src_df['sample'], src_df['total_reads']))
samples = list(src_df['sample'])

'''contig2len'''
'''NOTE: Outputs verified for section below. (KM 12/20/2022)'''
contig2len = {}
header_accounted = 0
for line in open(abund_dir + "/nodeIDs.txt" , 'r'):
    if header_accounted != 0:
        linep = line.strip("\n").split("\t")
        contig = linep[1] # primary key
        len = int(contig.split("_length_")[1].split("_cov_")[0])
        contig2len[contig] = {'len':len}
    else: 
        header_accounted += 1

'''assembly2len & assembly2contigs'''
'''assembly2len = assembly->total_len,N50_len,N50_actual_len'''
'''assembly2contigs = assembly->N50_len,contigs->len,weight,N50_contig(1/0)'''
'''NOTE: Outputs verified for section below. (KM 12/20/2022)'''
path_in = rrcpm_profile 
assembly2len = {}
assembly2contigs = {}
header_accounted = 0
cols = []
# add total assembly len and contig len 
file = open(path_in, 'r')
for line in file:
    linep = line.strip("\n").split("\t")
    if header_accounted == 0: 
        header_accounted += 1
        cols = linep
    else:
        dict_tmp = dict(zip(cols,linep))
        assembly = dict_tmp.get('assembly')
        contig = dict_tmp.get('contig')
        contig_len = contig2len.get(contig)['len']

        # update assembly2len
        total_len = 0
        if assembly2len.get(assembly) == None:
            assembly2len[assembly] = {'total_len':0}
        total_len = assembly2len.get(assembly)['total_len']
        new_total_len = total_len + contig_len
        assembly2len[assembly] = {'total_len':new_total_len}

        # add contig_len to assembly2contig
        if assembly2contigs.get(assembly) == None:
            assembly2contigs[assembly] = {'contigs':{contig:{'len':contig_len, 'N50_contig':0}}}
        else:
            assembly2contigs[assembly]['contigs'].update({contig:{'len':contig_len, 'N50_contig':0}})
# add N50_len to assembly2len
for assembly in list(assembly2len.keys()): 
    N50_len = int(assembly2len[assembly]['total_len']/2)
    assembly2len[assembly].update({'N50_len':N50_len})
# N50 contigs & weight term
for assembly in list(assembly2contigs.keys()): 
    # sort contigs in assembly by descending length
    sorted_contig_ls = sorted(assembly2contigs[assembly]['contigs'], key=lambda x:x[1], reverse=True)
    N50_len = assembly2len[assembly]['N50_len']
    total_len = assembly2len[assembly].get('total_len')
    len_accounted = 0
    for contig in sorted_contig_ls:
        if len_accounted <= N50_len:
            # label N50 contig
            assembly2contigs[assembly]['contigs'][contig]['N50_contig'] = 1
            contig_len = assembly2contigs[assembly]['contigs'][contig]['len']
            len_accounted += contig_len
    assert len_accounted < total_len
    # print("N50_len: ", N50_len, "len_accounted: ", len_accounted)
    # add N50_len_actual to assembly2len
    assembly2len[assembly].update({'N50_actual_len':len_accounted})
    # add weight to assembly2contigs
    for contig in sorted_contig_ls:
        contig_len = assembly2contigs[assembly]['contigs'][contig]['len']
        weight = contig_len/N50_len
        assembly2contigs[assembly]['contigs'][contig].update({'weight':weight})

'''assembly2profile'''
'''Weighted & Non-weighted Relative Redistributed Read Counts Per Million '''
'''assembly->contig->sample->wrrcpm ,rrcpm'''
'''wrrcpm  = weighted (and redistributed) relative counts per million'''
path_in = rrcpm_profile
path_out = abund_dir + "/weighted_relative_contig_abundance_profiles_crt_%s.tsv" % (mapping_tool)

assembly2profile = {}
header_accounted = 0
cols = []

with open(path_out, "w") as fp_out:
    # apply weight term to rrcpm 
    file = open(path_in, 'r')
    for line in file:
        linep = line.strip("\n").split("\t")
        if header_accounted == 0:
            cols = linep
            fp_out.write(line)
            header_accounted += 1
        else:
            lineDict_tmp = dict(zip(cols,linep))
            assembly_type = lineDict_tmp.get('assembly_type')
            assembly = lineDict_tmp.get('assembly')
            contig = lineDict_tmp.get('contig')
            crt = lineDict_tmp.get('crt')
            weight = assembly2contigs[assembly]['contigs'][contig].get('weight')
            sample2profile_tmp = {}
            for sample in samples:
                rrcpm = float(lineDict_tmp.get(sample))
                wrrcpm  = weight * rrcpm
                sample2profile_tmp.update({sample:{'wrrcpm ':wrrcpm , 'rrcpm':rrcpm}})
            # update assembly2profile with entire sample2profile_tmp dictionary
            if assembly2profile.get(assembly) == None:
                assembly2profile[assembly] = {'contig':{contig:{'sample':sample2profile_tmp}}}
            else:
                assembly2profile[assembly]['contig'][contig] = {'sample':sample2profile_tmp}
            # write contig wrrcpm  profile to file 
            new_row = [assembly_type, assembly, contig, crt] + [str(assembly2profile[assembly]['contig'][contig]['sample'][sample].get('wrrcpm ')) for sample in samples]
            new_row_string = "\t".join(new_row).strip("\t") + "\n"
            fp_out.write(new_row_string)               

'''assembly2A50'''
'''Average Weighted Abundance Profile of N50 contigs with respect to all metagenomic samples that were used in co-assembly'''
'''awrrcpm = average weighted (redistributed) relative counts per million across contigs for each sample'''
'''assembly->sample->awrrcpm'''
path_out = abund_dir + '/A50_concensus_profile_%s.tsv' % (mapping_tool)
assembly2A50 = {}
# write A50 consensus profiles to file 
with open(path_out, "w") as fp_out:
    # write header
    cols = ['assembly'] + samples
    header = "\t".join(cols).strip("\t") + "\n"
    fp_out.write(header)
    for assembly in list(assembly2profile.keys()):
        sample2awrrcpm_tmp = {}
        toggle = 0
        contigs = list(assembly2profile[assembly]['contig'].keys())
        for contig in contigs:
            N50_contig = assembly2contigs[assembly]['contigs'][contig].get('N50_contig')
            if N50_contig == 1:
                if toggle == 0:
                    for sample in samples:
                        wrrcpm  = assembly2profile[assembly]['contig'][contig]['sample'][sample].get('wrrcpm ')
                        sample2awrrcpm_tmp[sample] = {'awrrcpm':wrrcpm }
                toggle += 0
                for sample in samples:
                    old_awrrcpm = sample2awrrcpm_tmp[sample].get('awrrcpm')
                    next_wrrcpm  = assembly2profile[assembly]['contig'][contig]['sample'][sample].get('wrrcpm ')
                    awrrcpm = (old_awrrcpm + next_wrrcpm )/2
                    sample2awrrcpm_tmp[sample].update({'awrrcpm':awrrcpm})
        assembly2A50[assembly] = {'sample':sample2awrrcpm_tmp}
        new_row = [assembly] + [str(i) for i in [assembly2A50[assembly]['sample'][sample].get('awrrcpm') for sample in samples]]
        new_row_string = "\t".join(new_row).strip("\t") + "\n"
        fp_out.write(new_row_string)



# In[9]:
# ## Abundance 
# #### 1) Average Abundance = contig rrcpm values averaged across sample 
# #### 2) Average Weighted Abundance = contig wrrcpm values across sample 
# Note: The output of this section is only used for quick merging and metadata later on and distribution plots of awrrcpm and arrcpm 

'''assembly2abund'''
'''Average rrcpm and wrrcpm  across samples for each contig.'''
'''Each contig has 1 value to describe abundance (for weighted and non-weigthed rrcpm values respectively'''
'''average rrcpm across samples and average wrrcpm  across samples'''
'''assembly->contig->arrcpm, awrrcpm'''
assembly2abund = {}
path_out = abund_dir + '/contig_abund_%s.tsv' % (mapping_tool)
with open(path_out, "w") as fp_out:
    # write header
    cols = ['assembly','contig','N50_contig_boolean','arrcpm','awrrcpm']
    header = "\t".join(cols).strip("\t") + "\n"
    fp_out.write(header)
    for assembly in list(assembly2contigs.keys()):
        contig2abund_tmp = {}
        for contig in list(assembly2contigs[assembly]['contigs'].keys()):
            rrcpm_ls = []
            wrrcpm_ls = []
            sample_count = 0
            for sample in samples:
                rrcpm_ls.append(assembly2profile[assembly]['contig'][contig]['sample'][sample].get('rrcpm'))
                wrrcpm_ls.append(assembly2profile[assembly]['contig'][contig]['sample'][sample].get('wrrcpm '))
                sample_count += 1
            arrcpm = sum(rrcpm_ls)/sample_count
            awrrcpm = sum(wrrcpm_ls)/sample_count
            N50_contig_boolean = str(assembly2contigs[assembly]['contigs'][contig].get('N50_contig'))
            contig2abund_tmp.update({contig:{'arrcpm':arrcpm, 'awrrcpm':awrrcpm}})

            # write line to file 
            new_row = [assembly, contig, N50_contig_boolean, str(arrcpm), str(awrrcpm)]
            new_row_string = "\t".join(new_row).strip("\t") + "\n"
            fp_out.write(new_row_string)
            
        # update assembly2abund with entire contig2abund_tmp dictionary
        assembly2abund[assembly] = {'contig':contig2abund_tmp}


# In[10]:

# ## Cosine Distance
# #### 1) Cosine distance between all contigs in an assembly and their A50 profile 
# #### 2) Label contigs with outlier cosine distances from A50 profile based on N50 outlier standards

'''assembly2cos'''
'''A50_cosine_distance'''
'''assembly2cos = assembly->contig->A50_cosine_dist,outlier(0/1)'''
'''outliers are based on cosine distances of all contigs awrrcpm to A50 profiles (not all contigs)'''
'''To re-iterate, '''

assembly2cos = {}
path_out = abund_dir + '/A50_cosine_distance_%s.tsv' % (mapping_tool)
with open(path_out, "w") as fp_out:
    # write header
    cols = ['assembly','contig','A50_cosine_dist','outlier', 'N50_contig_boolean']
    header = "\t".join(cols).strip("\t") + "\n"
    fp_out.write(header)
    # cosine distance 
    for assembly in list(assembly2contigs.keys()):
        # A50 consensus profile for assembly (average weighted abundance profile of assembly)
        A50_tmp = [assembly2A50[assembly]['sample'][sample].get('awrrcpm') for sample in samples]
        # COSINE DISTANCE LOOP BELOW
        for contig in list(assembly2profile[assembly]['contig'].keys()):
            # relative abundance profile of contig 
            rrcpm_tmp = [assembly2profile[assembly]['contig'][contig]['sample'][sample].get('rrcpm') for sample in samples]
            # cosine distance 
            A50_cosine_dist = distance.cosine(A50_tmp, rrcpm_tmp)
            # N50_contig (0/1)
            N50_contig = assembly2contigs[assembly]['contigs'][contig].get('N50_contig')
            # add cos_dist to assembly2cos
            if assembly2cos.get(assembly) == None:
                assembly2cos[assembly] = {'contig':{contig:{'A50_cosine_dist':A50_cosine_dist}}}
            else:
                assembly2cos[assembly]['contig'].update({contig:{'A50_cosine_dist':A50_cosine_dist}})
    # outlier 
    for assembly in list(assembly2contigs.keys()):
        N50_contig_distances = []
        for contig in list(assembly2contigs[assembly]['contigs'].keys()):
            # outlier threshold based on N50 contigs only
            if assembly2contigs[assembly]['contigs'][contig].get('N50_contig') == 1:
                A50_cosine_dist = assembly2cos[assembly]['contig'][contig].get('A50_cosine_dist')
                N50_contig_distances.append(A50_cosine_dist)
        # calculate interquartile range (based on N50 contig distances)
        # q1 = np.percentile(N50_contig_distances, 25, interpolation='midpoint')
        # q3 = np.percentile(N50_contig_distances, 75, interpolation='midpoint')
        q1 = np.percentile(N50_contig_distances, 25)
        q3 = np.percentile(N50_contig_distances, 75)
        iqr = q3 - q1
        # Upper bound
        upper = q3+2*iqr
        # Lower bound
        lower = q1-2*iqr

        tot_assembly_len = assembly2len[assembly].get('total_len')
        contigs_sorted_by_dist = sorted(assembly2cos[assembly]['contig'].items(), key=lambda x:x[1]['A50_cosine_dist'], reverse=True)
        culled_len = 0
        running_len = 0
        total_so_far = 0
        for i in contigs_sorted_by_dist:
            contig = i[0]
            A50_cosine_dist = assembly2cos[assembly]['contig'][contig].get('A50_cosine_dist')
            current_contig_len = assembly2contigs[assembly]['contigs'][contig].get('len')
            total_so_far += current_contig_len
            if A50_cosine_dist > upper and (culled_len + current_contig_len) < tot_assembly_len*0.05:
            #if A50_cosine_dist < lower or A50_cosine_dist > upper:
                culled_len += current_contig_len
                assembly2cos[assembly]['contig'][contig].update({'outlier':1})
            else:
                running_len += current_contig_len
                assembly2cos[assembly]['contig'][contig].update({'outlier':0})
            assert culled_len + running_len == total_so_far
            outlier = assembly2cos[assembly]['contig'][contig].get('outlier')
            N50_contig_boolean = assembly2contigs[assembly]['contigs'][contig].get('N50_contig')
            new_row = [assembly, contig, str(A50_cosine_dist), str(outlier), str(N50_contig_boolean)]
            new_row_string = "\t".join(new_row).strip("\t") + "\n"
            fp_out.write(new_row_string)
        assert culled_len + running_len == tot_assembly_len
        assert culled_len <= tot_assembly_len*0.05
        assert running_len >= tot_assembly_len*0.95
    


# # FINAL OUTPUT
# #### weighted_relative_contig_abundance_profiles_crt.tsv
# #### A50_concensus_profile.tsv 
# #### contig_abund.tsv
# #### A50_cosine_distance.tsv: assembly/contig/A50_cosine_dist/IQR_outlier(0/1)

# %%
