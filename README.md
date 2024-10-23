# Hypervariable-region-aware-co-assembly-of-metagenomes


## Overview


Hypervariable regions in the genome are often lost during the assembly process for short read data. These regions can hold important information, such as, DNA from invading entities maintained by the CRISPR-Cas immune system. Co-assemlby of metagenome-assembled genomes (MAGs) takes advantage of sequencing data from multiple samples to assemble more complete taxa level metagenomes and better represent underrepresented species. Procurement of MAGs involves a method of sorting assembled contigs into various "bins" where each bin is representative of a taxon. Short read, hypervariable regions may not sort in a similar manner as other DNA native to the taxons. This pipeline allows for observation of the CRISPR artifacts, especially the hypervariable sapcer regions, in MAG procurement. 

The HyperVar pipeline takes Illumina short read data from many accessions (or environmental samples) and procures a cohort of MAGs using MetaWRAP. HyperVar then creates a relative abudance profile for each contig which serves as an accession-based fingerprint where each vector entry in the profile is with respect to an accession number. This profile is used to analyze the cohort of MAGs and contigs, especially in context of spacers, along with various other analyses outlined in the pipeline tutorial. 


## Software Requirements 


- Linux 
- Python 3.8
- R 3.6.3


## Installation and Setup


Note that HyperVar's bin.tar.gz file (below) is 2.5G zipped and 7G unzipped.

```
git clone https://github.com/katemortensen/HyperVar.git
wget -P ./HyperVar/. https://omics.informatics.indiana.edu/HyperVar/hypervar-pipeline/bin.tar.gz
tar -xzvf ./HyperVar/hypervar-pipeline/bin.tar.gz -C ./HyperVar/hypervar-pipeline/.
rm ./HyperVar/hypervar-pipeline/bin.tar.gz
```

## Conda Environment


```
conda env create -f ./HyperVar/hypervar-pipeline/scripts/hypervar_conda_env.yml 
conda activate hypervar_env
```


## Usage

First, prepare a directory for raw reads. Name this directory "INDIVIDUALS".

```
mkdir INDIVIDUALS
cd INDIVIDUALS
```

Each sub-directory in "INDIVIDUALS" must have the name of an accession number (ex: ./INDIVIDUALS/SRR1/, ./INDIVIDUALS/SRR2/, ..., ./INDIVIDUALS/SRR10. Download the forward and reverse fastq files for each accession to their appropriate directories. Note that the fastq files must be g-zipped such that their extension is ".gz". The name of the forward and reverse reads must be in the following format: <SRR_number_here>_1.fastq.gz, <SRR_number_here>_2.fastq.gz. 

```
python3 HyperVar/hypervar-pipeline/scripts/hypervar-pipeline.py \
	--accessions </path/to/INDIVIDUALS>
	--threads 8
	--metawrap 	
	--dastool	
	--magscot	
	--salmon
	--minimap
	--out </path/to/output/directory>
```


For more details, see the link below.

[pipeline tutorial](https://github.com/katemortensen/HyperVar/blob/main/PIPELINE_TUTORIAL.md)


