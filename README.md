# concall

![test-concall-on-mac](https://github.com/zztin/concall/workflows/test-concall-on-mac/badge.svg?branch=master)
![test-concall-ubuntu-linux](https://github.com/zztin/concall/workflows/test-concall-ubuntu-linux/badge.svg)

`concall` is a consensus calling pipeline based on open source deep learning based nanopore consensus sequences polishing "medaka", an experimental project developed and activately maintained by ONT Nanopore. We use the smolecule functionale to polish repeats of RCA product of short cfDNA fragments, aiming to provide high fidelity reads for clinical use. Accurate mutation calling is beneficial for early detection of cancer recurrence.

# Features
- snakemake pipeline management tool takes care of resource distribution and creating independent conda environments for jobs.


# Dependencies
- conda
- snakemake

# Preparation:
1. Reference file:
You need to provide an indexed reference genome to the snakemake pipeline in the configfiles. If you do not have a reference genome yet, you need to download them. For example, download hg19 human reference genome from here:
hg19
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/ 
using 

```
rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet .
gunzip hg19.p13.plusMT.no_alt_analysis_set.fa.gz
bwa index -a bwtsw hg19.p13.plusMT.no_alt_analysis_set.fa
```
This could take up to 2 hours.

2. Backbone files:
A rolling circle amplification product usually contains inserts and backbones. Provide the backbone sequence as fasta at ./concall/data/backbones/. BWA index of the backbone will be created at the same folder for mapping the reads to backbone to check for their existance.

3. Primer files:
Primers files need to be generated from the backbone sequence using ./concall/backbone_processing/get_primers.py. Create 3' and 5' primer with a user defined length, with the option to skip the first few bases if the sequences are too ambiguous.

4. Locate your fastq or fastq.gz files. Each batch can contain up to 40000 reads. If fastq is over this range, batch them with scripts/batch_fastq_gz.sh; scripts/batch_fastq.sh. If the fastq is too small, use scripts/concat_fastq_gz.py to combine them into bigger fastq.
```
scripts/batch_fastq_gz.sh <input_dir> <output_dir>
scripts/concat_fastq.py <input_dir> <output_dir: exact path> <batch-size> <prefix>
```
5. Create a configfile according to your data, and locate it at the ./concall/configfiles/ (example see ./concall/configfiles/config-TESTGITHUB.yaml) Update all parameters. Create a configfile named ./configfiles/config-XXXX.yaml. Elements in a configfile:
```
rawdir: the directory containing the fastq / fastq.gz files
SUP_SAMPLES: Prefix for output folder/files (sample name)
gz: if the fastq files is gzipped
genome: location of the reference genome fa file (make sure the genome has index files in the same folder)
backbone_fa: backbone sequence
min_insert_length: deprecated
max_insert_length:deprecated
3_prime: 3_prime_BBxxx.fa created by ./concall/backbone_processing/get_primers.py.
5_prime: 5_prime_BBxxx.fa created by ./concall/backbone_processing/get_primers.py.
mail: email address
sing: if using singularity image
```

# Usage
Once the preparation is done, concall can be run.

-  run locally

	cd concall/ # Where snakefile is located 
	snakemake --configfile configfiles/config-XXXX.yaml --snakefile tidehunter_only_sing.smk --use-conda --cores 4
	or
	snakemake --configfile configfiles/config-XXXX.yaml--use-conda --use-singularity --snakefile tidehunter_only_sing.smk --cores 4 # at the moment singularity image is stored locally at .sif file which is not included in this repo. Please contact authors for further information. 
	
- run on hpc: slurm
	

	snakemake --profile slurm --jobs 50 --latency-wait 240 --use-conda --use-singularity --rerun-incomplete --keep-going --restart-times 3 --configfile ./test/config/config-xxxx.yaml --snakefile tidehunter_only_sing.smk
	

# Development
1. Tidehunter program can be run locally on a Singularity container (tidehunter.sif) is used as a local file at the moment.
Recipe of the container is in the file Singularity. Currently Singularity Hub uses v2.5 but this image is built and tested on v3.5 machine. It gives errors if image is built by singularity hub and pull down to hpc. (this needs to fix. Before then, please contact the author to get the correct image of Singularity container.)
2. Another workflow make use of ONT Medaka is built in pipeline medaka.smk. medaka has dependency of py=3.6. Environment: create a conda environment and install all package in envs/concall-meta.yaml. Workflow CI still building.

#### Contributions are welcome by raising issues or creating pull request. 


 
