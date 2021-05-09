#!/bin/bash
FILE_FASTA=$1
FILE_SAM=$2
FILE_REF=$3
# hg19: /hpc/compgen/GENOMES/Cyclomics_reference_genome/version12/Homo_sapiens.GRCh37.GATK.illumina_cyclomics_backbone.fasta
# hg38: /hpc/cog_bioinf/GENOMES/1KP_GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa
THREADS=$4
bwa mem -t $THREADS -c 100 -M -R"@RG\\tID:$FILE_SAM\\tSM:$FILE_SAM\\tPL:NANOPORE\\tLB:$FILE_SAM" $FILE_REF $FILE_FASTA > ${FILE_SAM}.sam;
samtools view -h -F 256 ${FILE_SAM}.sam > ${FILE_SAM}.bam;
samtools sort -l 9 ${FILE_SAM}.bam -o ${FILE_SAM}.sorted.bam;
samtools index ${FILE_SAM}.sorted.bam;
# samtools view -F 256 <-- only primary reads
