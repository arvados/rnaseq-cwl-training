#!/bin/bash

# Based on
# https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/07_automating_workflow.html
#

# This script takes a fastq file of RNA-Seq data, runs FastQC and outputs a counts file for it.
# USAGE: sh rnaseq_analysis_on_input_file.sh <name of fastq file>

set -e

# initialize a variable with an intuitive name to store the name of the input fastq file
fq=$1

# grab base of filename for naming outputs
base=`basename $fq .subset.fq`
echo "Sample name is $base"

# specify the number of cores to use
cores=4

# directory with genome reference FASTA and index files + name of the gene annotation file
genome=rnaseq/reference_data
gtf=rnaseq/reference_data/chr1-hg19_genes.gtf

# make all of the output directories
# The -p option means mkdir will create the whole path if it
# does not exist and refrain from complaining if it does exist
mkdir -p rnaseq/results/fastqc
mkdir -p rnaseq/results/STAR
mkdir -p rnaseq/results/counts

# set up output filenames and locations
fastqc_out=rnaseq/results/fastqc
align_out=rnaseq/results/STAR/${base}_
counts_input_bam=rnaseq/results/STAR/${base}_Aligned.sortedByCoord.out.bam
counts=rnaseq/results/counts/${base}_featurecounts.txt

echo "Processing file $fq"

# Run FastQC and move output to the appropriate folder
fastqc $fq

# Run STAR
STAR --runThreadN $cores --genomeDir $genome --readFilesIn $fq --outFileNamePrefix $align_out --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard

# Create BAM index
samtools index $counts_input_bam

# Count mapped reads
featureCounts -T $cores -s 2 -a $gtf -o $counts $counts_input_bam
