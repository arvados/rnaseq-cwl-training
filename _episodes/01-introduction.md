---
title: "Introduction"
teaching: 0
exercises: 0
questions:
- "Key question (FIXME)"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

## Introduction

The goal of this training is to walk through the development of a
best-practices CWL workflow by translating an existing bioinformatics
shell script into CWL.  Specific knowledge of the biology of RNA-seq
is *not* a prerequisite for these lessons.

These lessons are based on "Introduction to RNA-seq using
high-performance computing (HPC)" lessons developed by members of the
teaching team at the Harvard Chan Bioinformatics Core (HBC).  The
original training, which includes additional lectures about the
biology of RNA-seq can be found here:

https://github.com/hbctraining/Intro-to-rnaseq-hpc-O2

## Background

RNA-seq is the process of sequencing RNA in a biological sample.  From
the sequence reads, we want to measure the relative number of RNA
molecules appearing in the sample that were produced by particular
genes.  This analysis is called "differential gene expression".

The entire process looks like this:

![](/assets/img/RNAseqWorkflow.png)

For this training, we are only concerned with the middle analytical
steps (skipping adapter trimming).

* Quality control (FASTQC)
* Alignment (mapping)
* Counting reads associated with genes

## Analysis shell script

This analysis is already available as a Unix shell script, which we
will refer to in order to build the workflow.

Some of the reasons to use CWL over a plain shell script: portability,
scalability, ability to run on platforms that are not traditional HPC.

rnaseq_analysis_on_input_file.sh

```
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
```


{% include links.md %}
