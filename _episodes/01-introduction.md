---
title: "Introduction"
teaching: 10
exercises: 0
questions:
- "What is CWL?"
- "What is the goal of this training?"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

# Introduction to Common Worklow Language

The Common Workflow Language (CWL) is an open standard for describing
analysis workflows and tools in a way that makes them portable and
scalable across a variety of software and hardware environments, from
workstations to cluster, cloud, and high performance computing (HPC)
environments. CWL is designed to meet the needs of data-intensive
science, such as Bioinformatics, Medical Imaging, Astronomy, High
Energy Physics, and Machine Learning.

# Introduction to this training

The goal of this training is to walk the student through the
development of a best-practices CWL workflow, starting from an
existing shell script that performs a common bioinformatics analysis.

Specific knowledge of the biology of RNA-seq is *not* a prerequisite
for these lessons.  CWL is not domain specific to bioinformatics.  We
hope that you will find this training useful even if you work in some
other field of research.

These lessons are based on [Introduction to RNA-seq using
high-performance computing
(HPC)](https://github.com/hbctraining/Intro-to-rnaseq-hpc-O2) lessons
developed by members of the teaching team at the Harvard Chan
Bioinformatics Core (HBC).  The original training, which includes
additional lectures about the biology of RNA-seq, can be found at that
link.

# Introduction to the example analysis

RNA-seq is the process of sequencing RNA present in a biological
sample.  From the sequence reads, we want to measure the relative
numbers of different RNA molecules appearing in the sample that were
produced by particular genes.  This analysis is called "differential
gene expression".

The entire process looks like this:

![](/assets/img/RNAseqWorkflow.png){: height="400px"}

For this training, we are only concerned with the middle analytical
steps (skipping adapter trimming).

* Quality control (FASTQC)
* Alignment (mapping)
* Counting reads associated with genes

In this training, we are not attempting to develop the analysis from
scratch, instead we we will be starting from an analysis written as a
shell script.  We will be using the following shell script as a guide to build
our workflow.

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
