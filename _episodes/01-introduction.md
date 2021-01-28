---
title: "Introduction"
teaching: 10
exercises: 0
questions:
- "What is CWL?"
- "What are the requirements for this training?"
- "What is the goal of this training?"
objectives:
- "Understand how the training will be motivated by an example analysis."
keypoints:
- "Common Workflow Language is a standard for describing data analysis workflows"
- "This training assumes some basic familiarity with editing text files, the Unix command line, and Unix shell scripts."
- "We will use an bioinformatics RNA-seq analysis as an example workflow, but does not require in-depth knowledge of biology."
- "After completing this training, you should be able to begin writing workflows for your own analysis, and know where to learn more."
---

# Introduction to Common Worklow Language

The Common Workflow Language (CWL) is an open standard for describing
automated, batch data analysis workflows.  Unlike many programming
languages, CWL is a declarative language.  This means it describes
_what_ should happen, but not _how_ it should happen.  This enables
workflows written in CWL to be portable and scalable across a variety
of software and hardware environments, from workstations to cluster,
cloud, and high performance computing (HPC) environments.  As a
standard with multiple implementations, CWL is particularly well
suited for research collaboration, publishing, and high-throughput
production data analysis.

# Introduction to this training

The goal of this training is to walk the student through the
development of a best-practices CWL workflow, starting from an
existing shell script that performs a simple RNA-seq bioinformatics
analysis.  At the conclusion of this training, you should have a grasp
of the essential components of a workflow, and have a basis for
learning more.

This training assumes some basic familiarity with editing text files,
the Unix command line, and Unix shell scripts.

Specific knowledge of the biology of RNA-seq is *not* a prerequisite
for these lessons.  Although orignally developed to solve big data
problems in genomics, CWL is not domain specific to bioinformatics,
and is used in a number of other fields including medical imaging,
astronomy, geospatial, and machine learning.  We hope that you will
find this training useful regardless of your area of research.

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
scratch, instead we we will be starting from an analysis already
written in a shell script, which will be supplied in lesson 2.

{% include links.md %}
