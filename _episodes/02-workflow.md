---
title: "Create a Workflow by Composing Tools"
teaching: 30
exercises: 10
questions:
- "What is the syntax of CWL?"
- "What are the key components of a workflow?"
objectives:
- "Write a workflow based on the source shell script, making use of existing tool wrappers."
keypoints:
- "CWL documents are written using a syntax called YAML."
- "The key components of the workflow are: the header, the inputs, the steps, and the outputs."
---

# Source shell script

In this lesson, we will develop an initial workflow inspired by the
following shell script.

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
{: .language-bash }

# CWL Syntax

CWL documents are written using a format called "YAML".  Here is a crash-course in YAML:

Data fields are written with the name, followed by a colon `:`, a space,
and then the value.

```
fieldName: value
```
{: .language-yaml }

The value is the remaining text to the end of the line.

Special characters in YAML include `:`, `{`, `}` `[`, `]`, `#`, `!`
and `%`.  If your text begins with any of these characters, you must
surround the string in single or double quotes.

```
fieldName: "#quoted-value"
```
{: .language-yaml }

You can write multi-line text by putting `|-` and writing an indented
block.  The leading whitespace will be removed from the actual value.

```
fieldName: |-
  This is a multi-
  line string.
  Horray!
```
{: .language-yaml }

Nested sections are indented:

```
section1:
  field1: value1
  field2: value2
```
{: .language-yaml }

Nested sections can _also_ be wrapped in curly brackets.  In this case, fields must be comma-separated.

```
section1: {field1: value1, field2, value2}
```
{: .language-yaml }

When each item is on its own line starting with a dash `-`, it is a list.

```
section2:
  - value1
  - value2
```
{: .language-yaml }

List can _also_ be wrapped in square brackets.  In this case, values must be comma-separated.

```
section2: [value1, value2]
```
{: .language-yaml }

Comments start with `#`.

```
# This is a comment about field3
field3: stuff

field4: stuff # This is a comment about field4
```
{: .language-yaml }

Finally, YAML is a superset of JSON.  Valid JSON is also valid YAML,
so you may sometimes see JSON format being used instead of YAML format
for CWL documents.

# Workflow header

Create a new file "main.cwl"

Let's start with the header.

```
cwlVersion: v1.2
class: Workflow
label: RNAseq CWL practice workflow
```
{: .language-yaml }

* cwlVersion - Every file must include this.  It declares the version of CWL in use.
* class - This is the type of CWL document.  We will see other types in future lessons.
* label - Optional title of your workflow.

# Workflow Inputs

The purpose of a workflow is to consume some input parameters, run a
series of steps, and produce output values.

For this analysis, the input parameters are the fastq file and the reference data required by STAR.

In the source shell script, the following variables are declared:

```
# initialize a variable with an intuitive name to store the name of the input fastq file
fq=$1

# directory with genome reference FASTA and index files + name of the gene annotation file
genome=rnaseq/reference_data
gtf=rnaseq/reference_data/chr1-hg19_genes.gtf
```
{: .language-bash }

In CWL, we will declare these variables in the `inputs` section.

The inputs section lists each input parameter and its type.  Valid
types include `File`, `Directory`, `string`, `boolean`, `int`, and
`float`.

In this case, the fastq and gene annotation file are individual files.  The STAR index is a directory.  We can describe these inputs in CWL like this:

```
inputs:
  fq: File
  genome: Directory
  gtf: File
```
{: .language-yaml }

# Workflow Steps

A workflow consists of one or more steps.  This is the `steps` section.

Now we need to describe the first step of the workflow.  In the source
script, the first step is to run `fastqc`.

```
# Run FastQC and move output to the appropriate folder
fastqc $fq
```
{: .language-bash }

A workflow step consists of the name of the step, the tool to `run`,
the input parameters to be passed to the tool in `in`, and the output
parameters expected from the tool in `out`.

The value of `run` references the tool file.  The tool file describes
how to run the tool (we will discuss how to write tool files in lesson
4).  If we look in `bio-cwl-tools` (which you should have imported
when setting up a practice repository in the initial setup
instructions) we find `bio-cwl-tools/fastqc/fastqc_2.cwl`.

Next, the `in` block is mapping of input parameters to the tool and
the workflow parameters that will be assigned to those inputs.  We
need to know what input parameters the tool accepts.

Let's open up the tool file and take a look:

Find the `inputs` section of `bio-cwl-tools/fastqc/fastqc_2.cwl`:

```
inputs:

  reads_file:
    type:
      - File
    inputBinding:
      position: 50
    doc: |
      Input bam,sam,bam_mapped,sam_mapped or fastq file
```
{: .language-yaml }

Now we know we need to provide an input parameter called `reads_file`.

Next, the `out` section is a list of output parameters from the tool
that will be used later in the workflow, or as workflow output.  We
need to know what output parameters the tool produces.  Find the
`outputs` section of `bio-cwl-tools/fastqc/fastqc_2.cwl`:

```
outputs:

  zipped_file:
    type:
      - File
    outputBinding:
      glob: '*.zip'
  html_file:
    type:
      - File
    outputBinding:
      glob: '*.html'
  summary_file:
    type:
      - File
    outputBinding:
      glob: |
        ${
          return "*/summary.txt";
        }
```
{: .language-yaml }

Now we know to expect an output parameter called `html_file`.

Putting this all together, the `fastq` step consists of a `run`, `in`
and `out` subsections, and looks like this:

```
steps:
  fastqc:
    run: bio-cwl-tools/fastqc/fastqc_2.cwl
    in:
      reads_file: fq
    out: [html_file]
```
{: .language-yaml }

# Running alignment with STAR

The next step is to run the STAR aligner.

```
# Run STAR
STAR --runThreadN $cores --genomeDir $genome --readFilesIn $fq --outFileNamePrefix $align_out --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard
```
{: .language-bash }

We will go through the same process as the first section.  We find
there is `bio-cwl-tools/STAR/STAR-Align.cwl`.  We will open the file
and look at the `inputs` section to determine what input parameters
correspond to the command line parmeters from our source script.
Command line flags generally appear appear in either the `arguments`
field, or the `prefix` field of the `inputBinding` section of an input
parameter declaration.  For example, this tells us that the
`GenomeDir` input parameter corresponds to the `--genomeDir` command
line parameter.

```
  GenomeDir:
    type: Directory
    inputBinding:
      prefix: "--genomeDir"
```
{: .language-yaml }

Sometimes we want to provide input values to a step without making
them as workflow-level inputs.  We can do this with `{default: N}`.
For example:

```
   in:
     RunThreadN: {default: 4}
```
{: .language-yaml }

> ## `Exercise`
>
> Look at `STAR-Align.cwl` and identify the other input parameters that
> correspond to the command line arguments used in the source script.
> Also identify the output parameter.  Use these to write the STAR
> step.
>
> > ## `Solution`
> >
> > ```
> >  STAR:
> >    run: bio-cwl-tools/STAR/STAR-Align.cwl
> >    in:
> >      RunThreadN: {default: 4}
> >      GenomeDir: genome
> >      ForwardReads: fq
> >      OutSAMtype: {default: BAM}
> >      OutSAMunmapped: {default: Within}
> >    out: [alignment]
> > ```
> > {: .language-yaml }
> {: .solution}
{: .challenge}

# Running samtools

The third step is to generate an index for the aligned BAM.

```
# Create BAM index
samtools index $counts_input_bam
```
{: .language-bash }

For this step, we need to use the output of a previous step as input
to this step.  We refer the output of a step by with name of the step
(STAR), a slash, and the name of the output parameter (alignment), e.g. `STAR/alignment`

This creates a dependency between steps.  This means the `samtools`
step will not run until the `STAR` step has completed successfully.

```
  samtools:
    run: bio-cwl-tools/samtools/samtools_index.cwl
    in:
      bam_sorted: STAR/alignment
    out: [bam_sorted_indexed]
```
{: .language-yaml }

# featureCounts

```
# Count mapped reads
featureCounts -T $cores -s 2 -a $gtf -o $counts $counts_input_bam
```
{: .language-bash }

As of this writing, the `subread` package that provides
`featureCounts` is not available in `bio-cwl-tools` (and if it has been
added since then, let's pretend that it isn't there.)  We will go over
how to write a CWL wrapper for a command line tool in lesson 4.  For
now, we will leave off the final step.

# Workflow Outputs

The last thing to do is declare the workflow outputs in the `outputs` section.

For each output, we need to declare the type of output, and what
parameter has the output value.

Output types are the same as input types, valid types include `File`,
`Directory`, `string`, `boolean`, `int`, and `float`.

The `outputSource` field refers the a step output in the same way that
the `in` block does, the name of the step, a slash, and the name of
the output parameter.

For our final outputs, we want the results from fastqc and the
aligned, sorted and indexed BAM file.

```
outputs:
  qc_html:
    type: File
    outputSource: fastqc/html_file
  bam_sorted_indexed:
    type: File
    outputSource: samtools/bam_sorted_indexed
```
{: .language-yaml }
