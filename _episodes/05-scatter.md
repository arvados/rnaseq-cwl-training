---
title: "Analyzing Multiple Samples"
teaching: 30
exercises: 30
questions:
- "How can you run the same workflow over multiple samples?"
objectives:
- "Modify the workflow to process multiple samples, then perform a joint analysis."
keypoints:
- "Separate the part of the workflow that you want to run multiple times into a subworkflow."
- "Use a scatter step to run the subworkflow over a list of inputs."
- "The result of a scatter is an array, which can be used in a combine step to get a single result."
---

In the previous lesson, we completed converting the function of the
original source shell script into CWL.  This lesson expands the scope
by demonstrating what changes to make to the workflow to be able to
analyze multiple samples in parallel.

# Subworkflows

In addition to running command line tools, a workflow step can also
execute another workflow.

First, copy `main.cwl` to `alignment.cwl`.

Next, open `main.cwl` for editing.  We are going to replace the `steps` and `outputs` sections.

Remove all the steps and replace them with a single `alignment` step
which invokes the `alignment.cwl` we just copied.

```
steps:
  alignment:
    run: alignment.cwl
    in:
      fq: fq
      genome: genome
      gtf: gtf
    out: [qc_html, bam_sorted_indexed, featurecounts]
```
{: .language-yaml }

In the `outputs` section, all the output sources are from the alignment step:

```
outputs:
  qc_html:
    type: File
    outputSource: alignment/qc_html
  bam_sorted_indexed:
    type: File
    outputSource: alignment/bam_sorted_indexed
  featurecounts:
    type: File
    outputSource: alignment/featurecounts
```
{: .language-yaml }

We also need add "SubworkflowFeatureRequirement" to tell the workflow
runner that we are using subworkflows:

```
requirements:
  SubworkflowFeatureRequirement: {}
```
{: .language-yaml }

> ## Running the workflow
>
> Run this workflow.  You should get exactly the same results as
> before, as all we have done so far is to wrap the inner workflow with
> an outer workflow.
>
{: .challenge }

# Scattering

The "wrapper" step lets us do something useful.  We can modify the
outer workflow to accept a list of files, and then invoke the inner
workflow step for every one of those files.  We will need to modify
the `inputs`, `steps`, `outputs`, and `requirements` sections.

First we change the `fq` parameter to expect a list of files:

```
inputs:
  fq: File[]
  genome: Directory
  gtf: File
```
{: .language-yaml }

Next, we add `scatter` to the alignment step.  The means we want to
run run `alignment.cwl` for each value in the list in the `fq`
parameter.

```
steps:
  alignment:
    run: alignment.cwl
    scatter: fq
    in:
      fq: fq
      genome: genome
      gtf: gtf
    out: [qc_html, bam_sorted_indexed, featurecounts]
```
{: .language-yaml }

Because the scatter produces multiple outputs, each output parameter
becomes a list as well:

```
outputs:
  qc_html:
    type: File[]
    outputSource: alignment/qc_html
  bam_sorted_indexed:
    type: File[]
    outputSource: alignment/bam_sorted_indexed
  featurecounts:
    type: File[]
    outputSource: alignment/featurecounts
```
{: .language-yaml }

We also need add "ScatterFeatureRequirement" to tell the workflow
runner that we are using scatter:

```
requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
```
{: .language-yaml }

# Input parameter lists

The `fq` parameter needs to be a list.  You write a list in yaml by
starting each list item with a dash.  Example `main-input.yaml`

```
fq:
  - class: File
    location: rnaseq/raw_fastq/Mov10_oe_1.subset.fq
    format: http://edamontology.org/format_1930
  - class: File
    location: rnaseq/raw_fastq/Mov10_oe_2.subset.fq
    format: http://edamontology.org/format_1930
  - class: File
    location: rnaseq/raw_fastq/Mov10_oe_3.subset.fq
    format: http://edamontology.org/format_1930
  - class: File
    location: rnaseq/raw_fastq/Irrel_kd_1.subset.fq
    format: http://edamontology.org/format_1930
  - class: File
    location: rnaseq/raw_fastq/Irrel_kd_2.subset.fq
    format: http://edamontology.org/format_1930
  - class: File
    location: rnaseq/raw_fastq/Irrel_kd_3.subset.fq
    format: http://edamontology.org/format_1930
genome:
  class: Directory
  location: hg19-chr1-STAR-index
gtf:
  class: File
  location: rnaseq/reference_data/chr1-hg19_genes.gtf
```
{: .language-yaml }

> ## Running the workflow
>
> Run this workflow.  You will now get results for each one of the
> input fastq files.
>
{: .challenge }

# Combining results

Each instance of the alignment workflow produces its own
`featurecounts.tsv` file.  However, to be able to compare results
easily, we would like single file with all the results.

We can modify the workflow to run `featureCounts` once at the end of
the workflow, taking all the bam files listed on the command line.

We will need to change a few things.

First, in `featureCounts.cwl` we need to modify it to accept either a
single bam file or list of bam files.

```
inputs:
  gtf: File
  counts_input_bam:
   - File
   - File[]
```
{: .language-yaml }

Second, in `alignment.cwl` we need to remove the `featureCounts` step from alignment.cwl, as well as the `featurecounts` output parameter.

Third, in `main.cwl` we need to remove `featurecounts` from the `alignment` step
outputs, and add a new step:

```
steps:
  alignment:
    run: alignment.cwl
    scatter: fq
    in:
      fq: fq
      genome: genome
      gtf: gtf
    out: [qc_html, bam_sorted_indexed]
  featureCounts:
    requirements:
      ResourceRequirement:
        ramMin: 500
    run: featureCounts.cwl
    in:
      counts_input_bam: alignment/bam_sorted_indexed
      gtf: gtf
    out: [featurecounts]
```
{: .language-yaml }

Last, we modify the `featurecounts` output parameter.  Instead of a
list of files produced by the `alignment` step, it is now a single
file produced by the new `featureCounts` step.

```
outputs:
  ...
  featurecounts:
    type: File
    outputSource: featureCounts/featurecounts
```
{: .language-yaml }

> ## Running the workflow
>
> Run this workflow.  You will still have separate results from fastq
> and and STAR, but now you will only have a single
> `featurecounts.tsv` file with a column for each bam file.
>
{: .challenge }
