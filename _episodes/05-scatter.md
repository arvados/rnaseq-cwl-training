---
title: " Analyzing multiple samples"
teaching: 0
exercises: 0
questions:
- "Key question (FIXME)"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

Analyzing a single sample is great, but in the real world you probably
have a batch of samples that you need to analyze and then compare.

# 1. Subworkflows

In addition to running command line tools, a workflow step can also
execute another workflow.

Let's copy "main.cwl" to "alignment.cwl".

Now, edit open "main.cwl" for editing.  We are going to replace the `steps` and `outputs` sections.

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

In the outputs section, all the output sources are from the alignment step:

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

We also need a little boilerplate to tell the workflow runner that we want to use subworkflows:

```
requirements:
  SubworkflowFeatureRequirement: {}
```

If you run this workflow, you will get exactly the same results as
before, we've just wrapped the inner workflow with an outer workflow.

# 2. Scattering

The wrapper lets us do something useful.  We can modify the outer
workflow to accept a list of files, and then invoke the inner workflow
step for every one of those files.  We will need to modify the
`inputs`, `steps`, `outputs`, and `requirements` sections.

First we change the `fq` parameter to expect a list of files:

```
inputs:
  fq: File[]
  genome: Directory
  gtf: File
```

Next, we add `scatter` to the alignment step.  The means it will
run `alignment.cwl` for each value in the list in the `fq` parameter.

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

Finally, we need a little more boilerplate to tell the workflow runner
that we want to use scatter:

```
requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
```

# 3. Running with list inputs

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

Now you can run the workflow the same way as in Lesson 2.

# 4. Combining results

Each instance of the alignment workflow produces its own featureCounts
file.  However, to be able to compare results easily, we need them a
single file with all the results.

The easiest way to do this is to run `featureCounts` just once at the
end of the workflow, with all the bam files listed on the command
line.

We'll need to modify a few things.

First, in `featureCounts.cwl` we need to modify it to accept either a
single bam file or list of bam files.

```
inputs:
  gtf: File
  counts_input_bam:
   - File
   - File[]
```

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

Run this workflow to get a single `featurecounts.tsv` file with a column for each bam file.
