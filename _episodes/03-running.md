---
title: "Running and Debugging a Workflow"
teaching: 15
exercises: 20
questions:
- "How do I provide input to run a workflow?"
- "What should I do if the workflow fails?"
objectives:
- "Write an input parameter file."
- "Execute the workflow."
- "Diagnose workflow errors."
keypoints:
- "The input parameter file is a YAML file with values for each input parameter."
- "A common reason for a workflow step fails is insufficient RAM."
- "Use ResourceRequirement to set the amount of RAM to be allocated to the job."
- "Output parameter values are printed as JSON to standard output at the end of the run."
---

# The input parameter file

CWL input values are provided in the form of a YAML or JSON file.
create a called file

This file gives the values for parameters declared in the `inputs`
section of our workflow.  Our workflow takes `fq`, `genome` and `gtf`
as input parameters.

When setting inputs, Files and Directories are given as an object with
`class: File` or `class: Directory`.  This distinguishes them from
plain strings that may or may not be file paths.

Note: if you don't have example sequence data or the STAR index files, see [setup](/setup.html).

<div>
{% tabs input %}

{% tab input generic %}
main-input.yaml
```
fq:
  class: File
  location: rnaseq/raw_fastq/Mov10_oe_1.subset.fq
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
> Type this into the terminal:
>
> ```
> cwl-runner main.cwl main-input.yaml
> ```
> {: .language-bash }
>
> This may take a few minutes to run, and will print some amount of
> logging.  The logging you see, how access other logs, and how to
> track workflow progress will depend on your CWL runner platform.
{: .challenge }

{% endtab %}

{% tab input arvados %}
main-input.yaml
```
fq:
  class: File
  location: keep:9178fe1b80a08a422dbe02adfd439764+925/raw_fastq/Mov10_oe_1.subset.fq
  format: http://edamontology.org/format_1930
genome:
  class: Directory
  location: keep:02a12ce9e2707610991bd29d38796b57+2912
gtf:
  class: File
  location: 9178fe1b80a08a422dbe02adfd439764+925/reference_data/chr1-hg19_genes.gtf
```
{: .language-yaml }

> ## Running the workflow
>
> If you are using VSCode with Arvados tasks, select `main.cwl` and
> then use the `Run CWL Workflow on Arvados` task.
>
{: .challenge }
{% endtab %}
{% endtabs %}
</div>

# Debugging the workflow

Depending on whether and how your workflow platform enforces memory
limits, your workflow may fail.  Let's talk about what to do when a
workflow fails.

A workflow can fail for many reasons: some possible reasons include
bad input, bugs in the code, or running out memory.  In our example,
the STAR workflow may fail with an out of memory error.

To help diagnose these errors, the workflow runner produces logs that
record what happened, either in the terminal or the web interface.

Some errors you might see in the logs that would indicate an out of
memory condition:

```
EXITING: fatal error trying to allocate genome arrays, exception thrown: std::bad_alloc
Possible cause 1: not enough RAM. Check if you have enough RAM 5711762337 bytes
Possible cause 2: not enough virtual memory allowed with ulimit. SOLUTION: run ulimit -v 5711762337
```

or

```
Container exited with code: 137
```

(Exit code 137 most commonly occurs when a container goes "out of memory" and is terminated by the operating system).

If this happens, you will need to request more RAM.

# Setting runtime RAM requirements

By default, a step is allocated 256 MB of RAM.  From the STAR error message:

> Check if you have enough RAM 5711762337 bytes

We can see that STAR requires quite a bit more RAM than 256 MB.  To
request more RAM, add a "requirements" section with
"ResourceRequirement" to the "STAR" step:

```
  STAR:
    requirements:
      ResourceRequirement:
        ramMin: 9000
    run: bio-cwl-tools/STAR/STAR-Align.cwl
	...
```
{: .language-yaml }

Resource requirements you can set include:

* coresMin: CPU cores
* ramMin: RAM (in megabytes)
* tmpdirMin: temporary directory available space
* outdirMin: output directory available space

> ## Running the workflow
>
> Now that you've fixed the workflow, run it again.
>
{: .challenge }

> ## Episode solution
> * <a href="{% link assets/answers/ep3/main.cwl %}">main.cwl</a>
{: .solution}

# Workflow results

The CWL runner will print a results JSON object to standard output.  It will look something like this (it may include additional fields).

<div>
{% tabs output %}

{% tab output generic %}
```
{
    "bam_sorted_indexed": {
        "location": "file:///home/username/rnaseq-cwl-training-exercises/Aligned.sortedByCoord.out.bam",
        "basename": "Aligned.sortedByCoord.out.bam",
        "class": "File",
        "size": 25370707,
        "secondaryFiles": [
            {
                "basename": "Aligned.sortedByCoord.out.bam.bai",
                "location": "file:///home/username/rnaseq-cwl-training-exercises/Aligned.sortedByCoord.out.bam.bai",
                "class": "File",
                "size": 176552,
            }
        ]
    },
    "qc_html": {
        "location": "file:///home/username/rnaseq-cwl-training-exercises/Mov10_oe_1.subset_fastqc.html",
        "basename": "Mov10_oe_1.subset_fastqc.html",
        "class": "File",
        "size": 383589
    }
}
```
{: .language-yaml }
{% endtab %}

{% tab output arvados %}
```
{
    "bam_sorted_indexed": {
        "basename": "Aligned.sortedByCoord.out.bam",
        "class": "File",
        "location": "keep:2dbaaef5aefd558e37f14280e47091a9+327/Aligned.sortedByCoord.out.bam",
        "secondaryFiles": [
            {
                "basename": "Aligned.sortedByCoord.out.bam.bai",
                "class": "File",
                "location": "keep:2dbaaef5aefd558e37f14280e47091a9+327/Aligned.sortedByCoord.out.bam.bai"
            }
        ],
        "size": 25370695
    },
    "qc_html": {
        "basename": "Mov10_oe_1.subset_fastqc.html",
        "class": "File",
        "location": "keep:2dbaaef5aefd558e37f14280e47091a9+327/Mov10_oe_1.subset_fastqc.html",
        "size": 383589
    }
}
```
{: .language-yaml }
{% endtab %}
{% endtabs %}
</div>

This has a similar structure as `main-input.yaml`.  The each output
parameter is listed, with the `location` field of each `File` object
indicating where the output file can be found.
