# Running and debugging a workflow

1. The input parameter file

CWL input values are provided in the form of a YAML or JSON file.
Create one by right clicking on the explorer, select "New File" and
create a called file "main-input.yaml".

This file gives the values for parameters declared in the `inputs`
section of our workflow.  Our workflow takes `fq`, `genome` and `gtf`
as input parameters.

When setting inputs, Files and Directories are given as an object with
`class: File` or `class: Directory`.  This distinguishes them from
plain strings that may or may not be file paths.


## Arvados

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
  location: keep:9178fe1b80a08a422dbe02adfd439764+925/reference_data/chr1-hg19_genes.gtf
```

## Generic

Note: if you don't have example sequence data or the STAR index files, see the Appendix below.

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

2. Running the workflow

## Arvados

In vscode, select "main.cwl" and then choose "Terminal -> Run task -> Run CWL workflow on Arvados"

## Generic

Type this into the terminal:

```
cwl-runner main.cwl main-input.yaml
```

3. Debugging the workflow

A workflow can fail for many reasons: some possible reasons include
bad input, bugs in the code, or running out memory.  In this case, the
STAR workflow might fail with an out of memory error.

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

4. Setting runtime RAM requirements

By default, a step is allocated 256 MB of RAM.  From the STAR error message:

> Check if you have enough RAM 5711762337 bytes

We can see that STAR requires quite a bit more RAM than that.  To
request more RAM, add a "requirements" section with
"ResourceRequirement" to the "STAR" step:

```
  STAR:
    requirements:
      ResourceRequirement:
        ramMin: 8000
    run: bio-cwl-tools/STAR/STAR-Align.cwl
```

Resource requirements you can set include:

* coresMin: CPU cores
* ramMin: RAM (in megabytes)
* tmpdirMin: temporary directory available space
* outdirMin: output directory available space

After setting the RAM requirements, re-run the workflow.

5. Workflow results

The CWL runner will print a results JSON object to standard output.  It will look something like this (it may include additional fields).


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

This has the same structure as `main-input.yaml`.  The each output
parameter is listed, with the `location` field of each `File` object
indicating where the output file can be found.

# Appendix

## Downloading sample and reference data

Start from your rnaseq-cwl-exercises directory.

```
mkdir rnaseq
cd rnaseq
wget --mirror --no-parent --no-host --cut-dirs=1 https://download.pirca.arvadosapi.com/c=9178fe1b80a08a422dbe02adfd439764+925/
```

## Downloading or generating STAR index

Running STAR requires index files generated from the reference.

This is a rather large download (4 GB).  Depending on your bandwidth, it may be faster to generate it yourself.

### Downloading

Go to the "Terminal" tab in the lower vscode panel.  If necessary, select `bash` from the dropdown list in the upper right corner.

```
mkdir hg19-chr1-STAR-index
cd hg19-chr1-STAR-index
wget --mirror --no-parent --no-host --cut-dirs=1 https://download.pirca.arvadosapi.com/c=02a12ce9e2707610991bd29d38796b57+2912/
```

### Generating

Create `chr1-star-index.yaml`:

```
InputFiles:
  - class: File
    location: rnaseq/reference_data/chr1.fa
    format: http://edamontology.org/format_1930
IndexName: 'hg19-chr1-STAR-index'
Gtf:
  class: File
  location: rnaseq/reference_data/chr1-hg19_genes.gtf
Overhang: 99
```

Next, go to the "Terminal" tab in the lower vscode panel.  If
necessary, select `bash` from the dropdown list in the upper right
corner.  Generate the index with your local cwl-runner.

```
cwl-runner bio-cwl-tools/STAR/STAR-Index.cwl chr1-star-index.yaml
```
