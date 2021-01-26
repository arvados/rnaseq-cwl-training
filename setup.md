---
title: Setup
---

# Setting up a practice repository

We will create a new git repository and import a library of existing
tool definitions that will help us build our workflow.

Create a new git repository to hold our workflow with this command:

```
git init rnaseq-cwl-training-exercises
```

On Arvados use this:

```
git clone https://github.com/arvados/arvados-vscode-cwl-template.git rnaseq-cwl-training-exercises
```

Next, import bio-cwl-tools with this command:

```
git submodule add https://github.com/common-workflow-library/bio-cwl-tools.git
```

# Downloading sample and reference data

Start from your rnaseq-cwl-exercises directory.

```
mkdir rnaseq
cd rnaseq
wget --mirror --no-parent --no-host --cut-dirs=1 https://download.pirca.arvadosapi.com/c=9178fe1b80a08a422dbe02adfd439764+925/
```

# Downloading or generating STAR index

Running STAR requires index files generated from the reference.

This is a rather large download (4 GB).  Depending on your bandwidth, it may be faster to generate it yourself.

## Downloading

```
mkdir hg19-chr1-STAR-index
cd hg19-chr1-STAR-index
wget --mirror --no-parent --no-host --cut-dirs=1 https://download.pirca.arvadosapi.com/c=02a12ce9e2707610991bd29d38796b57+2912/
```

## Generating

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

Generate the index with your local cwl-runner.

```
cwl-runner bio-cwl-tools/STAR/STAR-Index.cwl chr1-star-index.yaml
```


{% include links.md %}
