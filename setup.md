---
title: Setup
---

<div>
{% tabs setup %}
{% tab setup generic %}

# Setting up a practice repository

We will create a new git repository and import a library of existing
tool definitions that will help us build our workflow.

Create a new empty git repository to hold our workflow with this command:

```
git init rnaseq-cwl-training-exercises
```
{: .language-bash }

Next, import bio-cwl-tools with this command:

```
git submodule add https://github.com/common-workflow-library/bio-cwl-tools.git
```
{: .language-bash }

# Downloading sample and reference data

Start from your rnaseq-cwl-exercises directory.

```
mkdir rnaseq
cd rnaseq
wget --mirror --no-parent --no-host --cut-dirs=1 https://download.jutro.arvadosapi.com/c=9178fe1b80a08a422dbe02adfd439764+925/
```
{: .language-bash }

# Downloading or generating STAR index

Running STAR requires index files generated from the reference.

This is a rather large download (4 GB).  Depending on your bandwidth, it may be faster to generate it yourself.

## Downloading

```
mkdir hg19-chr1-STAR-index
cd hg19-chr1-STAR-index
wget --mirror --no-parent --no-host --cut-dirs=1 https://download.jutro.arvadosapi.com/c=02a12ce9e2707610991bd29d38796b57+2912/
```
{: .language-bash }

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
{: .language-yaml }

Generate the index with your local cwl-runner.

```
cwl-runner bio-cwl-tools/STAR/STAR-Index.cwl chr1-star-index.yaml
```
{: .language-bash }


{% endtab %}

{% tab setup arvados %}

# Setting up a practice repository

We will create a new git repository and import a library of existing
tool definitions that will help us build our workflow.

When using the recommended VSCode environment to develop on Arvados, start by forking this repository:
```
git clone https://github.com/arvados/arvados-vscode-cwl-template.git rnaseq-cwl-training-exercises
```
{: .language-bash }

Next, import bio-cwl-tools with this command:

```
git submodule add https://github.com/common-workflow-library/bio-cwl-tools.git
```
{: .language-bash }

# Downloading sample and reference data

> ## Note
>
> You may already have access to this collection.
>
> You can check by going to Workbench and pasting
> `9178fe1b80a08a422dbe02adfd439764+925` into the search box.  If you
> arrived at a collection page instead of a "not found" error, then
> you do not need to perform this download step.
{: .callout}

```
arv-copy --src jutro 9178fe1b80a08a422dbe02adfd439764+925
```
{: .language-bash }

# Downloading or generating STAR index

Running STAR requires index files generated from the reference.

This is a rather large download (4 GB).  Depending on your bandwidth, it may be faster to generate it yourself.

## Downloading

> ## Note
>
> As above, you can check by going to Workbench and pasting
> `02a12ce9e2707610991bd29d38796b57+2912` into the search box to see
> if you already have access to this collection.
{: .callout}

```
arv-copy --src jutro 02a12ce9e2707610991bd29d38796b57+2912
```
{: .language-bash }

## Generating

Create `chr1-star-index.yaml`:

```
InputFiles:
  - class: File
    location: keep:9178fe1b80a08a422dbe02adfd439764+925/reference_data/chr1.fa
    format: http://edamontology.org/format_1930
IndexName: 'hg19-chr1-STAR-index'
Gtf:
  class: File
  location: keep:9178fe1b80a08a422dbe02adfd439764+925/reference_data/chr1-hg19_genes.gtf
Overhang: 99
```
{: .language-yaml }

Generate the index with arvados-cwl-runner.

```
arvados-cwl-runner bio-cwl-tools/STAR/STAR-Index.cwl chr1-star-index.yaml
```
{: .language-bash }

{% endtab %}
{% endtabs %}
</div>

{% include links.md %}
