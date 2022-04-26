---
title: Setup
redirect: https://doc.arvados.org/rnaseq-cwl-training/setup.html
---

{% capture generic_tab_content %}

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

{% endcapture %}

{% capture arvados_tab_content %}

# Setting up a practice repository

We will create a new git repository and import a library of existing
tool definitions that will help us build our workflow.

When using the recommended [VSCode environment to develop on Arvados](https://doc.arvados.org/v2.3/user/cwl/arvados-vscode-training.html),
start by forking the
[arvados-vscode-cwl-template](https://github.com/arvados/arvados-vscode-cwl-template)
repository.

1. Vscode: On the left sidebar, choose `Explorer` ![](assets/img/Explorer.png)
1. Select `Clone Repository` and enter [https://github.com/arvados/arvados-vscode-cwl-template](https://github.com/arvados/arvados-vscode-cwl-template), then click `Open`
1. If asked `Would you like to open the cloned repository?` choose `Open`

Next, import the [bio-cwl-tools](https://github.com/common-workflow-library/bio-cwl-tools) repository:

1. Vscode: In the top menu, select `Terminal` &rarr; `New Terminal`
1. This will open a terminal window in the lower part of the screen
1. Run this command:
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

1. Go to [https://workbench2.jutro.arvadosapi.com](https://workbench2.jutro.arvadosapi.com) and sign in, this will create an account
2. Go to `Get an API token` under the user menu
3. Log into the shell node of your Arvados cluster
4. On the shell node, copy the host name and token for the `jutro` cluster into the file `~/.config/arvados/jutro.conf` as described on the page for [arv-copy](https://doc.arvados.org/user/topics/arv-copy.html).

Now, on shell node of your Arvados cluster, use `arv-copy` to copy the collection:

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

Use `arv-copy` to copy the collection:

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

{% endcapture %}

<div class="tabbed">
  <ul class="tab">
      <li><a href="#section-generic">generic</a></li>
      <li><a href="#section-arvados">arvados</a></li>
  </ul>

  <section id="section-generic">{{ generic_tab_content | markdownify}}</section>
  <section id="section-arvados">{{ arvados_tab_content | markdownify}}</section>
</div>

{% include links.md %}
