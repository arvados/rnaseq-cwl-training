---
title: "Writing a Tool Wrapper"
teaching: 20
exercises: 30
questions:
- "What are the key components of a tool wrapper?"
- "How do I use software containers to supply the software I want to run?"
objectives:
- "Write a tool wrapper for the featureCounts tool."
- "Find an software container that has the software we want to use."
- "Add the tool wrapper to our main workflow."
keypoints:
- "The key components of a command line tool wrapper are the header, inputs, baseCommand, arguments, and outputs."
- "Like workflows, CommandLineTools have `inputs` and `outputs`."
- "Use `baseCommand` and `arguments` to provide the program to run and the command line arguments to run it with."
- "Use `glob` to capture output files and assign them to output parameters."
- "Use DockerRequirement to supply the name of the Docker image that contains the software to run."
---

It is time to add the last step in the analysis.

```
# Count mapped reads
featureCounts -T $cores -s 2 -a $gtf -o $counts $counts_input_bam
```
{: .language-bash }

This will use the "featureCounts" tool from the "subread" package.

# File header

A CommandLineTool describes a single invocation of a command line
program.  It consumes some input parameters, runs a program, and
captures output, mainly in in the form of files produced by the
program.

Create a new file "featureCounts.cwl"

Let's start with the header.  This is very similar to the workflow, except that we use `class: CommandLineTool`.

```
cwlVersion: v1.2
class: CommandLineTool
label: featureCounts tool
```
{: .language-yaml }

# Command line tool inputs

The `inputs` section describes input parameters with the same form as
the Workflow `inputs` section.

> ## Exercise
>
> The variables used in the bash script are `$cores`, `$gtf`, `$counts` and `$counts_input_bam`.
>
> * $cores is the number of CPU cores to use.
> * $gtf is the input .gtf file
> * $counts is the name we will give to the output file
> * $counts_input_bam is the input .bam file
>
> Write the `inputs` section for the File inputs `gtf` and `counts_input_bam`.
>
> > ## Solution
> > ```
> > inputs:
> >   gtf: File
> >   counts_input_bam: File
> > ```
> > {: .language-yaml }
> {: .solution}
{: .challenge}

# Specifying the program to run

Give the name of the program to run in `baseCommand`.

```
baseCommand: featureCounts
```
{: .language-yaml }

# Command arguments

The easiest way to describe the command line is with an `arguments`
section.  This takes a comma-separated list of command line arguments.


```
arguments: [-T, $(runtime.cores),
            -a, $(inputs.gtf),
            -o, featurecounts.tsv,
            $(inputs.counts_input_bam)]
```
{: .language-yaml }

Input variables are included on the command line as
`$(inputs.name_of_parameter)`.  When the tool is executed, the
variables will be replaced with the input parameter values.

There are also some special variables.  The `runtime` object describes
the resources allocated to running the program.  Here we use
`$(runtime.cores)` to decide how many threads to request.

> ## `arguments` vs `inputBinding`
>
> You may recall from examining existing the fastqc and STAR tools
> wrappers in lesson 2, another way to express command line parameters
> is with `inputBinding` and `prefix` on individual input parameters.
>
> ```
> inputs:
>   parametername:
>     type: parametertype
>     inputBinding:
>       prefix: --some-option
> ```
> {: .language-yaml }
>
> We use `arguments` in the example simply because it is easier to see
> how it lines up with the source shell script.
>
> You can use both `inputBinding` and `arguments` in the same
> CommandLineTool document.  There is no "right" or "wrong" way, and
> one does not override the other, they are combined to produce the
> final command line invocation.
>
{: .callout}

# Outputs section

In CWL, you must explicitly identify the outputs of a program.  This
associates output parameters with specific files, and enables the
workflow runner to know which files must be saved and which files can
be discarded.

In the previous section, we told the featureCounts program the name of
our output files should be `featurecounts.tsv`.

We can declare an output parameter called `featurecounts` that will
have that output file as its value.

The `outputBinding` section describes how to determine the value of
the parameter.  The `glob` field tells it to search for a file in the
output directory called `featurecounts.tsv`

```
outputs:
  featurecounts:
    type: File
    outputBinding:
      glob: featurecounts.tsv
```
{: .language-yaml }

# Running in a container

In order to run the tool, it needs to be installed.
Using software containers, a tool can be pre-installed into a
compatible runtime environment, and that runtime environment (called a
container image) can be downloaded and run on demand.

Although plain CWL does not _require_ the use of containers, many
popular platforms that run CWL do require the software be supplied in
the form of a container image.

> ## Finding container images
>
> Many bioinformatics tools are already available as containers.  One
> resource is the BioContainers project.  Let's find the "subread" software:
>
>   1. Visit [https://biocontainers.pro/](https://biocontainers.pro/)
>   2. Click on "Registry"
>   3. Search for "subread"
>   4. Click on the search result for "subread"
>   5. Click on the tab "Packages and Containers"
>   6. Choose a row with type "docker", then on the right side of the "Full
> Tag" column for that row, click the "copy to clipboard" button.
>
> To declare that you want to run inside a container, add a section
> called `hints` to your tool document.  Under `hints` add a
> subsection `DockerRequirement`.  Under `DockerRequirement`, paste
> the text your copied in the above step.  Replace the text `docker
> pull` to `dockerPull:` ensure it is indented twice so it is a field
> of `DockerRequirement`.
>
> > ## Answer
> > ```
> > hints:
> >   DockerRequirement:
> >     dockerPull: quay.io/biocontainers/subread:1.5.0p3--0
> > ```
> > {: .language-yaml }
> {: .solution}
{: .challenge}

# Running a tool on its own

When creating a tool wrapper, it is helpful to run it on its own to test it.

The input to a single tool is the same kind of input parameters file
that we used as input to a workflow in the previous lesson.

`featureCounts.yaml`

```
counts_input_bam:
  class: File
  location: Aligned.sortedByCoord.out.bam
gtf:
  class: File
  location: rnaseq/reference_data/chr1-hg19_genes.gtf
```
{: .language-yaml }

> ## Running the tool
>
> Run the tool on its own to confirm it has correct behavior:
>
> ```
> cwl-runner featureCounts.cwl featureCounts.yaml
> ```
> {: .language-bash }
{: .challenge }

# Adding it to the workflow

Now that we have confirmed that the tool wrapper works, it is time to
add it to our workflow.

> ## Exercise
>
>   1. Add a new step called `featureCounts` that runs our tool
>      wrapper.  The new step should take input from
>      `samtools/bam_sorted_indexed`, and should be allocated a
>      minimum of 500 MB of RAM
>   2. Add a new output parameter for the workflow called
>      `featurecounts` The output source should come from the output
>      of the new `featureCounts` step.
>   3.  When you have an answer, run the updated workflow, which
>       should run the "featureCounts" step and produce "featurecounts"
>       output parameter.
>
> > ## Answer
> > ```
> > steps:
> >   ...
> >   featureCounts:
> >     requirements:
> >       ResourceRequirement:
> >         ramMin: 500
> >     run: featureCounts.cwl
> >     in:
> >       counts_input_bam: samtools/bam_sorted_indexed
> >       gtf: gtf
> >     out: [featurecounts]
> >
> > outputs:
> >   ...
> >   featurecounts:
> >     type: File
> >     outputSource: featureCounts/featurecounts
> > ```
> > {: .language-yaml }
> {: .solution}
{: .challenge}
