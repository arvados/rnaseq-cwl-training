# Writing a tool wrapper

It is time to add the last step in the analysis.

This will use the "featureCounts" tool from the "subread" package.

# Writing the tool wrapper

1. Create a new file "featureCounts.cwl"

2. Start with this header

```
cwlVersion: v1.2
class: CommandLineTool
```

3. Command line tool inputs

A CommandLineTool describes a single invocation of a command line program.

It consumes some input parameters, runs a program, and produce output
values.

Here is the original shell command:

```
featureCounts -T $cores -s 2 -a $gtf -o $counts $counts_input_bam
```

The variables used in the bash script are `$cores`, `$gtf`, `$counts` and `$counts_input_bam`.

The parameters

This gives us two file inputs, `gtf` and `counts_input_bam` which we can declare in our `inputs` section:

```
inputs:
  gtf: File
  counts_input_bam: File
```

4. Specifying the program to run

Give the name of the program to run in `baseCommand`.

```
baseCommand: featureCounts
```

5. Command arguments

The easiest way to describe the command line is with an `arguments`
section.  This takes a comma-separated list of command line arguments.

Input variables are included on the command line as
`$(inputs.name_of_parameter)`.  When the tool is executed, these input
parameter values are substituted for these variable.

Special variables are also available.  The runtime environment
describes the resources allocated to running the program.  Here we use
`$(runtime.cores)` to decide how many threads to request.

```
arguments: [-T, $(runtime.cores),
            -a, $(inputs.gtf),
			-o, featurecounts.tsv,
			$(inputs.counts_input_bam)]
```

6. Outputs section

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

7. Running in a container

In order to run the tool, it needs to be installed.
Using software containers, a tool can be pre-installed into a
compatible runtime environment, and that runtime environment (called a
container image) can be downloaded and run on demand.

Many bioinformatics tools are already available as containers.  One
resource is the BioContainers project.  Let's find the "subread" software:

   1. Visit https://biocontainers.pro/
   2. Click on "Registry"
   3. Search for "subread"
   4. Click on the search result for "subread"
   5. Click on the tab "Packages and Containers"
   6. Choose a row with type "docker", then on the right side of the "Full
Tag" column for that row, click the "copy to clipboard" button.

To declare that you want to run inside a container, create a section
called `hints` with a subsection `DockerRequirement`.  Under
`DockerRequirement`, paste the text your copied in the above step.
Replace the text `docker pull` to `dockerPull:` and indent it so it is
in the `DockerRequirement` section.

```
hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/subread:1.5.0p3--0
```

8. Running a tool on its own

When creating a tool wrapper, it is helpful to run it on its own to test it.

The input to a single tool is the same kind of input parameters file
that we used as input to a workflow in the previous lesson.

featureCounts.yaml:

```
counts_input_bam:
  class: File
  location: Aligned.sortedByCoord.out.bam
gtf:
  class: File
  location: rnaseq/reference_data/chr1-hg19_genes.gtf
```

The invocation is also the same:

```
cwl-runner featureCounts.cwl featureCounts.yaml
```

9. Adding it to the workflow

Now that we have confirmed that it works, we can add it to our workflow.
We add it to `steps`, connecting the output of samtools to
`counts_input_bam` and the `gtf` taking the workflow input of the same
name.

```
steps:
  ...
  featureCounts:
    requirements:
      ResourceRequirement:
        ramMin: 500
    run: featureCounts.cwl
	in:
      counts_input_bam: samtools/bam_sorted_indexed
	  gtf: gtf
	out: [featurecounts]
```

We will add the result from featurecounts to the output:

```
outputs:
  ...
  featurecounts:
    type: File
	outputSource: featureCounts/featurecounts

```

You should now be able to re-run the workflow and it will run the
"featureCounts" step and include "featurecounts" in the output.
