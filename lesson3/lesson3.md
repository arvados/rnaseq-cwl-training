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

Here's the original bash script

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

4. The base command

This one is easy.  This is the name of program to run:

```
baseCommand: featureCounts
```

5. The command arguments

The easiest way to describe the command line is with an `arguments`
section.  This takes a comma-separated list of command line arguments.

Input variables are included on the command line as
`$(inputs.name_of_parameter)`.  When the tool is executed, these input
parameter values are substituted for these variable.

Special variables are also available.  The runtime environment
describes the resources allocated to running the program.  Here we use
`$(runtime.cores)` to decide how many threads to request.

File variables can also yield a partial filename, by adding
`.nameroot`.  This is the filename with the final dot-extension
stripped off.

```
arguments: [-T, $(runtime.cores),
            -a, $(inputs.gtf),
			-o, $(inputs.counts_input_bam.nameroot)_featurecounts.txt,
			$(inputs.counts_input_bam)]
```

6. The outputs section

In CWL, you must explicitly identify the outputs of a program.  This
associates output parameters with specific files, and allows the
workflow runner to know which files must be saved and which files can
be discarded.

In the previous section, we told the featureCounts program the name of
our output files should be
`$(inputs.counts_input_bam.nameroot)_featurecounts.txt`.

We can declare an output parameter called `featurecounts` that will
have that output file as its value.

The `outputBinding` section describes how to determine the value of
the parameter.  The `glob` field tells it to search for a file in the
output directory with the
`$(inputs.counts_input_bam.nameroot)_featurecounts.txt`

```
outputs:
  featurecounts:
    type: File
	outputBinding:
	  glob: $(inputs.counts_input_bam.nameroot)_featurecounts.txt
```

N.

The most portable way to run a tool is to wrap it in a Docker
container.  (Some CWL platforms, such as Arvados, require it).  Many
bioinformatics tools are already available as containers.  One
resource is the BioContainers project.

Visit https://biocontainers.pro/

Click on "Registry"

Search for "subread"

Click on the search result for "subread"

Click on the tab "Packages and Containers"

Choose a row with type "docker", then click the "copy to clipboard"
button on the right side of the"Full Tag" column for that row.
