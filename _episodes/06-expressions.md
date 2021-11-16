---
title: "Dynamic Workflow Behavior"
teaching: 20
exercises: 10
questions:
- "What kind of custom logic can happen between steps?"
objectives:
- "Customize the STAR output filename to use the input filename."
- "Organize files into directories."
keypoints:
- "CWL expressions allow you to use custom logic to determine input parameter values."
- "CWL ExpressionTool can be used to reshape data, such as declaring directories that should contain output files."
---

# Expressions on step inputs

You might have noticed that the output bam files are all named
`Aligned.sortedByCoord.out.bam`.  This happens because because when we
call STAR, it gives the output a default file name.

During workflow execution, this is usually not a problem.  The
workflow runner is smart enough to know that these files are different
and keep them separate.  This can even make development easier by not
having to worry about assigning unique file names to every file.
Also, if we intend to discard the BAM files as intermediate results

However, it is a problem for humans interpreting the output.  We can
fix this by setting the parameter `OutFileNamePrefix` on STAR.  We
want the output filename to be based on the input filename.

In `alignment.cwl`, we can use `valueFrom` on the `OutFileNamePrefix`
input parameter to construct the output prefix from the input
filename.

```
requirements:
  StepInputExpressionRequirement: {}

steps:
  ...
  STAR:
    ...
    in:
      ForwardReads: fq
      ...
      OutFileNamePrefix: {valueFrom: "$(inputs.ForwardReads.nameroot)."}
```
{: .language-yaml }

The code between `$(...)` is called an "expression".  It is evaluated
when setting up the step to run, and the expression is replaced by the
result to get the parameter value.

An expression can refer to other inputs to the step that are either
directly connected to another value, or have a default value.  Here,
we refer to the input parameter ForwardReads, which is our fastq input
file.

ForwardReads is a File object, not a plain file name, so it has some
fields of its own.  The file object has a number of fields that we can
use.  These include `basename` (the name of the file, without a
directory), `size` (file size, in bytes), `nameext` (the last file
extension) and `nameroot` (the name with `nameext` removed).  Using

Finally, our expression is embedded in a string, so after replacing
the expression with the value of `inputs.ForwardReads.nameroot`, it
adds the remainder of the string, which just is a dot `.`.  This is to
separate the leading part of our filename from the "Aligned.bam"
extension that will be added by STAR.

# Organizing output files into Directories

You probably noticed that all the output files appear in the same
directory.  You might prefer that each file appears in its own
directory.  This will show you how to do that.

Unlike shell scripts, in CWL you cannot call `mkdir` and `mv` to
organize files into directories.  This is because the output files, at
this point, do not actually exist together in one directory.  They may
exist on different nodes, in different directories, or different cloud
buckets.  In CWL, instead of moving files around directly, you tell
the runner you want your directory to look like, and it will create it
for you.

We can use an "expression" to create a `Directory` object describing
each of our directories.  An expression is a piece of Javascript code
which is executed by the workflow runner to produce values that affect
the workflow.  These can be a simple as substituting the value of an
input variable, or as complex as en entire function that generates new
objects.

Javscript code must be bracketed inside `$(...)` or `${...}`. The
difference comes down to syntax.  The `$()` form is more compact but
can only include code that can appear on the right side of an
assignment (`=`), which cannot include control blocks like `if` or
`for`.  The `${}` form is a Javascript function, which can include
control blocks, and must end in a `return` statement.

Expressions can both appear in `valueFrom` fields as well as some
other fields, or in an `ExpressionTool` which, like `Workflow` or
`CommandLineTool` has explicitly defined `inputs` and `outputs`
sections.

The approach here is to define an expression tool which takes a

The `Directory` object has two fields, `basename` and `listing`.  The
`basename` is the name of the directory, and the `listing` is the
contents, which consists of other File and Directory objects.

Create `subdirs.cwl`:

```
cwlVersion: v1.2
class: ExpressionTool
requirements:
  InlineJavascriptRequirement: {}
inputs:
  fq: File[]
  bams: File[]
  qc: File[]
outputs:
  dirs: Directory[]
expression: |-
  ${
  var dirs = [];
  for (var i = 0; i < inputs.bams.length; i++) {
    dirs.push({
      "class": "Directory",
      "basename": inputs.fq[i].nameroot,
      "listing": [inputs.bams[i], inputs.qc[i]]
    });
  }
  return {"dirs": dirs};
  }
```
{: .language-yaml }

Then change `main.cwl`:

```
steps:
  ...
  output-subdirs:
    run: subdirs.cwl
    in:
      fq: fq
      bams: alignment/bam_sorted_indexed
      qc: alignment/qc_html
    out: [dirs]
outputs:
  dirs:
    type: Directory[]
    outputSource: output-subdirs/dirs
  featurecounts:
    type: File
    outputSource: featureCounts/featurecounts
```
{: .language-yaml }

> ## Running the workflow
>
> Run the workflow.  Look at the output.  The BAM and fastqc files
> should now be organized into directories, with better naming of the
> bam files.
>
{: .challenge }

> ## Episode solution
> * <a href="../assets/answers/ep6/main.cwl">main.cwl</a>
> * <a href="../assets/answers/ep6/alignment.cwl">alignment.cwl</a>
> * <a href="../assets/answers/ep6/featureCounts.cwl">featureCounts.cwl</a>
> * <a href="../assets/answers/ep6/subdirs.cwl">subdirs.cwl</a>
{: .solution}
