# Expressions

1. Expressions on step inputs

You might have noticed that the output bam files are all named
`Aligned.sortedByCoord.out.bam`.  This happens because because when we
call STAR, it gives the output a default file name.

Now, during workflow execution, this is usually not a problem.  The
workflow runner is smart enough to know that these files are different
and keep them separate.  This can even make development easier by not
having to worry about assigning unique file names to every file.

However, it is a problem for humans interpreting the output.  We can
fix this by setting the parameter `OutFileNamePrefix` on STAR.

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

2. Organizing output files into Directories

This is a more advanced example.

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

3. Other uses for expressions

- Creating configuration files on the fly.
