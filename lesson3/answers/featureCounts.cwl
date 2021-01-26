### 1. File header
cwlVersion: v1.2
class: CommandLineTool

### 2. Command line tool inputs
inputs:
  gtf: File
  counts_input_bam: File

### 3. Specifying the program to run
baseCommand: featureCounts

### 4. Command arguments
arguments: [-T, $(runtime.cores),
            -a, $(inputs.gtf),
            -o, featurecounts.tsv,
            $(inputs.counts_input_bam)]

### 5. Outputs section
outputs:
  featurecounts:
    type: File
      outputBinding:
      glob: featurecounts.tsv

### 6. Running in a container
hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/subread:1.5.0p3--0
