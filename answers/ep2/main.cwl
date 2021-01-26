### 1. File header
cwlVersion: v1.2
class: Workflow
label: RNAseq CWL practice workflow

### 2. Workflow Inputs
inputs:
  fq: File
  genome: Directory
  gtf: File

### 3. Workflow Steps
steps:
  fastqc:
    run: bio-cwl-tools/fastqc/fastqc_2.cwl
    in:
      reads_file: fq
    out: [html_file]

  ### 4. Running alignment with STAR
  STAR:
    requirements:
      ResourceRequirement:
        ramMin: 6000
    run: bio-cwl-tools/STAR/STAR-Align.cwl
    in:
      RunThreadN: {default: 4}
      GenomeDir: genome
      ForwardReads: fq
      OutSAMtype: {default: BAM}
      OutSAMunmapped: {default: Within}
    out: [alignment]

  ### 5. Running samtools
  samtools:
    run: bio-cwl-tools/samtools/samtools_index.cwl
    in:
      bam_sorted: STAR/alignment
    out: [bam_sorted_indexed]

### 7. Workflow Outputs
outputs:
  qc_html:
    type: File
    outputSource: fastqc/html_file
  bam_sorted_indexed:
    type: File
    outputSource: samtools/bam_sorted_indexed
