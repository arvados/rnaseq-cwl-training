#!/bin/bash

# Generate STAR genome index

STAR --runThreadN 4 \
--runMode genomeGenerate \
--genomeDir unix_lesson/reference_data \
--genomeFastaFiles unix_lesson/reference_data/chr1.fa \
--sjdbGTFfile unix_lesson/reference_data/chr1-hg19_genes.gtf \
--sjdbOverhang 99
