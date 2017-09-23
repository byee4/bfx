#!/bin/bash

### runs getfasta on introns.bed files from UCSC table browser to get intron sequences.
### intermediate step to be able to find ATAC introns within the gencode v19 reference.

bedtools getfasta \
-fi /projects/ps-yeolab/genomes/hg19/chromosomes/all.fa \
-bed /projects/ps-yeolab3/bay001/annotations/gencode.v19.annotations.introns.bed \
 > /home/bay001/projects/encode/analysis/atac_intron_analysis/atac_introns_from_ucsc/all_introns.fasta