#!/bin/bash

for FQ in $( ls inputs/R2_HCMM*.fqTr.fq.gz ); do
  STAR \
  --alignEndsType EndToEnd \
  --genomeDir /projects/ps-yeolab3/bay001/annotations/hg19/star_2_7_6a_gencode19_sjdb/ \
  --genomeLoad NoSharedMemory \
  --outBAMcompression 10 \
  --outFileNamePrefix ${FQ%.*} \
  --outFilterMultimapNmax 1 \
  --outFilterMultimapScoreRange 1 \
  --outFilterScoreMin 10 \
  --outFilterType BySJout \
  --outReadsUnmapped Fastx \
  --outSAMattrRGline ID:foo \
  --outSAMattributes All \
  --outSAMmode Full \
  --outSAMtype BAM Unsorted \
  --outSAMunmapped Within \
  --outStd Log \
  --readFilesIn ${FQ} \
  --runMode alignReads \
  --readFilesCommand zcat \
  --runThreadN 8;
done