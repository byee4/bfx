#!/bin/bash

for BAM in $( ls inputs/*.bam ); do
  makebigwigfiles \
  --bam ${BAM} \
  --genome /projects/ps-yeolab/genomes/mm9/mm9.chrom.sizes \
  --bw_pos ${BAM%.*}.norm.pos.bw \
  --bw_neg ${BAM%.*}.norm.neg.bw \
  --direction r
done
