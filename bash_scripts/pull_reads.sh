#!/usr/bin/env bash

# BAM=/home/bay001/projects/codebase/junction-counter/data/ENCFF468DJE.bam
# READNAME='DF8F08P1:365:C60G9ACXX:5:2111:10416:51891'
# OUTPUT_FILE=/home/bay001/projects/codebase/junction-counter/junction_counter/test/bams/inclusion_9.bam

samtools view -H $1 > tmp.sam

samtools view $1 | grep -m 2 $2 >> tmp.sam;
cat tmp.sam | samtools view -bS > $3;
samtools index $3