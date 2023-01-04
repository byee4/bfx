#!/bin/bash

for BAM in $( ls $@ ); do
  echo ${BAM};
  samtools view -cF 4 ${BAM}
done
