#!/bin/bash
#PBS -N input_normalization_
#PBS -o input_normalization_.sh.out
#PBS -e input_normalization_.sh.err
#PBS -V
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=1
#PBS -A yeo-group
#PBS -q home

# Go to the directory from which the script was called
cd $PBS_O_WORKDIR
perl /home/elvannostrand/data/clip/CLIPseq_analysis/scripts/LABshare/FULL_PIPELINE_WRAPPER.pl /home/bay001/projects/codebase/data/ALLDATASETS_submittedonly.trunc3.txt /home/bay001/projects/codebase/example_outputs/ hg19

