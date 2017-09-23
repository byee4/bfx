#!/bin/bash

python /home/bay001/projects/codebase/bfx/pyscripts/single_cell/filter_mtx.py \
--mtx /projects/ps-yeolab3/jobs/cellranger/yan_takeda_tph1_colon_gfp/results/temp/YS_1/outs/filtered_gene_bc_matrices/mm10gfp/matrix.mtx \
--genes /projects/ps-yeolab3/jobs/cellranger/yan_takeda_tph1_colon_gfp/results/temp/YS_1/outs/filtered_gene_bc_matrices/mm10gfp/genes.tsv \
--barcodes /projects/ps-yeolab3/jobs/cellranger/yan_takeda_tph1_colon_gfp/results/temp/YS_1/outs/filtered_gene_bc_matrices/mm10gfp/barcodes.tsv \
--filtered_mtx /home/bay001/projects/takeda_singlecell_20170802/data/yan_takeda_tph1_colon_gfp/filtered_gene_bc_matrices/mm10gfp_test/matrix.mtx \
--filtered_genes /home/bay001/projects/takeda_singlecell_20170802/data/yan_takeda_tph1_colon_gfp/filtered_gene_bc_matrices/mm10gfp_test/genes.txt \
--filtered_barcodes /home/bay001/projects/takeda_singlecell_20170802/data/yan_takeda_tph1_colon_gfp/filtered_gene_bc_matrices/mm10gfp_test/barcodes.txt