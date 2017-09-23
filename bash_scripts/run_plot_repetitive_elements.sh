#!/bin/bash

python /home/bay001/projects/codebase/bfx/pyscripts/clipseq/plot_repetitive_elements.py \
--sam /home/bay001/projects/mpileup_maps_20161108/temporary_data/ALLCLIP_repmapping_20170201/204_01_RBFOX2.combined_w_uniquemap.rmDup.sam \
--genelist /home/bay001/projects/codebase/bfx/pyscripts/data/genelist.txt \
--outdir /home/bay001/projects/codebase/temp/RNU/ \
--reference /home/bay001/projects/codebase/bfx/pyscripts/data/MASTER_filelist.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.fa