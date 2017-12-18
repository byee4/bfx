#!/bin/bash

python /home/bay001/projects/codebase/rbp-maps/maps/plot_peak.py \
-i /home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_FINALforpapers_20170325/321_02.basedon_321_02.peaks.l2inputnormnew.bed.compressed.bed \
-o /projects/ps-yeolab3/bay001/gabe_qc_20170612/permanent_data/integrated_paper_maps/321_02.basedon_321_02.peaks.l2inputnormnew.bed.compressed.png \
-m /projects/ps-yeolab3/bay001/maps/current_annotations/as_miso_renamed/HNRNPUL1-BGHLV26-HepG2-included-upon-knockdown \
/projects/ps-yeolab3/bay001/maps/current_annotations/as_miso_renamed/HNRNPUL1-BGHLV26-HepG2-excluded-upon-knockdown \
/projects/ps-yeolab3/bay001/maps/current_annotations/as_miso_renamed/HepG2-constitutive-exons.miso \
/projects/ps-yeolab3/bay001/maps/current_annotations/as_miso_renamed/HepG2-native-cassette-exons.miso \
/projects/ps-yeolab3/bay001/maps/current_annotations/as_miso_renamed/HepG2-native-included-exons.miso \
/projects/ps-yeolab3/bay001/maps/current_annotations/as_miso_renamed/HepG2-native-excluded-exons.miso \
-bgnum 4 \
-p 3 \
-f 0