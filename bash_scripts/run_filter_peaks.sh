#!/bin/bash

RUNNER=/home/bay001/projects/codebase/clip_analysis/clip_analysis/filter_input_norm.py

python ${RUNNER} \
--input inputs/PARP13.WTSS1P.---.r-.fqTrTrU-SoMaSoCpSoMeV2ClN-C-Fc3Pv3.bed \
--output outputs/PARP13.WTSS1P.---.r-.fqTrTrU-SoMaSoCpSoMeV2ClN-C-Fc3Pv3.bed.minus_KOSS1P.bed \
--l10p 3 \
--l2fc 3 \
--ko_peaks inputs/PARP13.KOSS1P_caution.---.r-.fqTrTrU-SoMaSoCpSoMeV2ClN-C-Fc3Pv3.bed ;

python ${RUNNER} \
--input inputs/PARP13.WTSS2P.---.r-.fqTrTrU-SoMaSoCpSoMeV2ClN-C-Fc3Pv3.bed \
--output outputs/PARP13.WTSS2P.---.r-.fqTrTrU-SoMaSoCpSoMeV2ClN-C-Fc3Pv3.bed.minus_KOSS2P.bed \
--l10p 3 \
--l2fc 3 \
--ko_peaks inputs/PARP13.KOSS2P.---.r-.fqTrTrU-SoMaSoCpSoMeV2ClN-C-Fc3Pv3.bed ;

python ${RUNNER} \
--input inputs/PARP13.WT3P1P.---.r-.fqTrTrU-SoMaSoCpSoMeV2ClN-C-Fc3Pv3.bed \
--output outputs/PARP13.WT3P1P.---.r-.fqTrTrU-SoMaSoCpSoMeV2ClN-C-Fc3Pv3.bed.minus_KO3P1P.bed \
--l10p 3 \
--l2fc 3 \
--ko_peaks inputs/PARP13.KO3P1P.---.r-.fqTrTrU-SoMaSoCpSoMeV2ClN-C-Fc3Pv3.bed ;

### Note: This is using the input norm that I ran separately due to overwriting of the real sample ###
python ${RUNNER} \
--input inputs/PARP13.WT3P2P.---.r-.fqTrTrU-SoMaSoCpSoMeV2ClN-C-Fc0Pv0.bed \
--output outputs/PARP13.WT3P2P.---.r-.fqTrTrU-SoMaSoCpSoMeV2ClN-C-Fc3Pv3.bed.minus_KO3P2P.bed \
--l10p 3 \
--l2fc 3 \
--ko_peaks inputs/PARP13.KO3P2P_caution.---.r-.fqTrTrU-SoMaSoCpSoMeV2ClN-C-Fc3Pv3.bed ;
###

python ${RUNNER} \
--input inputs/WT3P.01v02.IDR.out.0102merged.bed \
--output outputs/WT3P.01v02.IDR.out.0102merged.bed.minus_KO3P.bed \
--l10p 3 \
--l2fc 3 \
--ko_peaks inputs/KO3P.01v02.IDR.out.0102merged.bed ;

python ${RUNNER} \
--input inputs/WTSS.01v02.IDR.out.0102merged.bed \
--output outputs/WTSS.01v02.IDR.out.0102merged.bed.minus_KOSS.bed \
--l10p 3 \
--l2fc 3 \
--ko_peaks inputs/KOSS.01v02.IDR.out.0102merged.bed ;

### Note: This is re-running RBFOX2 data for consistency checks ###

python ${RUNNER} \
--input inputs/204_01.basedon_204_01.peaks.l2inputnormnew.bed.compressed.bed \
--output outputs/204_01.basedon_204_01.peaks.l2inputnormnew.bed.compressed.Fc3Pv3.bed \
--l10p 3 \
--l2fc 3;

python ${RUNNER} \
--input inputs/204_02.basedon_204_02.peaks.l2inputnormnew.bed.compressed.bed \
--output outputs/204_02.basedon_204_02.peaks.l2inputnormnew.bed.compressed.Fc3Pv3.bed \
--l10p 3 \
--l2fc 3;

###
