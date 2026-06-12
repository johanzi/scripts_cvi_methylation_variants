#!/bin/bash
#===============================================================================#
# FILE: ld_Chr5_15047549.sh
# DESCRIPTION: ---
# AUTHOR: Emmanuel Tergemina
# ORGANIZATION: Department of Plant Developmental Biology, Max Planck Institute for Plant Breeding Research
# VERSION: 1.0
# CREATED: 20221222
#===============================================================================#

path2plink=/srv/netscratch/irg/grp_hancock/Manu/software/plink
wd=/srv/netscratch/irg/grp_hancock/johan/GWAS/dna_methylation/GWAS_83_SA/20221222/

${path2plink} --file ${wd}subset_83_only_chr_biallelic_only_alt_DP3_GQ25_wo_singletons \
--r2 \
--ld-snp 5:15047549 \
--ld-window-kb 1000 \
--ld-window 99999 \
--ld-window-r2 0 \
--out ${wd}ld/Chr5_15047549