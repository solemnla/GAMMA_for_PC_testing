#!/bin/bash
set -e

# ------------------------------------------------------------------------
#  Input
# ------------------------------------------------------------------------
CONFIG=$1
GAMMA_HOME=$(eval echo $(yq .input.GAMMA_HOME "${CONFIG}"))
SCRIPT_DIR=$(eval echo $(yq .script.path "${CONFIG}"))
GWAS_DATA=$(eval echo $(yq .input.gwas "${CONFIG}"))
trait_name=$(eval echo $(yq .input.trait "${CONFIG}"))
OUTPUT=$(eval echo $(yq .input.output "${CONFIG}"))

# ------------------------------------------------------------------------
#  Clumping results
# ------------------------------------------------------------------------
gencode=$(eval echo $(yq .gene.gencode "${CONFIG}"))

env=$(eval echo $(yq .environment.R_421 "${CONFIG}"))
source activate $env

# ----
awk 'NR==1 || FNR>1' ${OUTPUT}/Clumping/detail/${trait_name}_chr*.clumped > ${OUTPUT}/Clumping/summary/${trait_name}.clumped 

if [ -f "${OUTPUT}/Clumping/summary/${trait_name}.clumped" ]; then

    Rscript ${SCRIPT_DIR}/GWAS/Clumping_step2_results.R ${trait_name} ${OUTPUT}

    Rscript ${SCRIPT_DIR}/GWAS/Clumping_step2_results_locus_window.R ${trait_name} ${OUTPUT} 1000000 ${gencode}
    Rscript ${SCRIPT_DIR}/GWAS/Clumping_step2_results_locus_window.R ${trait_name} ${OUTPUT} 750000 ${gencode}
    Rscript ${SCRIPT_DIR}/GWAS/Clumping_step2_results_locus_window.R ${trait_name} ${OUTPUT} 500000 ${gencode}
    Rscript ${SCRIPT_DIR}/GWAS/Clumping_step2_results_locus_window.R ${trait_name} ${OUTPUT} 250000 ${gencode}

fi
