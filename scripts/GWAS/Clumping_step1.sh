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
clump_field=$(eval echo $(yq .clumping.field "${CONFIG}"))

if [ "$clump_field" = "null" ]; then
  clump_field="p"
fi
# clump_field=`awk 'NR==1 {print $7}' ${GWAS_DATA}`

mkdir -p ${OUTPUT}/Clumping/detail
mkdir -p ${OUTPUT}/Clumping/summary

# ------------------------------------------------------------------------
#  Clumping analysis
# ------------------------------------------------------------------------

 
plink1_9=$(eval echo $(yq .software.plink1_9 "${CONFIG}"))
REFERENCE=$(eval echo $(yq .reference.reference_bfile "${CONFIG}"))

# ----
chr=$(eval echo $(yq .input.chr "${CONFIG}"))
i=${chr}
# i=${SLURM_ARRAY_TASK_ID}
${plink1_9} \
    --bfile ${REFERENCE}_chr${i} \
    --chr ${i} \
    --maf 0.01 \
    --clump ${GWAS_DATA} \
    --clump-p1 5e-8 \
    --clump-p2 5e-8 \
    --clump-r2 0.05 \
    --clump-kb 1000 \
    --clump-snp-field SNP \
    --clump-field ${clump_field} \
    --out ${OUTPUT}/Clumping/detail/${trait_name}_chr${i}


