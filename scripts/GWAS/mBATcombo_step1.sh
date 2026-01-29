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


mkdir -p ${OUTPUT}/mBATcombo/detail
mkdir -p ${OUTPUT}/mBATcombo/summary

# ------------------------------------------------------------------------
#  mBATcombo analysis
# ------------------------------------------------------------------------
GCTA=$(eval echo $(yq .software.gcta "${CONFIG}"))
REFERENCE=$(eval echo $(yq .reference.reference_bfile "${CONFIG}"))
gene_list=$(eval echo $(yq .mBAT.gene_list "${CONFIG}"))

# ----
chr=$(eval echo $(yq .input.chr "${CONFIG}"))
i=${chr}
# i=${SLURM_ARRAY_TASK_ID}
${GCTA} --bfile ${REFERENCE}_chr${i} \
	--mBAT-combo ${GWAS_DATA} \
	--mBAT-gene-list ${gene_list} \
	--mBAT-print-all-p \
	--chr ${i} \
	--out ${OUTPUT}/mBATcombo/detail/${trait_name}_mBATcombo_chr${i}

