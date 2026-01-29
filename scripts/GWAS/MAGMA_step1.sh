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

mkdir -p ${OUTPUT}/MAGMA/detail
mkdir -p ${OUTPUT}/MAGMA/summary

# ------------------------------------------------------------------------
#  MAGMA analysis
# ------------------------------------------------------------------------
MAGMA=$(eval echo $(yq .software.magma "${CONFIG}"))
REFERENCE=$(eval echo $(yq .reference.reference_bfile_GRCh37 "${CONFIG}"))
gene_annot=$(eval echo $(yq .magma.gene_annot "${CONFIG}"))

# ----
chr=$(eval echo $(yq .input.chr "${CONFIG}"))
i=${chr}
# i=${SLURM_ARRAY_TASK_ID}
${MAGMA} \
	--bfile ${REFERENCE}_chr${i} \
	--gene-annot ${gene_annot} \
	--pval ${GWAS_DATA} ncol=N \
	--gene-model snp-wise=mean \
	--out ${OUTPUT}/MAGMA/detail/${trait_name}_chr${i}


