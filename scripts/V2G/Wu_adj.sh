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


mkdir -p ${OUTPUT}/Wu_adj/LD
mkdir -p ${OUTPUT}/Wu_adj/Result
mkdir -p ${OUTPUT}/Wu_adj/summary


# ------------------------------------------------------------------------
# Wu_adj analysis
# ------------------------------------------------------------------------
reference_bfile=$(eval echo $(yq .reference.reference_bfile "${CONFIG}"))
reference_all_bim=$(eval echo $(yq .reference.reference_all_bim "${CONFIG}"))
reference_freq=$(eval echo $(yq .reference.reference_freq "${CONFIG}"))

plink1_9=$(eval echo $(yq .software.plink1_9 "${CONFIG}"))
plink2=$(eval echo $(yq .software.plink2 "${CONFIG}"))

env=$(eval echo $(yq .environment.R_421 "${CONFIG}"))
source activate $env

# --------------------------
# No matter Clumping or COJO, we only take the CHR/POS/SNP information.
# if input is COJO:

# COJO_or_Clumping_file="${OUTPUT}/Clumping/summary/${trait_name}.clumped"
# COJO_or_Clumping_file="${OUTPUT}/COJO/summary/${trait_name}.jma.cojo"
locus_path="${OUTPUT}/Clumping/summary/${trait_name}.locus"

Rscript ${SCRIPT_DIR}/V2G/Wu_adj.R \
	${GWAS_DATA} \
	${trait_name} \
	${OUTPUT} \
	${locus_path} \
	${plink1_9} \
	${reference_all_bim} \
	${reference_freq} \
	${reference_bfile}
	
