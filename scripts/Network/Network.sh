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

mkdir -p ${OUTPUT}/Network/summary
mkdir -p ${OUTPUT}/Network/score
mkdir -p ${OUTPUT}/Network/feature
mkdir -p ${OUTPUT}/Network/plot


# ------------------------------------------------------------------------
#  Network analysis
# ------------------------------------------------------------------------
R_functions=$(eval echo $(yq .software.R_functions "${CONFIG}"))
gamma_gene=$(eval echo $(yq .gene.gamma_gene "${CONFIG}"))

env=$(eval echo $(yq .environment.R_421 "${CONFIG}"))
source activate $env
Rscript  ${SCRIPT_DIR}/Network/Network_summary.R \
	${trait_name} \
	${OUTPUT} \
	${gamma_gene} \
	${R_functions}
