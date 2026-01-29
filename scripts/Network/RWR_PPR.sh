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

mkdir -p ${OUTPUT}/RWR_PPR/summary

# ------------------------------------------------------------------------
#  RWR and PPR analysis
# ------------------------------------------------------------------------
env=$(eval echo $(yq .environment.python_rwr_ppr "${CONFIG}"))
source activate $env
ppi_file=$(eval echo $(yq .rwr_ppr.ppi_file "${CONFIG}"))


python ${SCRIPT_DIR}/Network/RWR_PPR.py \
	${trait_name} \
	${OUTPUT} \
	${ppi_file}
