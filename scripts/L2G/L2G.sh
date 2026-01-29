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

mkdir -p ${OUTPUT}/L2G/summary
mkdir -p ${OUTPUT}/L2G/score
mkdir -p ${OUTPUT}/L2G/feature
mkdir -p ${OUTPUT}/L2G/plot

# ------------------------------------------------------------------------
#  L2G summary analysis
# ------------------------------------------------------------------------
L2G_column_file=$(eval echo $(yq .l2g.L2G_column_file "${CONFIG}"))
QTL_name_mapping_file=$(eval echo $(yq .l2g.QTL_name_mapping_file "${CONFIG}"))
gsmr_pQTL_file=$(eval echo $(yq .l2g.gsmr_pQTL_file "${CONFIG}"))
pQTL_epi_file=$(eval echo $(yq .l2g.pQTL_epi_file "${CONFIG}"))

gamma_gene_file=$(eval echo $(yq .gene.gamma_gene "${CONFIG}"))
R_functions_file=$(eval echo $(yq .software.R_functions "${CONFIG}"))

env=$(eval echo $(yq .environment.R_421 "${CONFIG}"))
source activate $env

Rscript ${SCRIPT_DIR}/L2G/L2G_summary.R \
	${trait_name} \
	${OUTPUT} \
	${QTL_name_mapping_file} \
	${gsmr_pQTL_file} \
	${pQTL_epi_file} \
	${gamma_gene_file} \
	${R_functions_file} \
	${L2G_column_file}

