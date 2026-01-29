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

mkdir -p ${OUTPUT}/GAMMA/score
mkdir -p ${OUTPUT}/GAMMA/feature
mkdir -p ${OUTPUT}/GAMMA/plot
mkdir -p ${OUTPUT}/GWAS/manhattan_plot

R_functions_file=$(eval echo $(yq .software.R_functions "${CONFIG}"))
gamma_gene_file=$(eval echo $(yq .gene.gamma_gene "${CONFIG}"))
Pharmaprojects_data_file=$(eval echo $(yq .gamma.pharmaprojects "${CONFIG}"))
reference_all_bim_file=$(eval echo $(yq .reference.reference_all_bim "${CONFIG}"))

MeSH_id=$(eval echo $(yq .input.mesh_id "${CONFIG}"))

# T2D
echo ${trait_name}
echo ${OUTPUT}
echo ${Pharmaprojects_data_file}
echo ${gamma_gene_file}
echo ${R_functions_file}

env=$(eval echo $(yq .environment.R_421 "${CONFIG}"))
source activate ${env}
echo ${MeSH_id}

echo "------------"
Rscript ${SCRIPT_DIR}/GAMMA/GAMMA_summary.R \
	${trait_name} \
	${OUTPUT} \
	${Pharmaprojects_data_file} \
	${gamma_gene_file} \
	${R_functions_file} \
	${MeSH_id}

Rscript ${SCRIPT_DIR}/GAMMA/manhattan_plot.R \
	${trait_name} \
	${GWAS_DATA} \
	${OUTPUT} \
	${reference_all_bim_file}


