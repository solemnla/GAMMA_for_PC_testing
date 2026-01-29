#!/bin/bash
set -e

# ------------------------------------------------------------------------
#  Input
# ------------------------------------------------------------------------
CONFIG=$1
GAMMA_HOME=$(eval echo $(yq .input.GAMMA_HOME "${CONFIG}"))
omim_path=$(eval echo $(yq .model_data.omim "${CONFIG}"))
clinvar_path=$(eval echo $(yq .model_data.clinvar "${CONFIG}"))
mgi_path=$(eval echo $(yq .model_data.mgi "${CONFIG}"))
gene_features=$(eval echo $(yq .model_data.gene_features "${CONFIG}"))
uniprot_path=$(eval echo $(yq .model_data.uniprot "${CONFIG}"))
trait_name=$(eval echo $(yq .input.trait "${CONFIG}"))
SCRIPT_DIR=$(eval echo $(yq .script.path "${CONFIG}"))
MESH_ID=$(eval echo $(yq .input.mesh_id "${CONFIG}"))
pharmap_path=$(eval echo $(yq .gamma.pharmaprojects_0507_shufeng "${CONFIG}"))
model_path=$(eval echo $(yq .model_model.model "${CONFIG}"))
scaler_path=$(eval echo $(yq .model_model.scaler "${CONFIG}"))
OUTPUT=$(eval echo $(yq .input.output "${CONFIG}"))

mkdir -p ${OUTPUT}/GAMMA/score

env=$(eval echo $(yq .environment.python_ml "${CONFIG}"))
source activate $env

python ${SCRIPT_DIR}/GAMMA_ML/0.0_backend.py \
    --gamma ${OUTPUT}/GAMMA/feature/${trait_name}_GAMMA.feature \
    --mesh_id ${MESH_ID} \
    --output ${OUTPUT}/GAMMA/score/ML_score.csv \
    --omim_path $omim_path \
    --clinvar_path $clinvar_path \
    --mgi_path $mgi_path \
    --gene_features $gene_features \
    --model_path $model_path \
    --scaler_path $scaler_path \
    --pharmap_path $pharmap_path

python ${SCRIPT_DIR}/GAMMA_ML/0.1_gamma_ml_result_analysis.py \
    --gamma ${OUTPUT}/GAMMA/feature/${trait_name}_GAMMA.feature \
    --mesh_id ${MESH_ID} \
    --output ${OUTPUT}/GAMMA/score/ML_score.csv \
    --omim_path $omim_path \
    --clinvar_path $clinvar_path \
    --mgi_path $mgi_path \
    --gene_features $gene_features \
    --uniprot_path $uniprot_path \
    --pharmap_path $pharmap_path
