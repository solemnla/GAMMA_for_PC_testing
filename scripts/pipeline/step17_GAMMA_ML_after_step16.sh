#!/bin/bash
set -e

CONFIG=$1
GAMMA_HOME=$(eval echo $(yq .input.GAMMA_HOME "${CONFIG}"))
SCRIPT_DIR=$(eval echo $(yq .script.path "${CONFIG}"))

bash ${SCRIPT_DIR}/GAMMA_ML/GAMMA_ML.sh ${CONFIG}
