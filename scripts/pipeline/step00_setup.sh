#!/bin/bash
set -e

CONFIG=$1
GAMMA_HOME=$(eval echo $(yq .input.GAMMA_HOME "${CONFIG}"))
OUTPUT=$(eval echo $(yq .input.output "${CONFIG}"))

mkdir -p ${OUTPUT}
echo "Result directory created at ${OUTPUT}!"

INPUT=$(eval echo $(yq .input.gwas_raw "${CONFIG}"))
echo "Input GWAS: ${INPUT}"

