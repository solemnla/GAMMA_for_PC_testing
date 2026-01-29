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


# ------------------------------------------------------------------------
#  mBATcombo results
# ------------------------------------------------------------------------
awk 'NR==1 || FNR>1' ${OUTPUT}/mBATcombo/detail/${trait_name}_mBATcombo_chr*.gene.assoc.mbat > ${OUTPUT}/mBATcombo/summary/${trait_name}.gene.assoc.mbat
