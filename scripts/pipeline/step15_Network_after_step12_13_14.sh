#!/bin/bash
set -e

CONFIG=$1
GAMMA_HOME=$(eval echo $(yq .input.GAMMA_HOME "${CONFIG}"))
trait_name=$(eval echo $(yq .input.trait "${CONFIG}"))
SCRIPT_DIR=$(eval echo $(yq .script.path "${CONFIG}"))

bash ${SCRIPT_DIR}/Network/Network.sh ${CONFIG}
