#!/bin/bash
set -e

# Set up config files
CONFIG=/home/sunshufeng/GAMMA_for_PC_testing/deploy/GAMMA.yaml
GAMMA_HOME=$(eval echo $(yq .input.GAMMA_HOME "${CONFIG}"))
cd ${GAMMA_HOME}

# Run pipeline
bash scripts/pipeline/step00_setup.sh ${CONFIG}
bash scripts/pipeline/step01_Clumping.sh ${CONFIG} 
bash scripts/pipeline/step02_MAGMA.sh ${CONFIG}
bash scripts/pipeline/step03_mBATcombo.sh ${CONFIG}
bash scripts/pipeline/step04_Wu_adj_after_step01.sh ${CONFIG}
bash scripts/pipeline/step05_V2G_after_step04.sh ${CONFIG}
