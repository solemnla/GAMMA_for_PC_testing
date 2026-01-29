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
bash scripts/pipeline/step06_SMR.sh ${CONFIG}
bash scripts/pipeline/step07_xMAGIC_after_step06.sh ${CONFIG}
bash scripts/pipeline/step08_COLOC.sh ${CONFIG}
bash scripts/pipeline/step09_FUSION.sh ${CONFIG}
bash scripts/pipeline/step10_GSMR_Rcode.sh ${CONFIG}
bash scripts/pipeline/step11_L2G_after_step06_07_08_09_10.sh ${CONFIG}
bash scripts/pipeline/step12_PoPS_after_step02.sh ${CONFIG}
bash scripts/pipeline/step13_DEPICT.sh ${CONFIG}
bash scripts/pipeline/step14_RWR_PPR.sh ${CONFIG}
bash scripts/pipeline/step15_Network_after_step12_13_14.sh ${CONFIG}
bash scripts/pipeline/step16_GAMMA_after_step05_11_15.sh ${CONFIG}
bash scripts/pipeline/step17_GAMMA_ML_after_step16.sh ${CONFIG}
