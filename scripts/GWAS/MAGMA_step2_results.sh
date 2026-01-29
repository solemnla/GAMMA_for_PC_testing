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
#  MAGMA results
# ------------------------------------------------------------------------
# head -n 1 ${OUTPUT}/MAGMA/detail/${trait_name}_chr1.genes.out > ${OUTPUT}/MAGMA/summary/${trait_name}_chrALL.genes.out
# for chr in {1..22}
# do
# 	tail -n +2 ${OUTPUT}/MAGMA/detail/${trait_name}_chr${chr}.genes.out >> ${OUTPUT}/MAGMA/summary/${trait_name}_chrALL.genes.out 
# done

# head -n 2 ${OUTPUT}/MAGMA/detail/${trait_name}_chr1.genes.raw > ${OUTPUT}/MAGMA/summary/${trait_name}_chrALL.genes.raw
# for chr in {1..22}
# do
# 	tail -n +3 ${OUTPUT}/MAGMA/detail/${trait_name}_chr${chr}.genes.raw >> ${OUTPUT}/MAGMA/summary/${trait_name}_chrALL.genes.raw 
# done


chr=$(eval echo $(yq .input.chr "${CONFIG}"))
head -n 1 ${OUTPUT}/MAGMA/detail/${trait_name}_chr${chr}.genes.out > ${OUTPUT}/MAGMA/summary/${trait_name}_chrALL.genes.out
# for chr in {1..22}
# do
tail -n +2 ${OUTPUT}/MAGMA/detail/${trait_name}_chr${chr}.genes.out >> ${OUTPUT}/MAGMA/summary/${trait_name}_chrALL.genes.out 
# done

head -n 2 ${OUTPUT}/MAGMA/detail/${trait_name}_chr${chr}.genes.raw > ${OUTPUT}/MAGMA/summary/${trait_name}_chrALL.genes.raw
# for chr in {1..22}
# do
tail -n +3 ${OUTPUT}/MAGMA/detail/${trait_name}_chr${chr}.genes.raw >> ${OUTPUT}/MAGMA/summary/${trait_name}_chrALL.genes.raw 
# done

