#!/bin/bash
set -e

# ------------------------------------------------------------------------
#  Input
# ------------------------------------------------------------------------
CONFIG=$1
SCRIPT_DIR=`yq .script.path "${CONFIG}"`
GWAS_DATA=`yq .input.gwas "${CONFIG}"`
trait_name=`yq .input.trait "${CONFIG}"`
OUTPUT=`yq .input.output "${CONFIG}"`


mkdir -p ${OUTPUT}/COLOC/detail
mkdir -p ${OUTPUT}/COLOC/summary
mkdir -p ${OUTPUT}/COLOC/tmp


# ------------------------------------------------------------------------
#  COLOC analysis
# ------------------------------------------------------------------------
reference_freq=`yq .reference.reference_freq "${CONFIG}"`
QTL_list=`yq .smr.QTL_list "${CONFIG}"`

# qtl_i=${SLURM_ARRAY_TASK_ID}
# qtl_num=`cat $QTL_list | wc -l`
# only consider eQTL, sQTL, pQTL
qtl_num=4
for qtl_i in $(seq 1 $qtl_num); do

qtl_name=`head -n ${qtl_i} ${QTL_list} | tail -n1 | awk -F "\t" '{print $1}'`
qtl_data=`head -n ${qtl_i} ${QTL_list} | tail -n1 | awk -F "\t" '{print $2}'`
qtl_chr=`head -n ${qtl_i} ${QTL_list} | tail -n1 | awk -F "\t" '{print $3}'`
qtl_n=`head -n ${qtl_i} ${QTL_list} | tail -n1 | awk -F "\t" '{print $4}'`

echo "Processing QTL $qtl_i / $qtl_num ..."
echo "QTL name: $qtl_name"

# COLOC=`yq .software.coloc "${CONFIG}"`
# COLOC="${SCRIPT_DIR}/L2G/COLOC.R"
COLOC="${SCRIPT_DIR}/L2G/COLOC_within_locus.R"

R_functions=`yq .software.R_functions "${CONFIG}"`
QTL_dir=`yq .coloc.QTL_dir "${CONFIG}"`
env=`yq .environment.R_421 "${CONFIG}"`
source activate ${env}
# ----
# for i in $(seq 1 22); do
    chr=`yq .input.chr "${CONFIG}"`
    i=${chr}
    if [ "$qtl_chr" = "TRUE" ]; then
        QTL_data="${qtl_data}${i}"
    else
        QTL_data="${qtl_data}"
    fi

SECONDS=0

if [ -f "${QTL_dir}/QTL_data/${qtl_name}_chr${i}.txt" ]; then
echo "Running COLOC for chromosome $i ..."
Rscript ${COLOC} \
	--gwas ${GWAS_DATA} \
	--reference_freq ${reference_freq} \
	--qtl_query  ${QTL_dir}/QTL_data/${qtl_name}_chr${i}.txt \
	--qtl_number ${qtl_n} \
	--R_functions ${R_functions} \
	--locus_file ${OUTPUT}/Clumping/summary/${trait_name}.locus \
	--out ${OUTPUT}/COLOC/detail/${trait_name}_${qtl_name}_chr${i}.coloc
fi

elapsed_time=$SECONDS
result_file="${OUTPUT}/COLOC/detail/${trait_name}_${qtl_name}_chr${i}.coloc"
current_time=$(date "+%Y-%m-%d %H:%M:%S")

echo "COLOC analysis results have been written into [$result_file]"
echo "COLOC analysis completed: $current_time"
echo "COLOC analysis computational time: $(($elapsed_time / 3600)):$((($elapsed_time / 60) % 60)):$(($elapsed_time % 60))"

# done

awk 'NR==1 || FNR>1' ${OUTPUT}/COLOC/detail/${trait_name}_${qtl_name}_chr*.coloc > ${OUTPUT}/COLOC/summary/${trait_name}_${qtl_name}_chrALL.coloc 

done

