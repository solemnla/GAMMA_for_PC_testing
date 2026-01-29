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

mkdir -p ${OUTPUT}/MR/GSMR/detail
mkdir -p ${OUTPUT}/MR/GSMR/summary
mkdir -p ${OUTPUT}/MR/GSMR/results

# ------------------------------------------------------------------------
# Other MR analysis
# ------------------------------------------------------------------------
MR_snp_min=$(eval echo $(yq .mr.MR_snp_min "${CONFIG}"))
default_MR_snp_min=3
MR_snp_min=${MR_snp_min:-$default_MR_snp_min}

default_batch_num=10
batch_num=$(eval echo $(yq .mr.batch_num "${CONFIG}"))
batch_num=${batch_num:-$default_batch_num}

QTL_clumped_list=$(eval echo $(yq .mr.GSMR_QTL_clumped_list "${CONFIG}"))

# -------------------------------
# QTL clumped results
# qtl_i=${SLURM_ARRAY_TASK_ID}

# For demo, only run 1 dataset.
qtl_i=6

# qtl_name=`awk -F "\t" -v row=$qtl_i 'NR==row {print $1}' $QTL_clumped_list`
# qtl_clumpled_file=`awk -F "\t" -v row=$qtl_i 'NR==row {print $2}' $QTL_clumped_list`
qtl_name=$(awk -F "\t" -v row=$qtl_i 'NR==row {print $1}' $QTL_clumped_list)
qtl_clumpled_file=$(eval echo $(awk -F "\t" -v row=$qtl_i 'NR==row {print $2}' $QTL_clumped_list))

if [ -z "$qtl_name" ]; then
    echo "no QTL in line $qtl_i, skip"
    exit 0
fi

env=$(eval echo $(yq .environment.R_421 "${CONFIG}"))
source activate $env


pids=()
for ((process_num=1; process_num<=${batch_num}; process_num=process_num+1)); do

Rscript ${SCRIPT_DIR}/L2G/MR/GSMR_MR.R \
	${SCRIPT_DIR} \
	${trait_name} \
	${GWAS_DATA} \
	${qtl_name} \
	${qtl_clumpled_file} \
	${MR_snp_min} \
	${OUTPUT} \
	${process_num} &

	pids+=($!)
	if (( ${process_num} % ${batch_num} == 0 )); then
		for pid in "${pids[@]}"; do
			wait $pid || [ $? -eq 99 ] && true || exit $?
		done
		pids=()
	fi
done

for pid in "${pids[@]}"; do
    wait $pid || [ $? -eq 99 ] && true || exit $?
done

set +e
awk 'NR==1 || FNR>1' ${OUTPUT}/MR/GSMR/detail/${trait_name}_${qtl_name}_*_MR_comparison.txt > ${OUTPUT}/MR/GSMR/summary/${trait_name}_${qtl_name}_MR_comparison.txt
# rm ${OUTPUT}/MR/detail/${trait_name}_${qtl_name}_*__MR_comparison.txt
awk 'NR==1 || FNR>1' ${OUTPUT}/MR/GSMR/summary/${trait_name}_*_MR_comparison.txt > ${OUTPUT}/MR/GSMR/results/${trait_name}.gsmr
set -e


