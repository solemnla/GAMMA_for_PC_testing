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

mkdir -p ${OUTPUT}/DEPICT/GWAS_input
mkdir -p ${OUTPUT}/DEPICT/cfg
mkdir -p ${OUTPUT}/DEPICT/output


# ------------------------------------------------------------------------
# DEPICT analysis
# ------------------------------------------------------------------------
DEPICT=$(eval echo $(yq .software.depict "${CONFIG}"))
snp_loc_GRCh37=$(eval echo $(yq .depict.snp_loc_GRCh37 "${CONFIG}"))
cfg_demo=$(eval echo $(yq .depict.cfg_demo "${CONFIG}"))
env=$(eval echo $(yq .environment.python_depict "${CONFIG}"))
source activate ${env}
plink1_9=$(eval echo $(yq .software.plink1_9 "${CONFIG}"))

wk -v OFS='\t' '
    NR == FNR {
        gwas[$1] = $0; 
        next
    } 
    ($1 in gwas) {
        print gwas[$1], $2, $3
    }
' ${GWAS_DATA} ${snp_loc_GRCh37} > ${OUTPUT}/DEPICT/GWAS_input/${trait_name}_depict.txt
awk -v OFS='\t' '{if($9 != "Y" && $9 != "X" && $9 != "" && $10 != "" && $7 != "NA" && $1 != "NA") print $1, $9,$10,$7}' ${OUTPUT}/DEPICT/GWAS_input/${trait_name}_depict.txt > ${OUTPUT}/DEPICT/GWAS_input/${trait_name}_depict.clean.txt
sed -i $'1 i \tSNP\tChr\tPos\tp' ${OUTPUT}/DEPICT/GWAS_input/${trait_name}_depict.clean.txt


cp ${cfg_demo} ${OUTPUT}/DEPICT/cfg/${trait_name}.cfg 
cfg_file=${OUTPUT}/DEPICT/cfg/${trait_name}.cfg
sed -i "7c analysis_path: ${OUTPUT}/DEPICT/output" ${cfg_file}
sed -i "40c plink_executable: ${plink1_9}" ${cfg_file}
sed -i "13c gwas_summary_statistics_file: ${OUTPUT}/DEPICT/GWAS_input/${trait_name}_depict.clean.txt" ${cfg_file}
sed -i "19c label_for_output_files: ${trait_name}_clean" ${cfg_file}


python ${DEPICT} ${cfg_file}

