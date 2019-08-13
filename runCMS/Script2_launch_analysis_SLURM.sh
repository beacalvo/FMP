#!/bin/sh

#SBATCH --mem=6000
##SBATCH -o runCMS/test/output%a.log
##SBATCH -e runCMS/test/error%a.log

module load Python/3.5.2-foss-2016b
module load R/3.4.2-foss-2016b

pipeline_directory=$1
analysis_directory=$2
genotypes_file=$3
summary_file=$4
SNP_number=$5
block_size=$6

i=${SLURM_ARRAY_TASK_ID}

SNP1=`echo "${block_size}*(${i}-1)+1" | bc -l`
SNP2=`echo $((${block_size}*${i}<${SNP_number}?${block_size}*${i}:${SNP_number}))`

python3 ${pipeline_directory}/others/Script2.1_create_data_file.py ${pipeline_directory} ${analysis_directory} ${genotypes_file} ${summary_file} ${i} ${SNP1} ${SNP2}
