# ---- Script 1

# Folder path for pipeline
arg1=/home/isglobal.lan/bcalvo/data/WS_HELIX/HELIX_analyses/metabQLTs_BC/CMS_algorithm/runCMS_approach/runCMS-master

# Folder for analysis
arg2=/home/isglobal.lan/bcalvo/data/WS_HELIX/HELIX_analyses/metabQLTs_BC/CMS_algorithm/runCMS_approach/urine_genotyped/urine_genotyped_MAF0.05_PC_adjusted

# Input phenotypes file
arg3=/home/isglobal.lan/bcalvo/data/WS_HELIX/HELIX_analyses/metabQLTs_BC/CMS_algorithm/runCMS_approach/input_files/phenotype_data_in_log_PCs/phenotypes_urine.txt

# Input summary file
arg4=/home/isglobal.lan/bcalvo/data/WS_HELIX/HELIX_analyses/metabQLTs_BC/CMS_algorithm/runCMS_approach/input_files/summary_urine.csv

# Covariates number
arg5=30

# Threshold on outcome variance explained by covariates
arg6=0.7

Rscript /home/isglobal.lan/bcalvo/data/WS_HELIX/HELIX_analyses/metabQLTs_BC/CMS_algorithm/runCMS_approach/runCMS-master/Script1_data_preparation.R $arg1 $arg2 $arg3 $arg4 $arg5 $arg6



# ---- Script 2

# Folder path for pipeline
arg1=/homes/users/bcalvo/CMS/runCMS_backup

# Folder for analysis
arg2=/homes/users/bcalvo/CMS/urine_genotyped/urine_genotyped_MAF0.05_PC_adjusted

# Genotypes file
arg3=/homes/users/bcalvo/shared_data/GWAS_data/MAF_0.05/HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT_MAF0.05.raw

# Input summary file
arg4=/homes/users/bcalvo/CMS/input_data/summary_urine_PCs.csv

# SNP number
arg5=283704

# Number of SNPs per block
arg6=1000

sbatch --array=1-284 /homes/users/bcalvo/CMS/runCMS_backup/Script2_launch_analysis_SLURM.sh $arg1 $arg2 $arg3 $arg4 $arg5 $arg6




# ---- Script 3

# Folder path for pipeline
arg1=/home/isglobal.lan/bcalvo/data/WS_HELIX/HELIX_analyses/metabQLTs_BC/CMS_algorithm/runCMS_approach/example_test_runCMS/runCMS_new

# Folder for analysis
arg2=/home/isglobal.lan/bcalvo/data/WS_HELIX/HELIX_analyses/metabQLTs_BC/CMS_algorithm/runCMS_approach/urine_genotyped/urine_genotyped_MAF0.05_PC_adjusted
  
# Genotypes file
arg3=/home/isglobal.lan/bcalvo/data/WS_HELIX/HELIX_analyses/metabQLTs_BC/CMS_algorithm/runCMS_approach/urine_genotyped/HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT_MAF0.05_MAF.bim

# Input summary file
arg4=/home/isglobal.lan/bcalvo/data/WS_HELIX/HELIX_analyses/metabQLTs_BC/CMS_algorithm/runCMS_approach/input_files/summary_urine_PCs.csv

# Blocks number
arg5=284

# Significance threshold

arg6=0.000000001522778    # 1.522778e-09

python3 /home/isglobal.lan/bcalvo/data/WS_HELIX/HELIX_analyses/metabQLTs_BC/CMS_algorithm/runCMS_approach/example_test_runCMS/runCMS_new/Script3_results_analysis.py $arg1 $arg2 $arg3 $arg4 $arg5 $arg6
