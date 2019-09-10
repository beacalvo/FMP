# ---- Script 1

#Folder path for pipeline
arg1=/home/isglobal.lan/bcalvo/data/WS_HELIX/HELIX_analyses/metabQLTs_BC/CMS_algorithm/runCMS_approach/runCMS-master

#Folder for analysis
arg2=/home/isglobal.lan/bcalvo/data/WS_HELIX/HELIX_analyses/metabQLTs_BC/CMS_algorithm/runCMS_approach/FADS_runCMS/urine/log2_typification

#Input phenotypes file
arg3=/home/isglobal.lan/bcalvo/data/WS_HELIX/HELIX_analyses/metabQLTs_BC/CMS_algorithm/runCMS_approach/input_files/phenotype_data_in_log/phenotypes_urine.txt

#Input summary file
arg4=/home/isglobal.lan/bcalvo/data/WS_HELIX/HELIX_analyses/metabQLTs_BC/CMS_algorithm/runCMS_approach/input_files/summary_urine.csv

#Covariates number
arg5=30

#Threshold on outcome variance explained by covariates
arg6=0.7

Rscript /home/isglobal.lan/bcalvo/data/WS_HELIX/HELIX_analyses/metabQLTs_BC/CMS_algorithm/runCMS_approach/runCMS-master/Script1_data_preparation.R $arg1 $arg2 $arg3 $arg4 $arg5 $arg6


# ---- Script 2

#Folder path for pipeline
arg1=/homes/users/bcalvo/CMS/runCMS

#Folder for analysis
arg2=/homes/users/bcalvo/CMS/AGXT2/urine/log2_typification

#Genotypes file
arg3=/homes/users/bcalvo/CMS/AGXT2/HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT_MAF0.05_AGXT2.raw

#Input summary file
arg4=/homes/users/bcalvo/CMS/input_data/summary_urine.csv

#SNP number
arg5=1065

#Number of SNPs per block
arg6=22

sbatch --array=1-50 /homes/users/bcalvo/CMS/runCMS/Script2_launch_analysis_SLURM.sh $arg1 $arg2 $arg3 $arg4 $arg5 $arg6


# ---- Script 3

#Folder path for pipeline
arg1=/home/isglobal.lan/bcalvo/data/WS_HELIX/HELIX_analyses/metabQLTs_BC/CMS_algorithm/runCMS_approach/example_test_runCMS/runCMS

#Folder for analysis
arg2=/home/isglobal.lan/bcalvo/data/WS_HELIX/HELIX_analyses/metabQLTs_BC/CMS_algorithm/runCMS_approach/AGXT2/urine/log2_typification

#Genotypes file
arg3=/home/isglobal.lan/bcalvo/data/WS_HELIX/HELIX_analyses/metabQLTs_BC/CMS_algorithm/runCMS_approach/AGXT2/HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT_MAF0.05_AGXT2_MAF.bim

#Input summary file
arg4=/home/isglobal.lan/bcalvo/data/WS_HELIX/HELIX_analyses/metabQLTs_BC/CMS_algorithm/runCMS_approach/input_files/summary_urine.csv

#Blocks number
arg5=50


#Significance threshold

arg6=0.00000000056 ## 5.594062e-10

python3 /home/isglobal.lan/bcalvo/data/WS_HELIX/HELIX_analyses/metabQLTs_BC/CMS_algorithm/runCMS_approach/example_test_runCMS/runCMS/Script3_results_analysis.py $arg1 $arg2 $arg3 $arg4 $arg5 $arg6