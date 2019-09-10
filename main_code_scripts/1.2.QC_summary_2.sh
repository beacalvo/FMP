# Creating histogram of MAF after QC

## Creating file containing all SNPs selected after post-imputation QC
## and their MAF

cd ~/data/WS_HELIX/HELIX_preproc/gwas/HRC_imp/QC_post_imp

## Selecting the SNPs name and the MAFs columns
awk '{print $1, $7}' good_SNPs/QC_lowfreq_Michigan_rsq05_MAF001_INMA_2019-01-29_Chr*.txt > post_imp_QC_BC/MAF_All.txt

cd post_imp_QC_BC/
  
mv MAF_All.txt MAF_All_after_filtering0.01.txt