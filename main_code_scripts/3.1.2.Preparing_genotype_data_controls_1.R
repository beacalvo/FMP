###### FADS ######

##############################################################################
# Recoding the genotype data from the binary format to the .raw format and 
# obtaining the FRQ file
plink --bfile HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT_MAF0.05 --from-kb 61300 --to-kb 61800  --chr 11 --recodeA --freq --out HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT_MAF0.05_FADS

# Obtaining the BIM file
plink --bfile HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT_MAF0.05 --from-kb 61300 --to-kb 61800  --chr 11 --make-bed --out HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT_MAF0.05_FADS


###### AGXT2 ######

##############################################################################
# Recoding the genotype data from the binary format to the .raw format and 
# obtaining the FRQ file
plink --bfile HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT_MAF0.05 --from-kb 34800 --to-kb 35200  --chr 5 --recodeA --freq --out HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT_MAF0.05_AGXT2

# Obtaining the BIM file
plink --bfile HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT_MAF0.05 --from-kb 34800 --to-kb 35200  --chr 5 --make-bed --out HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT_MAF0.05_AGXT2
