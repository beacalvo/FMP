##############################################################################
# Recoding the genotype data from the binary format to the .raw format

plink --bfile HELIX_GWAS_FINAL_genotyped_XYchr --recodeA --not-chr 24,26 --maf 0.05 --out HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT_MAF0.05 # MAF 0.05


# Obtaining the BIM file
plink --bfile HELIX_GWAS_FINAL_genotyped_XYchr --make-bed --not-chr 24,26 --maf 0.05 --out HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT_MAF0.05

# Obtaining the FRQ file, which contains the MAF values
plink --bfile HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT_MAF0.05 --freq --out HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT_MAF0.05
