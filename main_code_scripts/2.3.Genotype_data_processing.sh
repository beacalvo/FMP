##############################################################################
# This command is used to change the codification of the sexual chromosomes 
# and mitochondrial DNA from the chrX, chrY and MT nomenclature to the
# numerical nomenlcature in PLINK (23, 24 and 26, respectively)

plink --bfile /home/isglobal.lan/bcalvo/data/WS_HELIX/HELIX_preproc/gwas/Final_data/v2/HELIX_GWAS_FINAL  --make-bed --out HELIX_GWAS_FINAL_genotyped

# Moving the variants located in the PAR1 and PAR2 regions to a different
# chromosome (chromosome number 25 in PLINK)
plink --bfile HELIX_GWAS_FINAL_genotyped --split-x 'hg19' --make-bed --out HELIX_GWAS_FINAL_genotyped_XYchr


# Removing variants located in mitochondrial DNA and in the non-canonical PAR region of chromosome Y
plink --bfile HELIX_GWAS_FINAL_genotyped_XYchr --make-bed --not-chr 24,26 --out HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT
