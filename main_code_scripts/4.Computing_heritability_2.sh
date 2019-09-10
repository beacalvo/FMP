# Creating PLINK files of only the 996 individuals for which there is 
# phenotypic information
plink --bfile ../HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT \
--keep ID_to_keep.txt --make-bed \
--out HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT_996_individuals