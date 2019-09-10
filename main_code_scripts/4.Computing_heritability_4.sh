##############################################################################
# ---- GCTA-GRM: calculating the genetic relationship matrix (GRM) from all 
#      the autosomal SNPs
gcta64 --bfile files_gcta/HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT_996_individuals \
--autosome --maf 0.05 \
--make-grm --out HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT_996_individuals \
--thread-num 10


# ---- GCTA-GREML analysis: estimating the variance explained by the SNPs

# Loop running GCTA-GREML for the 44 urine metabolites
for i in {1..44}
do gcta64 \
--grm HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT_996_individuals \
--pheno urine.phen --mpheno $i \
--covar urine.covar \
--qcovar urine.qcovar --reml \
--out HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT_996_individuals_$i \
--thread-num 10
done
