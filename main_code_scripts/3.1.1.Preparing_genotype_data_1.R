# Modifying the BIM file in order to add the MAF
maf <- read.table("HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT_MAF0.05.frq", header=TRUE)
maf_col <- maf$MAF
write.table(maf_col, file ="maf_column_MAF0.05.txt", row.names = FALSE, col.names = FALSE) 
