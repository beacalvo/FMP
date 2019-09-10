load("files_gtca.RData")
identical(as.character(all_ordered$ID),as.character(plink$V2))
all_ordered[,1] <- all_ordered[,2] 
all_ordered$fam_id <- ID   
colnames(all_ordered)[1:2] <- c("FIID","ID")

# File containing the phenotypes (metabolite levels)
pheno_file <- all_ordered[,c(1:2,204:247)]
write.table(pheno_file, file = "urine.phen", quote=FALSE, 
            row.names=FALSE, col.names=FALSE, sep = "\t")

# Covar: file containing categorical variables (sex and sample type)
covar_file <- all_ordered[,c(1:3,5:6)]
write.table(covar_file, file = "urine.covar", quote=FALSE, 
            row.names=FALSE, col.names=FALSE, sep = "\t")

# Qcovar: file containing quantitive variables (age and 20 PCs)
qcovar_file <- all_ordered[,c(1:2,4,7:26)]
write.table(qcovar_file, file = "urine.qcovar", quote=FALSE, 
            row.names=FALSE, col.names=FALSE, sep = "\t")

save.image("files_gtca.RData")
