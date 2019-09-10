## Histogram MAF after QC
maf_filtered <- read.table("MAF_All_after_filtering0.01.txt", 
                           header= T , sep = " ")
maf_numeric <- as.numeric(as.character(maf_filtered$MAF))
maf_numeric_NA <- maf_numeric[!is.na(maf_numeric)]

pdf("hist_MAF_after_post_imputation_QC.pdf")
hist(maf_numeric_NA, xlab = "MAF", col = "red", 
     main = "Histogram MAF after post-imp QC")

# ---- Histogram of R2

# Joining all the values of R2 of all chromosomes into a single vector
library(data.table)
R2 <- NULL
for (i in 1:22){
  dataInfo <- fread(paste("~/data/WS_HELIX/HELIX_preproc/gwas/HRC_imp/Final_data_Michigan/chr",
                          i,".info",sep=""),
                    data.table=F,stringsAsFactors = F)
  R2 <- c(R2,dataInfo$Rsq)
}
pdf("hist_R2_accuracy_imputation.pdf")
hist(R2, xlab = "R2", col = "gray", main = "Histogram of the R2 accuracy of imputation")



# ---- Plot correlation MAF and R2

# Only 10000 values of MAF and R2 were used for the plot in order to see more
# clearly if there was any correlation
R2 <- NULL
MAF <- NULL
for (i in 1:22){
  dataInfo <- fread(paste("~/data/WS_HELIX/HELIX_preproc/gwas/HRC_imp/Final_data_Michigan/chr",
                          i,".info",sep=""),data.table=F,stringsAsFactors = F)
  R2 <- c(R2,dataInfo$Rsq)
  MAF <- c(MAF,dataInfo$MAF)
}
vec <- round(runif(10000, min =1,  max = length(MAF)))
png("plot_MAF_vs_R2.png")
plot(R2[vec],MAF[vec], main = "Plot of imputation accuracy (R2) vs MAF",
     xlab = "Imputation accuracy (R2)", ylab = "MAF", pch = 20)
