##############################################################################
################ R script to summarize results of variant QC #################
##############################################################################
#
# 25/02/2019
# Beatriz Calvo

##############################################################################
############ Counting the number of variants excluded during QC ##############

# Loading packages
library(data.table)

##############################################################################
# Start of the script
for (i in 1:22){
  
# Reading data for each chromosome (22 chromosomes in total) from the files 
# generated with vcftools and from ".info" file containing information about 
# the imputation accuracy
  
dataFreq <- fread(
  paste("~/data/WS_HELIX/HELIX_preproc/gwas/HRC_imp/QC_post_imp/MAF_filtering/freq_chr",
        i,".frq",sep=""))

dataInfo <- fread(
  paste("~/data/WS_HELIX/HELIX_preproc/gwas/HRC_imp/Final_data_Michigan/chr",
        i,".info",sep=""),data.table=F,stringsAsFactors = F)

dataHWE <- fread(
  paste("~/data/WS_HELIX/HELIX_preproc/gwas/HRC_imp/QC_post_imp/HWE_pvalue/hwe_chr",
        i,".hwe",sep=""),stringsAsFactors =F)

##############################################################################
# Editing "Freq" File

## Naming columns
names(dataFreq)[names(dataFreq)=="V1"] <- "CHR"
names(dataFreq)[names(dataFreq)=="CHROM"] <- "Pos"
names(dataFreq)[names(dataFreq)=="{ALLELE:FREQ}"] <- "Freq"

## Make dataframe of wanted columns
dataFreq <- data.frame(dataFreq$CHR,dataFreq$Pos,dataFreq$Freq)

## Calculate MAF from counts
ALT_ALLELE_FREQ <- as.numeric(
  (unlist(strsplit(as.character(dataFreq$dataFreq.Freq),
                   ":")))[c(FALSE,TRUE)])

ALT_ALLELE <- as.character(
  (unlist(strsplit(as.character(dataFreq$dataFreq.Freq),":")))[c(TRUE,FALSE)])

MAF <- ifelse(as.numeric(ALT_ALLELE_FREQ)>0.5,1-as.numeric(ALT_ALLELE_FREQ),
              as.numeric(ALT_ALLELE_FREQ))

## Make CHR:POS column
SNP <-as.character(paste(dataFreq$dataFreq.CHR,":",
                         dataFreq$dataFreq.Pos,sep=""))

## Store all information in dataframe
dataFreqClean <- data.frame(SNP,dataFreq$dataFreq.CHR,dataFreq$dataFreq.Pos,
                            ALT_ALLELE_FREQ,MAF,ALT_ALLELE)
dataFreqClean$SNP <- as.character(dataFreqClean$SNP)
colnames(dataFreqClean) <- c("SNP","CHR","POS","ALT_ALLELE_FREQ",
                             "MAF","ALT_ALLELE")

##############################################################################
# Editing "Info" File

## Make dataframe of wanted columns
dataInfoClean <- data.frame(dataInfo$SNP,dataInfo[,"REF(0)"],
                            dataInfo[,"ALT(1)"],dataInfo$Rsq,
                            dataInfo$Genotyped)
names(dataInfoClean)[names(dataInfoClean)=="dataInfo.SNP"] <- "SNP"
names(dataInfoClean)[names(dataInfoClean)=="dataInfo.Rsq"] <- "Rsq"

dataInfoClean$SNP <- as.character(dataInfoClean$SNP)
colnames(dataInfoClean) <- c("SNP","REF_ALLELE","ALT_ALLELE","Rsq","SOURCE")

##############################################################################
# Editing "HWE" File

## Make dataframe of wanted columns
SNP <- as.character(paste(dataHWE$CHR,":",dataHWE$POS,sep=""))

dataHWEClean <- data.frame(SNP,dataHWE$P_HWE)

##############################################################################
# Merging "Freq" and "Info" Files
dataMerge1 <- merge(dataInfoClean,dataFreqClean,by = c("SNP","ALT_ALLELE"))

##############################################################################
# Deleting multiallelic SNPs
dataMerge1sort <- dataMerge1[order(SNP),]

dataMerge2 <- data.frame(dataMerge1sort$SNP,dataMerge1sort$CHR,
                         dataMerge1sort$POS,dataMerge1sort$REF_ALLELE,
                         dataMerge1sort$ALT_ALLELE,
                         dataMerge1sort$ALT_ALLELE_FREQ,
                         dataMerge1sort$MAF,dataMerge1sort$Rsq,
                         dataMerge1sort$SOURCE)

colnames(dataMerge2)<- c("SNP","CHR","POS","REF_ALLELE","ALT_ALLELE",
                         "ALT_ALLELE_FREQ","MAF","Rsq","SOURCE")

toRemove <-names(which(table(dataMerge2$SNP)>1))

dataMerge1Clean <-dataMerge2[ !dataMerge2$SNP %in% toRemove ,  ]

dataHWEClean <- dataHWEClean[ !dataHWEClean$SNP %in% toRemove,]

##############################################################################
# Merging combined "Freq" & "Info" File with "HWE" File
dataTotal <- merge(dataMerge1Clean, dataHWEClean, by = "SNP")

colnames(dataTotal)<- c("SNP","CHR","POS","REF_ALLELE","ALT_ALLELE",
                        "ALT_ALLELE_FREQ","MAF","Rsq","SOURCE","HWE_PValue")

##############################################################################
# Filtering Files

## Setting Date
date <- Sys.Date()

dataMAF001_R2_0.25 <- dataTotal[(dataTotal$Rsq > 0.25) & 
                                  (dataTotal$MAF > 0.01) &
                                  (dataTotal$HWE_PValue > 10e-6),]

dataMAF001_R2_0.5 <- dataTotal[(dataTotal$Rsq > 0.5) &
                                 (dataTotal$MAF > 0.01) &
                                 (dataTotal$HWE_PValue > 10e-6),]

dataMAF001_R2_0.80 <- dataTotal[(dataTotal$Rsq > 0.80) &
                                  (dataTotal$MAF > 0.01) &
                                  (dataTotal$HWE_PValue > 10e-6),]

resume_table <- data.frame("Initial_variants" = dim(dataFreq)[1], 
                           "Multiallelic_SNPs"= length(toRemove) , 
                           "R2_0.25_included" = sum(dataTotal$Rsq > 0.25), 
                           "R2_0.25_excluded" = sum(dataTotal$Rsq <= 0.25), 
                           "R2_0.5_included" = sum(dataTotal$Rsq > 0.5), 
                           "R2_0.5_excluded" = sum(dataTotal$Rsq <= 0.5), 
                           "R2_0.80_included" = sum(dataTotal$Rsq > 0.80), 
                           "R2_0.80_excluded" = sum(dataTotal$Rsq <= 0.80), 
                           "MAF0.01_included" = sum(dataTotal$MAF > 0.01), 
                           "MAF0.01_excluded" = sum(dataTotal$MAF <= 0.01), 
                           "HWE_10e-6_included" = 
                             sum(dataTotal$HWE_PValue > 10e-6), 
                           "HWE_10e-6_excluded" = 
                             sum(dataTotal$HWE_PValue <= 10e-6), 
                           "Final_variantsR2_0.25" = 
                             dim(dataMAF001_R2_0.25)[1],
                           "Final_variantsR2_0.5" = dim(dataMAF001_R2_0.5)[1],
                           "Final_variantsR2_0.80" = 
                             dim(dataMAF001_R2_0.80)[1]
                           )


##############################################################################
# Creating a table for each chromosome with the number of variants that 
# are filtered and selected

write.table(
  resume_table,
  paste("~/data/WS_HELIX/HELIX_preproc/gwas/HRC_imp/QC_post_imp/post_imp_QC_BC/number_variantsR0.80",
        date,"_Chr",i,".txt",sep=""),quote=F,row.names=F,col.names=T)
}

##############################################################################
# Merging the information found in separate files for each chromosome 
# into a single file
total_initial <- 0
multiallelic_excluded <- 0
R2_0.5_excluded <- 0
R2_0.5_included <- 0
R2_0.25_excluded <- 0
R2_0.25_included <- 0
R2_0.80_excluded <- 0
R2_0.80_included <- 0
HWE_excluded <- 0
HWE_included <- 0
MAF_excluded <- 0
MAF_included <- 0
total_finalR2_0.25 <- 0
total_finalR2_0.5 <- 0
total_finalR2_0.80 <- 0

for (i in 1:22){
  x <- read.table(
    paste("~/data/WS_HELIX/HELIX_preproc/gwas/HRC_imp/QC_post_imp/post_imp_QC_BC/number_variantsR0.802019-03-02_Chr",
          i,".txt",sep=""), header = T)
  total_initial <- total_initial + x[1,1]
  multiallelic_excluded <- multiallelic_excluded + x[1,2]
  R2_0.25_included <- R2_0.25_included + x[1,3]
  R2_0.25_excluded <- R2_0.25_excluded + x[1,4]
  R2_0.5_included <- R2_0.5_included + x[1,5]
  R2_0.5_excluded <- R2_0.5_excluded + x[1,6]
  R2_0.80_included <- R2_0.80_included + x[1,7]
  R2_0.80_excluded <- R2_0.80_excluded + x[1,8]
  MAF_included <- MAF_included + x[1,9]
  MAF_excluded <- MAF_excluded + x[1,10]
  HWE_included <- HWE_included + x[1,11]
  HWE_excluded <- HWE_excluded + x[1,12]
  total_finalR2_0.25 <- total_finalR2_0.25 + x[1,13]
  total_finalR2_0.5 <- total_finalR2_0.5 + x[1,14]
  total_finalR2_0.80 <- total_finalR2_0.80 + x[1,15]
}

resume_all <- data.frame(total_initial, multiallelic_excluded, 
                         R2_0.25_included, R2_0.25_excluded, 
                         R2_0.5_included, R2_0.5_excluded, 
                         R2_0.80_included, R2_0.80_excluded, 
                         MAF_included, MAF_excluded, HWE_included, 
                         HWE_excluded, total_finalR2_0.25, 
                         total_finalR2_0.5, total_finalR2_0.80)

# Writing the file that contains all the summarized information
write.table(resume_all, 
            file = "~/data/WS_HELIX/HELIX_preproc/gwas/HRC_imp/QC_post_imp/post_imp_QC_BC/resume_all_differentR_withR0.80.txt",
            quote=F,row.names=F,col.names=T)



##############################################################################
###################### Creating plots of R2 and MAF ########################

# ---- Histogram of MAF

# Creating histogram of MAF before QC
pdf("hist_MAF_before_post_imputation_QC.pdf") 
hist(MAF, xlab = "MAF", col = "red", main ="Histogram MAF before post-imp QC")

# Storing the MAF values in a file for future uses
write.table(MAF, "MAF_before_post_imputation_QC.txt")
