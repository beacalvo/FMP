##############################################################################
############# Comparison of our results with results in adults ###############
##############################################################################

##############################################################################
# ---- Combining results of urine papers ----

# -- Suhre et al (2011)
suhre <- read.table("/significant_snps_paper_Suhre.txt",
                    sep = "\t", header = TRUE)

# Selecting only those associations that are not with ratios
suhre <- suhre[suhre$Nominator == "ONE",]
suhre_final <- suhre

# Converting alleles to + strand
suhre_final$Allele.A[suhre$Allele.A == "A" & suhre$Strand == "-"] <- "T"
suhre_final$Allele.A[suhre$Allele.A == "T" & suhre$Strand == "-"] <- "A"
suhre_final$Allele.A[suhre$Allele.A == "C" & suhre$Strand == "-"] <- "G"
suhre_final$Allele.A[suhre$Allele.A == "G" & suhre$Strand == "-"] <- "C"

suhre_final$Allele.B[suhre$Allele.B == "A" & suhre$Strand == "-"] <- "T"
suhre_final$Allele.B[suhre$Allele.B == "T" & suhre$Strand == "-"] <- "A"
suhre_final$Allele.B[suhre$Allele.B == "C" & suhre$Strand == "-"] <- "G"
suhre_final$Allele.B[suhre$Allele.B == "G" & suhre$Strand == "-"] <- "C"

# Final dataset Suhre
suhre_final_2 <- suhre_final[,c("SNP.ID","CHR", "BP", "Enumerator", "BETA",
                                "P","Allele.A","Allele.B")]

alleles <- paste(suhre_final_2$Allele.A, suhre_final_2$Allele.B, sep = "/")

suhre_final_3 <- cbind(suhre_final_2[,c(1:6)], alleles)
write.table(suhre_final_3, 
            file = "../urine_papers_results/urine_results_Suhre_etal_2011.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")

# -- Rueedi et al (2014): directly copied into the excel file

# -- Raffler et al (2015): used the data in The GWAS metabolome database
#    for this paper (single metabolites and targeted and urine (SHIP))

raffle <- read.table("results_single_targeted_for_gwasserver_download.txt",
                     sep = "\t", header = TRUE)
raffle_final <- raffle[raffle$P < 5E-8,]

# Adding the HMDB ID
IDs <- read.table("IDs_equivalences_GWAS_metab_database.txt", sep = "\t",
                  header = TRUE)

IDs_simplifyed <- IDs[,c(1,2,4)]

raffle_final_HMDBid <- merge(raffle_final, IDs_simplifyed, by.x = "METID",
                             by.y = "Chenomx")

raffle_final_HMDBid <- raffle_final_HMDBid[,c("SNP", "CHR", "POS", "MET", 
                                              "BETA", "P", "EA", "HMDB")]

write.table(raffle_final_HMDBid, 
            file = "results_single_targeted_pvalue_filtered.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")