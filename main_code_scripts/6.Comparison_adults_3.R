##############################################################################
# ---- Comparison of our results with results in adults from other studies


# -- Defining function that will be used to find common SNP-metabolite 
#    associations

# By-SNP

combine_results_by_snp <- function(results_snp,db) {
  
  # "results": our results // "db": results of other studies
  
  setkey(db, "SNP.ID_db")
  results_combined <- NULL
  results_combined <- as.data.frame(results_combined)
  met_combined <- NULL
  met_combined <- as.data.frame(met_combined)
  
  for (i in 1:nrow(results)) {
    snp_db <- as.data.frame(db[.(results$SNP[i]),])
    snp_db <- unique(snp_db)  # In "db" there are repeated rows
    
    if (is.na(snp_db$POS_db)) { # SNP not found in other results
      met_combined <- cbind(results[i,],"SNP not found in database",rep(NA,13))              
      names(met_combined) <- c(colnames(results), 
                               colnames(snp_db),
                               "Other_metabolite_db")
      results_combined <- rbind(results_combined, met_combined)
      met_combined <- NULL
      print(i)
    }
    
    else { 
      mets <- unlist(strsplit(as.character(results$HMDB[i]), ","))
      
      if (is.na(mets)) {
        met_combined <- cbind(results[i,], 
                              t(c("Missing metabolite ID in results", 
                                  rep(NA,13))))
        names(met_combined) <- c(colnames(results), 
                                 colnames(snp_db), 
                                 "Other_metabolite_db")
        results_combined <- rbind(results_combined, met_combined)
        met_combined <- NULL
        
      }
      
      else {
        
        for (met in mets) {
          index <- grep(met, snp_db$HMDB_db)
          
          if (length(index) == 0) {
            met_combined <- cbind(results[i,], snp_db[,1:3], "Other", 
                                                   snp_db[,5:13], 
                                                   snp_db$Metabolite_db)
            met_combined[j,17] <- as.character(snp_db[,"CHR_db"])
            met_combined[j,21] <- as.character(snp_db[,"P_db"])
            met_combined[j,24] <- as.character(snp_db[,"Region_db"])
            
          }
          
          else if (length(index) == 1) {
            met_combined <- cbind(results[i,], snp_db[index,], NA)
            met_combined[j,17] <- as.character(snp_db[index,"CHR_db"])
            met_combined[j,21] <- as.character(snp_db[index,"P_db"])
            met_combined[j,24] <- as.character(snp_db[index,"Region_db"])
          }
          
          else if (length(index) > 1) {
            for (j in 1:length(index)) {
              met_combined[j,1:(ncol(db)+ncol(results)+1)] <- as.data.frame(
                cbind(results[i,], snp_db[index[j],], NA))
              met_combined[j,1] <- as.character(results[i,1])
              met_combined[j,17] <- as.character(snp_db[index[j],"CHR_db"])
              met_combined[j,21] <- as.character(snp_db[index[j],"P_db"])
              met_combined[j,24] <- as.character(snp_db[index[j],"Region_db"])
            }
            
          }      
          names(met_combined) <- c(colnames(results), 
                                   colnames(snp_db), 
                                   "Other_metabolite_db")


          results_combined <- rbind(results_combined, as.data.frame(met_combined))
          met_combined <- NULL
          met_combined <- as.data.frame(met_combined)
          print(i)
        }
      }
    }    
  }
  
  return(results_combined)

}


# First, the column names of the data frame with the results of other studies 
# need to be changed so that they are differ from the column names of the
# data frame with our results


library(data.table)
db <- fread("urine_papers_results/Results_all_papers_urine_combined.txt", 
            header = T)

colnames(db) <- c("SNP.ID_db", "CHR_db", "POS_db", "Metabolite_db", "BETA_db",
                  "P_db", "Alleles_db", "HMDB_db", "Region_db", "MAF_db", 
                  "Author_db", "EA_db", "NEA_db")

results_snp <- read.table("summary_significant_snps_HMDB.txt", header = T)

results_combined_snp <- combine_results_by_snp(results_snp,db)
results_combined_snp <- results_combined_snp[!is.na(results_combined_snp$
                                                      Phenotype), ]


write.table(results_combined_snp, file= ".results_combined_snp_MAF0.05_.txt",
            col.names = T, quote = F,sep = "\t", row.names = F)


# By region

combine_results_by_region <- function(results_region,db) {
  
  # "results": our results // "db": results of other studies
  
  setkey(db, "Region_db")
  results_combined <- NULL
  results_combined <- as.data.frame(results_combined)
  met_combined <- NULL
  met_combined <- as.data.frame(met_combined)
  
  for (i in 1:nrow(results_region)) {
    region_db <- as.data.frame(db[.(results_region$Region[i]),])
    region_db <- unique(region_db)  # In "db" there are repeated rows
    
    if (is.na(region_db$POS_db)) { # Region not found in other results
      met_combined <- cbind(results_region[i,], "Region not found in database", 
                            rep(NA,13))
      names(met_combined) <- c(colnames(results_region), 
                               colnames(region_db),
                               "Other_metabolite_db",
                               "Number_sign_SNPs_in_region_db")
      results_combined <- rbind(results_combined, met_combined)
      met_combined <- NULL
      print(i)
    }
    
    else { 
      mets <- unlist(strsplit(as.character(results_region$HMDB[i]), ","))
      
      if (is.na(mets)) {
        met_combined <- cbind(results_region[i,],
                              t(c("Missing metabolite ID in results", 
                                                      rep(NA,14))))
        names(met_combined) <- c(colnames(results_region), 
                                 colnames(region_db), 
                                 "Other_metabolite_db",
                                 "Number_sign_SNPs_in_region_db")
        results_combined <- rbind(results_combined, met_combined)
        met_combined <- NULL
        
      }
      
      else {
        
        for (met in mets) {
          index <- grep(met, region_db$HMDB_db)
          
          if (length(index) == 0) {
            met_combined <- cbind(results_region[i,], region_db[,1:3], "Other", 
                                  region_db[,5:13], 
                                  region_db$Metabolite_db)
            met_combined[j,17] <- as.character(region_db[,"CHR_db"])
            met_combined[j,21] <- as.character(region_db[,"P_db"])
            met_combined[j,24] <- as.character(region_db[,"Region_db"])
            
          }
          
          else {
            met_combined <- cbind(results_region[i,], region_db[index[1],], 
                                  NA, length(index))
            met_combined[j,17] <- as.character(region_db[index[1],"CHR_db"])
            met_combined[j,21] <- as.character(region_db[index[1],"P_db"])
            met_combined[j,24] <- as.character(region_db[index[1],"Region_db"])
            met_combined[j,30] <- as.character(length(index))
          }
          
          names(met_combined) <- c(colnames(results_region), 
                                   colnames(region_db), 
                                   "Other_metabolite_db",
                                   "Number_sign_SNPs_in_region_db")
          
          results_combined <- rbind(results_combined, 
                                    as.data.frame(met_combined))
          met_combined <- NULL
          met_combined <- as.data.frame(met_combined)
          print(i)
        }
      }
    }    
  }
  
  return(results_combined)
  
}


results_region <- read.table("summary_significant_regions_HMDB.txt", header = T)
results_combined_region_2 <- combine_results_by_region(results_region,db)
results_combined_region <- results_combined_region[!is.na(results_combined_region$Phenotype), ]

write.table(results_combined_region, file= "../results_combined_region_MAF0.05.txt", 
            col.names = T, quote = F,sep = "\t", row.names = F)
write.table(results_combined_region_2, file= "../results_combined_region_MAF0.05_2.txt", 
            col.names = T, quote = F,sep = "\t", row.names = F)
