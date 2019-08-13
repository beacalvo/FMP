# R script to run CMS

CMS_analysis <- function(pipeline_directory, phenotypes_file, summary_file, covariates_file, temp_genotypes_file, results_file) {
  
  # Import CMS function
  source(paste(pipeline_directory, "/others/Script2.3_CMS.R", sep=""))
  
  # Number of results
  results_number = 6 # betaMA, pvalMA, betaMC, pvalMC, Yrsq, Ncov
  
  #################################################################################################
  # INPUT
  
  phenotypes = read.table(file = phenotypes_file, na.strings=NA, sep = "\t", header=TRUE) # phenotypes
  covariates = as.matrix(read.table(file = covariates_file, na.strings=NA, sep = "\t", header=TRUE)) # covariates
  genotypes = read.table(file = temp_genotypes_file, na.strings=NA, sep = " ", header=TRUE) # genotypes
  
  summary = as.matrix(read.table(file = summary_file, header = TRUE, na.strings = NA, sep = ","))
  
  # Integrate information from summary file
  confounding_factors = summary[which(summary[,"Conf"] == 1), "Label"]
  outcomes = summary[which(summary[,"Outcome"] == 1), "Label"]
  
  SNPs =  colnames(genotypes)[-1] # SNPs rs numbers

  #################################################################################################
  # RUN CMS on each outcome and each SNP
  
  DATAMAT = merge(genotypes,phenotypes,by.x="IID", by.y = colnames(phenotypes)[1]) # merge data by IDs
  
  results_matrix = matrix(NA,length(SNPs),length(outcomes)*results_number+1) # results that will be written in the result file (one per block)
  results_matrix[,1] = SNPs # first column = SNP
  
  o = 1
  s = 1
  
  for (out in outcomes) { # Analyse each outcome
    
    covariates_list = covariates[,out]
    
    for (snp in SNPs) { # And each SNP
      
      results_CMS = MC(DATAMAT,out,snp,confounding_factors,covariates_list,NA,0,1,100,2,8,0.05,0)
      
      results_matrix[s,results_number*(o-1)+2] = as.numeric(results_CMS[1]) #betaMA
      results_matrix[s,results_number*(o-1)+3] = as.numeric(results_CMS[3]) #pvalMA
      results_matrix[s,results_number*(o-1)+4] = as.numeric(results_CMS[4]) #betaMC
      results_matrix[s,results_number*(o-1)+5] = as.numeric(results_CMS[6]) #pvalMC
      results_matrix[s,results_number*(o-1)+6] = as.numeric(results_CMS[7]) #Yrsq
      results_matrix[s,results_number*(o-1)+7] = as.numeric(results_CMS[8]) #Ncov
      
      s = s + 1
    }
    o = o + 1
    s = 1
  }
  
  #################################################################################################
  # WRITE RESULTS in text file
  write.table(results_matrix, file=results_file, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  
}
