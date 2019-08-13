###########################################################################################################################################
# Function to compute pairwise r2 between phenotypes                                                                                      #
# INPUT: Matrix of phenotypes (One phenotype per column, one line per individual), threshold to filter pairwise r2                        #
# OUTPUT: a matrix with pairwise r2 values for each pair of phenotypes, with NA values on the diagonal and for each value above threshold #
###########################################################################################################################################
r2_computation <- function(phenos, outcomes, covariates, excluded_covs, r2_threshold) {
  
  r2 = matrix(NA, length(covariates), length(outcomes)) # matrix of pairwise r2
  colnames(r2) = outcomes
  rownames(r2) = covariates
  
  for (i in covariates) {
    
    for (j in outcomes) {
      
      if (i != j) {
        to_exclude = strsplit(as.character(excluded_covs[,j]), ";")[[1]]
        if (!(i %in% to_exclude)) {
          reg = lm(phenos[,j] ~ phenos[,i]) # perform linear regression for each pair 
          rsq = summary(reg)$r.squared # get r2 value from this regression
      
          if(rsq < r2_threshold) { # keep only r2 values below the threshold (others will be NA)
            r2[i,j] = rsq
          }
        }
      }
    }
  }
  
  return(r2)
}


################################################################################################################################
# Function to pre-select the covariates                                                                                        #
# INPUT: Matrix of phenotypes (One phenotype per column, one line per individual), matrix of pairwise r2, covariates number   #
# OUTPUT: Dataframe with one column per phenotype with the pre-selected covariates                                             #
################################################################################################################################
covariates_pre_selection <- function(phenos, outcomes, covariates, r2, covariates_number) {
  
  covs = matrix(NA, covariates_number, length(outcomes)) # data frame of pre-selected covariates
  colnames(covs) = outcomes # one column per outcome
  
  for (i in outcomes) {

    aic = rep(NA,covariates_number) # AIC scores
    
    covs[1, i] = names(which.max(r2[, i])) # first pre-selected covariate = the most correlated
    aic[1] = AIC(lm(phenos[, i] ~ phenos[, covs[1, i]])) # first AIC score
    
    # Function to compute the AIC difference between last pre-selected covariates and remaining phenotypes
    f_aic <- function(phe, c) { 
      return(aic[c - 1] - AIC(lm(as.formula(paste("phenos[, i] ~ ", paste("phenos[, covs[", 1:(c - 1), ", i]]", sep = "", collapse = "+"), "+ phe")))))
    }
    
    if (covariates_number > 1) {
      for (c in 2:covariates_number) {
        
        # Exclude phenotypes indexes corresponding to already pre-selected covariates and NA value in r2 matrix
        remaining_covariates = covariates[!(covariates %in% unique(c(covs[1:(c - 1), i], names(which(is.na(r2[, i]))))))]
        
        # Choose covariate which bring the more new information (maximal AIC difference between current model with and without this variable)
        covs[c, i] = names(which.max(apply(phenos[, remaining_covariates], 2, f_aic, c)))
        
        # Compute new AIC score
        aic[c] = AIC(lm(as.formula(paste("phenos[, i] ~ ", paste("phenos[, covs[", 1:c, ", i]]", sep = "", collapse = "+")))))
        
      }
    }
  }
  
  return(covs)
}

