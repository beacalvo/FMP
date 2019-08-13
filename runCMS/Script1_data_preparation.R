# R script to prepare data before CMS analysis
# PHENOTYPES: perform rank-transformation to normality and pre-select covariates

#################################################################################################
# Get parameters from command line

args = commandArgs(trailingOnly=TRUE)
if (length(args) == 6) {
  pipeline_directory = args[1]
  analysis_directory = args[2]
  input_phenotypes_file = args[3]
  input_summary_file = args[4]
  covariates_number = as.integer(args[5])
  r2_threshold = as.numeric(args[6])
  
  if (!(dir.exists(pipeline_directory))) {
    stop("The pipeline folder does not exist")
  }
  if (!(dir.exists(analysis_directory))) {
    stop("The analysis folder does not exist")
  }
  if (!(file.exists(input_phenotypes_file))) {
    stop("The phenotypes file does not exist")
  } 
  if (!(file.exists(input_summary_file))) {
    stop("The summary file does not exist")
  }
  
} else {
  stop("The number of arguments is not correct.
Make sure you launch this script with these 7 arguments:
  1) = Folder path for pipeline  
  2) = Folder for analysis  
  3) = Input phenotypes file  
  4) = Input summary file  
  5) = Covariates number  
  6) = Threshold on outcome variance explained by covariates")
}

source(paste(pipeline_directory, "/others/Script1_functions.R", sep = ""))

output_phenotypes_file = paste(analysis_directory, "/phenotypes_for_analysis.txt", sep="")
output_covariates_file = paste(analysis_directory, "/covariates.txt"             , sep="")
log_file               = paste(analysis_directory, "/log_Script1.txt"            , sep="")
file.create(log_file)

#################################################################################################
# Phenotypes treatment

# Read files
phenos  = read.table(file = input_phenotypes_file, header = TRUE, na.strings = NA, sep = "\t")
summary = read.table(file = input_summary_file   , header = TRUE, na.strings = NA, sep = ",")


# Get classification from summary file
conf         = as.matrix(summary)[which(summary$Conf         == 1), "Label"]
outcomes     = as.matrix(summary)[which(summary$Outcome      == 1), "Label"]
covariates   = as.matrix(summary)[which(summary$Covariate    == 1), "Label"]
excluded     = as.matrix(summary)[which(summary$Outcome      == 1), "Excluded"]

excluded_covs = data.frame(t(excluded))
colnames(excluded_covs) = outcomes

write("Confounding factors:"            , file = log_file)
write(paste(conf, collapse = ", ")      , file = log_file, append = TRUE)
write("\nOutcomes:"                     , file = log_file, append = TRUE)
write(paste(outcomes, collapse = ", ")  , file = log_file, append = TRUE)
write("\nCovariates:"                   , file = log_file, append = TRUE)
write(paste(covariates, collapse = ", "), file = log_file, append = TRUE)

# Check if some variables are confounding factors AND outcomes
pb_1 = intersect(conf, outcomes)
if(!setequal(pb_1,c())) {
  stop(paste("These variables cannot be confounding factors and outcomes at the same time:\n", paste(pb_1, collapse = "\n"), sep = ""))
}

# Check if some variables are confounding factors AND covariates
pb_2 = intersect(conf, covariates)
if(!setequal(pb_2, c())) {
  stop(paste("These variables cannot be confounding factors and covariates at the same time:\n", paste(pb_2, collapse = "\n"), sep = ""))
}

# Standardization of all phenotypes
for (i in unique(c(conf, outcomes, covariates))) {
  phenos[,i] = (phenos[,i] - mean(phenos[,i] , na.rm = TRUE)) / sqrt(var(phenos[, i] , na.rm = TRUE))
}

# Exportation
write.table(phenos, file = output_phenotypes_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write("\nStandardized phenotypes exported to file \"phenotypes_for_analysis.txt\"", file = log_file, append = TRUE)

#################################################################################################
# Covariates pre-selection

print("Covariates pre-selection...")
write(paste("\nSelection of", covariates_number, "covariates per outcome...", sep = " "), log_file, append = TRUE)

r2 = r2_computation(phenos, outcomes, covariates, excluded_covs, r2_threshold)
pre_selected_covariates = covariates_pre_selection(phenos, outcomes, covariates, r2, covariates_number)

print("Done")
write.table(pre_selected_covariates, file = output_covariates_file, row.names = FALSE, col.names = TRUE, quote = FALSE, na = "NA", sep = "\t")

write("Covariates exported to file \"covariates.txt\"", file = log_file, append = TRUE)

#################################################################################################
# Prepare folder analysis

dir.create(file.path(analysis_directory, "tmp"), showWarnings = FALSE)
dir.create(file.path(analysis_directory, "results_per_block"), showWarnings = FALSE)
dir.create(file.path(analysis_directory, "results_summary"), showWarnings = FALSE)

