##############################################################################
# ---- Creating Phenotypes summary file for urine

phenotypes_summary_urine <- as.data.frame(matrix(0, ncol = 5,
                                                 nrow =  dim(phenotypes_file_urine)[2]))

colnames(phenotypes_summary_urine) <- c("Label", "Conf", "Outcome", 
                                        "Covariate", "Excluded")

phenotypes_summary_urine$Label <- colnames(phenotypes_file_urine)

# Confusing variables
phenotypes_summary_urine$Conf[1:24] <- 1

# Outcome variables
phenotypes_summary_urine$Outcome[202:245] <- 1

# Covariate variables
phenotypes_summary_urine$Covariate[25:245] <- 1

# Adding the "Excluded" metabolites: for metabolites found both in serum and 
# in urine, the metabolite in urine or serum, respectively, is excluded (not 
# considered as a possible covariate)

phenotypes_summary_urine$Excluded <- ""

metabolites_names <- c(gsub(".serum", "", phenotypes_summary_urine$Label[25:201]), 
                       gsub(".urine", "", phenotypes_summary_urine$Label[202:245]))

common_metab <- metabolites_names[duplicated(metabolites_names)]

for (i in common_metab) {
  position <- which(phenotypes_summary_urine$Label == paste(common_metab,
                                                            "urine", sep ="."))
  
  phenotypes_summary_urine[position, "Excluded"] <- paste(common_metab, 
                                                          "serum", sep =".") 
}

# Excluding "X3.aminoisobutyrate.urine" from being a possible covariate of 
# "X3.hydroxybutyrate.3.aminoisobutyrate.urine" and vice versa, since 
# "X3.hydroxybutyrate.3.aminoisobutyrate.urine" represents the concentrations 
# of both "X3.aminoisobutyrate.urine" and "3.hydroxybutyrate"

position <- which(phenotypes_summary_urine$Label == "X3.aminoisobutyrate.urine")
phenotypes_summary_urine[position, "Excluded"] <- 
  "X3.hydroxybutyrate.3.aminoisobutyrate.urine" 

position <- which(phenotypes_summary_urine$Label == 
                    "X3.hydroxybutyrate.3.aminoisobutyrate.urine")
phenotypes_summary_urine[position, "Excluded"] <- "X3.aminoisobutyrate.urine"


write.csv(phenotypes_summary_urine, file = "summary_urine.csv", 
          row.names = FALSE, quote = FALSE)

save(phenotypes_summary_urine, file = "summary_urine.RData")

