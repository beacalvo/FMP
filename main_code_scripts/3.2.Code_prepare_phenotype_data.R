##############################################################################
########### Script to prepare the data for the runCMS analysis ###############
##############################################################################
#
# 30/04/2019
# Beatriz Calvo

##############################################################################
# ---- Dummy variables

# Loading packages
library(fastDummies)

# Creating dummy variables from sex variable
dummy_sex <- dummy_columns(pdata_final_metadata_EUR$e3_sex.ser)
dummy_sex <- dummy_sex[,2]

# Creating dummy variables from sampling type variable
dummy_sample_type <- dummy_columns(pdata_final_metadata_EUR$Urine.SamplingType.uri)
dummy_sample_type <- dummy_sample_type[,2:3]

# Adding the dummy variables
pdata_final_metadata_EUR_dummy <- pdata_final_metadata_EUR
pdata_final_metadata_EUR_dummy$e3_sex.ser <- dummy_sex
colnames(pdata_final_metadata_EUR_dummy)[10] <- "sex.data_male"

colnames(dummy_sample_type) <- paste("sample_type", 
                                     colnames(dummy_sample_type), sep ="")

pdata_final_metadata_EUR_dummy <- cbind(pdata_final_metadata_EUR_dummy[,1:29], 
                                        dummy_sample_type, 
                                        pdata_final_metadata_EUR_dummy[,31:58])



# Testing individuals are in the same order in metadata and in expression data

identical(rownames(pdata_final_metadata_EUR_dummy),
          rownames(exprs_helixID_ordered)) #TRUE


# Obtaining "phenotypes file"
pdata_final_metadata_EUR_dummy_urine <- cbind("HelixID.ser" = 
                                                pdata_final_metadata_EUR$HelixID.ser, 
                                              "sex.data_male" = 
                                                pdata_final_metadata_EUR_dummy$sex.data_male,
                                              "age_sample_years" = 
                                                pdata_final_metadata_EUR_dummy$age_sample_years.ser, 
                                              dummy_sample_type, 
                                              pdata_final_metadata_EUR[,30:49])  

pdata_final_metadata_EUR_dummy_urine$HelixID.ser <- as.character(
  pdata_final_metadata_EUR_dummy_urine$HelixID.ser)


identical(pdata_final_metadata_EUR_dummy_urine$HelixID.ser, 
          rownames(exprs_helixID_ordered))  # TRUE: Same Helix_ID order

# Joining metabolites level data to the metadata: PHENOTYPES file

phenotypes_file_urine <- cbind(pdata_final_metadata_EUR_dummy_urine, 
                               exprs_helixID_ordered)


write.table(phenotypes_file_urine, "phenotypes_urine.txt", sep = "\t", quote =FALSE, row.names = FALSE)

# Saving workspace
save.image(file = "dummy_and_phenotypes_file_workspace.RData")

