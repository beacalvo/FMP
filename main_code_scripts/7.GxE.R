# ---- Creating dataset with dietary variables to be used

expos <- expos(exppostnatal_raw)

expos_short <- cbind(helixID = sampleNames(exppostnatal_raw), 
                     meat = expos$hs_total_meat_Ter, 
                     beverages = expos$hs_beverages_Ter, 
                     dairy = expos$hs_dairy_Ter, 
                     fish = expos$hs_total_fish_Ter)
	
pNames <- phenotypeNames(exppostnatal_raw)

helixID_sampleID <- cbind(helixID=pdata$HelixID, sampleID= pdata$SampleID,
                          obesity = pdata$hs_bmicat_None)

expos_helixID <- merge(helixID_sampleID, expos_short, by.x = "sampleID",
                       by.y = "helixID")
expos_helixID <- expos_helixID[,-1]
expos_helixID <- expos_helixID[complete.cases(expos_helixID),]

# Creating dummy variables from dietary intake frequency variables in tertiles
obesity <- dummy_cols(expos_helixID$obesity, remove_first_dummy = TRUE)
meat <- dummy_cols(expos_helixID$meat, remove_first_dummy = TRUE)      
beverages <- dummy_cols(expos_helixID$beverages, remove_first_dummy = TRUE)
dairy <- dummy_cols(expos_helixID$dairy, remove_first_dummy = TRUE)        
fish <- dummy_cols(expos_helixID$fish, remove_first_dummy = TRUE)

colnames(obesity) <- gsub(".","obesity_",colnames(obesity), fixed=TRUE)
colnames(meat) <- gsub(".","meat_",colnames(meat), fixed=TRUE)
colnames(beverages) <- gsub(".","beverages_",colnames(beverages), fixed=TRUE)
colnames(dairy) <- gsub(".","dairy_",colnames(dairy), fixed=TRUE)
colnames(fish) <- gsub(".","fish_",colnames(fish), fixed=TRUE)

expos_helixID_dummies <- cbind(helixID  = expos_helixID$helixID, obesity[,-1],
                               meat[,-1], beverages[,-1],dairy[,-1],fish[,-1])

# Merging dietary intake frequency dummy variables with genotype data and confounding variables
met_exp <- merge(phenotypes_file_urine2, expos_helixID_dummies,
                 by.x = "row.names", by.y = "helixID")
met_exp_gen <- merge(met_exp, gen, by.x = "Row.names", by.y = "IID")


# Example of computation of one GxE model, in this case interaction between meat intake and SNP rs12576623 for taurine
taurine_meat <- lm(Taurine.urine~ rs12576623_C*meat_data_3+
                     rs12576623_C*meat_data_2+sex.data_male+age_sample_years+
                     sample_type.data_CUX+sample_type.data_MUX+PC1+PC2+PC3+
                     PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+
                     PC16+PC17+PC18+PC19+PC20, data = met_exp_gen)

# Null model
taurine_meat_null <- lm(Taurine.urine~ rs12576623_C + meat_data_3+
                          meat_data_2+sex.data_male+age_sample_years+
                          sample_type.data_CUX+sample_type.data_MUX+PC1+PC2+
                          PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+
                          PC14+PC15+PC16+PC17+PC18+PC19+PC20,
                        data = met_exp_gen)

# Comparison of models with ANOVA
anova(taurine_meat_null,taurine_meat)
