##############################################################################
################### R script to process metabolomic data #####################
##############################################################################
#
# 19/03/2019
# Beatriz Calvo


# Loading packages
library(minfi)

# Loading the serum and urine datasets
load("/home/isglobal.lan/bcalvo/data/WS_HELIX/HELIX_preproc/metabolome/Final_data/metab_urine_subcohort_v3.RData")
load("/home/isglobal.lan/bcalvo/data/WS_HELIX/HELIX_preproc/metabolome/Final_data/metab_serum_subcohort_v3.RData")

# Extracting metabolites levels
exprs_urine <- exprs(metab_urine_subcohort)
exprs_serum <- exprs(metab_serum_subcohort) 


# Obtaining the transposed to have the samples ID as the row names
exprs_serum_t <- t(exprs_serum)
exprs_urine_t <- t(exprs_urine)

dim(exprs_serum_t) # 1208 samples and  177 metabolites
dim(exprs_urine_t) # 1211 samples and  44 metabolites


##############################################################################
# ---- Merging the two metabolite levels datasets
exprs <- merge(exprs_serum_t, exprs_urine_t, by = "row.names")
rownames(exprs) <- exprs$Row.names
exprs <- exprs[,-1]
dim(exprs) # 1202 samples and 221 metabolites


##############################################################################
# ---- Merging phenodata datasets of the serum and urine ExpressionSet objects
pdata_serum <- pData(metab_serum_subcohort)
pdata_urine <- pData(metab_urine_subcohort)

# Merging the phenodata of serum and urine based on the sample IDs
pdata <- merge(pdata_serum, pdata_urine, by = "row.names")

# The HelixIDs of the urine and serum samples match
identical(pdata$HelixID.ser, pdata$HelixID.uri)                                                                                                                                

# Counting samples with NAs in the HelixID column and removing them
sum(is.na(pdata$HelixID.ser)) # There are 2 NAs
pdata <- pdata[!is.na(pdata$HelixID.ser),]
dim(pdata)

# Deleting the "row.names" column created and assigning the HelixIDs as the 
# actual row names
rownames(pdata) <- pdata$HelixID.ser
pdata <- pdata[,-1]

# Deleting duplicated columns
pdata_final <- pdata[,!duplicated(t(pdata))]

##############################################################################
# ---- Merging the metabolomic data with the phenodata of individuals
exprs_helixID <- merge(exprs, pdata_final, by = "row.names")

# Setting helixIDs as rownames
row.names(exprs_helixID) <- exprs_helixID$HelixID.ser 

# Removing the SampleID and HelixID columns
exprs_helixID <- exprs_helixID[,c(-1,-223,-224)] 

write.csv(exprs_helixID, "expression_serum_urine_helixID.csv")
save(exprs_helixID, file = "expression_serum_urine_helixID.RData")


##############################################################################
# Obtaining reference names for metabolites
p.urine <- pData(featureData(metab_urine_subcohort))

p.serum <- pData(featureData(metab_serum_subcohort))

# Common metabolites in serum and urine
common_chebi <- intersect(p.serum$CHEBI,p.urine$CHEBI) 

common_serum <- p.serum[p.serum$CHEBI%in%common_chebi,]
common_urine <- p.urine[p.urine$CHEBI%in%common_chebi,]

#The common metabolites in serum and urine with equal CHEBI number but
# different reference nomenclature are mostly aminoacids, which are 
# abbreviated (e.g. Lys) in serum but not in urine (e.g. Lysine)


##############################################################################
# --- SERUM
identical(p.serum$Rvar.log, names(exprs_helixID)[1:177]) # TRUE

# Same order in annotation file than in expression file, so we can assign
# the reference metabolite names by:
names(exprs_helixID)[1:177] <- p.serum$Rvar

# Position of the common metabolites according to CHEBI ID in serum
common_serum_pos <- which(names(exprs_helixID)[1:177]%in%common_serum$Rvar) 

identical(names(exprs_helixID)[common_serum_pos],common_serum$Rvar) # TRUE

## NOTE: The metabolites names that are going to be used are the ones in urine  
## which are not abbreviated

# Ordering the urine names following the order of appearance in serum
common_serum_ordered <- common_urine[match(common_serum$CHEBI, 
                                           common_urine$CHEBI),]  

# Assigning the metabolite names to the common metabolites
names(exprs_helixID)[common_serum_pos] <- common_serum_ordered$Rvar

#Adding ".serum" to all serum metabolites
names(exprs_helixID)[1:177] <- paste(names(exprs_helixID)[1:177], ".serum", 
                                     sep = "")


##############################################################################
# ---- URINE

identical(names(exprs_helixID)[178:221], p.urine$Rvar)
# TRUE: urine metabolites already in the reference nomenclature

# Only the ".urine" suffixe has to be added:

names(exprs_helixID)[178:221] <- paste(names(exprs_helixID)[178:221],
                                       ".urine", sep = "")

#Saving file
save(exprs_helixID, file = "expression_helixID_metabolites_names_OK.RData")


##############################################################################
# ---- Computing urine metabolite ratios
c <- 1
p <- 1
name <- NULL
df <- matrix(nrow = 996,ncol = 946)
for (i in 178:220) {
  for (j in (c+178):221) {
    ratio <- exprs_helixID[,i]-exprs_helixID[,j]
    df[,p] <- ratio
    name[p] <- paste(colnames(exprs_helixID)[i], 
                     colnames(exprs_helixID)[j], 
                     sep ="_")
    p <- p + 1
  }
  c <- c + 1
}
colnames(df) <- name

exprs_helixID_ratios <- data.frame(exprs_helixID, df)
write.table(exprs_helixID_ratios, "phenotypes_urine_PCs_plus_ratios.txt",
            sep = "\t", quote =FALSE, row.names = FALSE)
