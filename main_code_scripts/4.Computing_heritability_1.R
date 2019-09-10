##############################################################################
############################## GCTA analysis #################################
##############################################################################

##############################################################################
# ---- Input files to run GCTA analysis

fam <- read.table("HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT.fam", 
                  sep = " ")
head(fam)
id <- fam[,1:2]

load("dummy_and_phenotypes_file_workspace.RData")

colnames(id) <- c("fam_id","ID") 

all <- merge(phenotypes_file_urine_plus_ratios, id, 
             by.x="HelixID.ser", by.y="ID")

all_ordered <- all[match(id$ID, 
                         phenotypes_file_urine_plus_ratios$HelixID.ser),]

match(id$ID, phenotypes_file_urine_plus_ratios$HelixID.ser)[
  !is.na(match(id$ID,phenotypes_file_urine_plus_ratios$HelixID.ser))]

all_ordered <- all[order,]

id_keep <- all_ordered[,c(2,1)] 
ID <- id_keep$ID
fam_id <- id_keep$fam_id

id_keep[,1] <- fam_id
id_keep[,2] <- as.character(ID)

write.table(id_keep, file = "ID_to_keep.txt", sep = "\t", quote=FALSE,
            row.names = FALSE, col.names = FALSE)

save.image("files_gtca.RData")