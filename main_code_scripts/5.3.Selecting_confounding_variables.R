##############################################################################
################### Script to select confounding variables ###################
##############################################################################
#
# 13/04/2019
# Beatriz Calvo

# ---- Testing the existance of correlation between batch variables and 
#      metabolite levels

load("dummy_and_phenotypes_file_workspace.RData")
cov <- pdata_final_metadata_EUR

# Correlation between urine metabolite levels and Run Order variable (batch)
urine_RunOrder <- NULL
urine_RunOrder <- as.list(urine_RunOrder)
for (i in 178:221) {
  res <- lm(expres_helixID[,i]~ cov$RunOrder.uri)
  urine_RunOrder[[i]] <- data.frame("Pval" = summary(res)$coefficients[2,4],
                                    "adj.r.squared" = summary(res)$adj.r.squared)
}
all.data <- do.call(rbind, urine_RunOrder)
cor_urine_RunOrder <- cbind(all.data, 
                            "P.val.adj" = p.adjust(all.data$Pval, 
                                                   method = "BH"))
write.table(cor_urine_RunOrder, file = "correlation_urine_RunOrder")
sum(cor_urine_RunOrder$P.val.adj < 0.05)
length(cor_urine_RunOrder$P.val.adj)


# Correlation between time sampling and cohort
hourvsCohort <- lm(cov$blood_sam4.ser~ cov$cohort)

# PCA of metabolites coloured by cohort
library("RColorBrewer")
pca <- prcomp(expres_helixID, scale = TRUE, center = TRUE)
color <- as.character(factor(cov$cohort, 
                             levels = c("BIB", "EDEN", "KANC", "MOBA",
                                        "RHEA","SAB"), 
                             labels = brewer.pal(n = 6, name = "RdBu")))
xl <- paste("PC1 (", round((pca$sdev[1]/sum(pca$sdev))*100, 2), "%)", sep="")
yl <- paste("PC2 (", round((pca$sdev[2]/sum(pca$sdev))*100, 2), "%)", sep="")
pdf("PCA_metabolite_cohorts.pdf")
plot(pca$x[, 1], pca$x[, 2], xlab=xl, ylab=yl, 
     main="PCA metabolite levels (serum + urine)", 
     cex.axis=0.8, cex.lab=0.8, col = color ,pch=19)
legend("bottomright", legend = c("BIB", "EDEN", "KANC", "MOBA", "RHEA","SAB"),
       col= brewer.pal(n = 6, name = "RdBu"), pch = 19)
dev.off()


# PCA of urine metabolites coloured by SamplingType (batch)
pca2 <- prcomp(expres_helixID[,178:221] , scale = TRUE, center = TRUE)
xl <- paste("PC1 (", round((pca2$sdev[1]/sum(pca2$sdev))*100,2), "%)", sep="")
yl <- paste("PC2 (", round((pca2$sdev[2]/sum(pca2$sdev))*100,2), "%)", sep="")
color <- as.character(factor(cov$Urine.SamplingType.uri, 
                             levels = c("CUX", "MUX","NUX"), 
                             labels = brewer.pal(n = 3, name = "Set2")))
pdf("PCA_metabolite_urine_SamplingType.pdf")
plot(pca2$x[, 1], pca2$x[, 2], xlab=xl, ylab=yl, 
     main="PCA urine metabolite levels by Sampling type", 
     cex.axis=0.8, cex.lab=0.8, col = color, pch = 19)
legend("bottomright", legend = c("CUX", "MUX","NUX"), 
       fill= brewer.pal(n = 3, name = "Set2"))

res.aov.SamplingType <- aov(pca2$x[,1]~cov$Urine.SamplingType.uri)
summary(res.aov.SamplingType)

res.aov.SamplingType2 <- aov(pca2$x[,2]~as.factor(cov$Urine.SamplingType.uri))
summary(res.aov.SamplingType2)

# ---- Testing the existance of correlation between cohort and PCs and 
#      metabolite levels

# Association between urine levels and cohort
p_val_all_coh <- NULL
for (m in 178:221) {
  p_val_coh <- summary(aov(exprs_helixID[,m]~cov$cohort.ser))[[1]][["Pr(>F)"]][1]
  p_val_all_coh <- c(p_val_all_coh,p_val_coh)
}

p_val_all_coh <- p.adjust(p_val_all_coh, method = "BH")
sum(p_val_all_coh < 0.05) 
# 39 of 44 metabolites are statistically significantly associated with cohort

# Association between urine levels and 20 first PCs
p_val_all_PC <- NULL
for (n in 178:221) {
  p_val_PC <- summary(lm(exprs_helixID[,n]~cov$PC1+cov$PC2+cov$PC3+cov$PC4+
                           cov$PC5+cov$PC6+cov$PC7+cov$PC8+cov$PC9+cov$PC10+
                           cov$PC11+cov$PC12+cov$PC13+cov$PC14+cov$PC15+
                           cov$PC16+cov$PC17+cov$PC18+cov$PC19+
                           cov$PC20))$coefficients[,4]
  p_val_all_PC <- rbind(p_val_all_PC,p_val_PC)
}

p_val_all_PC_adj <- apply(p_val_all_PC, 2,
                          function(x) p.adjust(x, method = "BH"))
sum(p_val_all_PC < 0.05) 
# 87 of 880 metabolite-PC associations are statistically significant


# Association between cohort variable and 20 first PCs

PCAvsCohort_res <- matrix(nrow = 20, ncol = 2)
colnames(PCAvsCohort_res) <- c("Adj_r2","P-value")
for(i in 1:20) {
  PCAvsCohort <- lm(cov[,i+6]~ cov$cohort.ser)
  PCAvsCohort_res[i,] <- c(summary(PCAvsCohort)$adj.r.squared,
                           summary(PCAvsCohort)$coefficients[2,4])
}