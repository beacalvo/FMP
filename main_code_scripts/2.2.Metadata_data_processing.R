##############################################################################
# ---- Selecting only individuals with European predicted ancestry
exprs_helixID_ratios_EUR <- exprs_helixID_ratios[exprs_helixID_ratios$
                                                   FINAL_ancestry=="EUR",]
dim(exprs_helixID_ratios_EUR)
