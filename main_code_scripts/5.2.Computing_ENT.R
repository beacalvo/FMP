##############################################################################
######################### Script to compute ENT ##############################
##############################################################################
#
# 12/05/2019
# Beatriz Calvo


##############################################################################
# ---- Computing ENT

urine <- read.table("phenotypes_urine.txt", header=T)

TEF<-function(data){
  library(lsr)
  
  a1<-99
  
  ### Correlation Matrix 
  for (i in 1:dim(data)[2]) {
    ex_i<-colnames(data[i])   
    ty_i <- class(data[,i])
    for (j in 1:dim(data)[2]) {
      ex_j<-colnames(data[j])    
      ty_j <- class(data[,j])
      if(ty_i == "numeric" & ty_j == "numeric") {
        a<-cor(x = data[ , c(ex_i, ex_j)])[1, 2]
      } else if(ty_i == "factor" & ty_j == "factor") {
        vi2<-as.numeric(data[,i])
        vj2<-as.numeric(data[,j])
        if (cor(vi2,vj2)>0){
          a<-cramersV(x = table(data[ , c(ex_i, ex_j)]))
        }
        else  if (cor(vi2,vj2)<0){
          a<-(-cramersV(x = table(data[ , c(ex_i, ex_j)])))
        }  
      } else {
        vi2<-as.numeric(data[,i])
        vj2<-as.numeric(data[,j])
        fm <- paste(ex_i, "~", ex_j)
        if(ty_i == "factor") {
          fm <- paste(ex_j, "~", ex_i)
        }
        fm <- lm(formula = fm, data = data[ , c(ex_i, ex_j)])
        if  (cor(vi2,vj2)>0){
          a<-sqrt(summary(fm)$r.squared)
        } else  if  (cor(vi2,vj2)<0){
          a<-(-sqrt(summary(fm)$r.squared))
        }
      }
      a1<-cbind(a1,a)
    }
    #a1<-cbind(a1,a)
  }
  
  a2<-matrix(a1[-1],nrow=dim(data)[2],ncol=dim(data)[2])
  
  colnames(a2) <- colnames(data)
  rownames(a2) <- colnames(data)
  
  # TEF outside rexposome
  M <- ncol(a2)
  lambdas <- eigen(a2)$values
  Vobs <- sum(((lambdas - 1)^2)) / (M - 1)
  Meff <- M - sum((lambdas>1)*(lambdas-1))
  alpha_corrected <- 1 - (1 - 0.05)^(1 / Meff)
  return(c(Meff,alpha_corrected))
}


res <- TEF(urine[,203:246])

5e-8/res[1]
# [1] 1.522778e-09

