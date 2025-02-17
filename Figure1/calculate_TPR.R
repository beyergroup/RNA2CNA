library(ggplot2)
library(ggpubr)
library(plotROC)

# Function for estimating ROC curves with two CNA calling runs
# Input:	CNAsDNA = dataframe with DNA based CNAs (absolut number of copies, not ratios); rows are genes, cols are cells
#		PredCNAs1 = dataframe of first CNA calling variation; rows are genes, cols are cells
#		PredCNAs2 = dataframe of second CNA calling variation; rows are genes, cols are cells
#		PredNamePredCNAs1 = Name of first CNA calling variation
#		PredNamePredCNAs2 = Name of second CNA calling variation
# Output:	Plot with ROC curves for calling of deletions and amplifications

calulcate_TPR <- function(CNAsDNA, PredCNAs1, PredCNAs2, PredCNAs3, NamePredCNAs1, NamePredCNAs2, NamePredCNAs3,colors){
  # Prepare CNVs from DNA
  CNAsDNA_Del <- CNAsDNA
  CNAsDNA_Del[CNAsDNA_Del < 2] <- 1
  CNAsDNA_Del[CNAsDNA_Del >= 2] <- 0
  CNAsDNA_Amp <- CNAsDNA
  CNAsDNA_Amp[CNAsDNA_Amp <= 2] <- 0
  CNAsDNA_Amp[CNAsDNA_Amp > 2] <- 1
  
  # Intersecting Peptides and Cells
  PredCNAs1 <- PredCNAs1[intersect(row.names(PredCNAs1), row.names(CNAsDNA)),
                         intersect(colnames(PredCNAs1), colnames(CNAsDNA))]
  cat(nrow(PredCNAs1), "genes and", ncol(PredCNAs1), "cells remaining in PredCNAs1.\n")
  PredCNAs2 <- PredCNAs2[intersect(row.names(PredCNAs2), row.names(CNAsDNA)),
                         intersect(colnames(PredCNAs2), colnames(CNAsDNA))]
  cat(nrow(PredCNAs2), "genes and", ncol(PredCNAs2), "cells remaining in PredCNAs2.\n")
  
  PredCNAs3 <- PredCNAs3[intersect(row.names(PredCNAs3), row.names(CNAsDNA)),
                         intersect(colnames(PredCNAs3), colnames(CNAsDNA))]
  cat(nrow(PredCNAs3), "genes and", ncol(PredCNAs3), "cells remaining in PredCNAs3.\n")
  
  # Getting data for plotting
  PredDel1 <- EstimateTPR(CNAsDNA_Del, PredCNAs1, "Del")
  PredAmp1 <- EstimateTPR(CNAsDNA_Amp, PredCNAs1, "Amp")
  PredDel2 <- EstimateTPR(CNAsDNA_Del, PredCNAs2, "Del")
  PredAmp2 <- EstimateTPR(CNAsDNA_Amp, PredCNAs2, "Amp")
  PredDel3 <- EstimateTPR(CNAsDNA_Del, PredCNAs3, "Del")
  PredAmp3 <- EstimateTPR(CNAsDNA_Amp, PredCNAs3, "Amp")

  df<-data.frame("TPR"=c(PredAmp1,PredDel1,PredAmp2,PredDel2,PredAmp3,PredDel3),"Method"=c(NamePredCNAs1,NamePredCNAs1,NamePredCNAs2,NamePredCNAs2,NamePredCNAs3,NamePredCNAs3),"CNA"=rep(c("Amplification","Deletion"),3))
    
  return(df)
  cat("Done!\n")
}

EstimateTPR <- function(CNAs, PredCNAs, type){
  CNAs <- as.vector(as.matrix(CNAs[row.names(PredCNAs), colnames(PredCNAs)]))
  if(type == "Del"){
    PredCNAs <- as.vector(as.matrix(PredCNAs *(-1)))
  }
  if(type == "Amp"){
    PredCNAs <- as.vector(as.matrix(PredCNAs))
  }
  df <- data.frame(ref = CNAs,
                   Pred = PredCNAs)
  myroc <- roc(response=df$ref, predictor=df$Pred)
  target.fpr <- c(0.1)
  target.sp <- 1 - target.fpr
  coords(myroc, x = target.sp, input = "specificity", ret = c("se", "sp"))
  tpr <- coords(myroc, x = target.sp, input = "specificity", ret = "se")[1,]
  return(tpr)
}

