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

calculate_TPR_Bonder<- function(CNAsDNA1,CNAsDNA2,CNAsDNA3, PredCNAs1, PredCNAs2, PredCNAs3, NamePredCNAs1, NamePredCNAs2, NamePredCNAs3,colors){
  # Prepare CNVs from DNA
  CNAsDNA1_Del <- CNAsDNA1
  CNAsDNA1_Del[CNAsDNA1_Del >= 0] <- 0
  CNAsDNA1_Del[CNAsDNA1_Del < 0] <- 1
  CNAsDNA1_Amp <- CNAsDNA1
  CNAsDNA1_Amp[CNAsDNA1_Amp <=0] <- 0
  CNAsDNA1_Amp[CNAsDNA1_Amp > 0] <- 1
  CNAsDNA2_Del <- CNAsDNA2
  CNAsDNA2_Del[CNAsDNA2_Del >= 0] <- 0
  CNAsDNA2_Del[CNAsDNA2_Del < 0] <- 1
  CNAsDNA2_Amp <- CNAsDNA2
  CNAsDNA2_Amp[CNAsDNA2_Amp <= 0] <- 0
  CNAsDNA2_Amp[CNAsDNA2_Amp > 0] <- 1
  CNAsDNA3_Del <- CNAsDNA3
  CNAsDNA3_Del[CNAsDNA3_Del >= 0] <- 0
  CNAsDNA3_Del[CNAsDNA3_Del < 0] <- 1
  CNAsDNA3_Amp <- CNAsDNA3
  CNAsDNA3_Amp[CNAsDNA3_Amp <= 0] <- 0
  CNAsDNA3_Amp[CNAsDNA3_Amp > 0] <- 1
  
  # Intersecting Peptides and Cells
  PredCNAs1 <- PredCNAs1[intersect(row.names(PredCNAs1), row.names(CNAsDNA1)),
                         intersect(colnames(PredCNAs1), colnames(CNAsDNA1))]
  cat(nrow(PredCNAs1), "genes and", ncol(PredCNAs1), "cells remaining in PredCNAs1.\n")
  PredCNAs2 <- PredCNAs2[intersect(row.names(PredCNAs2), row.names(CNAsDNA2)),
                         intersect(colnames(PredCNAs2), colnames(CNAsDNA2))]
  cat(nrow(PredCNAs2), "genes and", ncol(PredCNAs2), "cells remaining in PredCNAs2.\n")
  
  PredCNAs3 <- PredCNAs3[intersect(row.names(PredCNAs3), row.names(CNAsDNA3)),
                         intersect(colnames(PredCNAs3), colnames(CNAsDNA3))]
  cat(nrow(PredCNAs3), "genes and", ncol(PredCNAs3), "cells remaining in PredCNAs3.\n")
  
  # Getting data for plotting
  PredDel1 <- EstimateTPR(CNAsDNA1_Del, PredCNAs1, "Del")
  PredAmp1 <- EstimateTPR(CNAsDNA1_Amp, PredCNAs1, "Amp")
  PredDel2 <- EstimateTPR(CNAsDNA2_Del, PredCNAs2, "Del")
  PredAmp2 <- EstimateTPR(CNAsDNA2_Amp, PredCNAs2, "Amp")
  PredDel3 <- EstimateTPR(CNAsDNA3_Del, PredCNAs3, "Del")
  PredAmp3 <- EstimateTPR(CNAsDNA3_Amp, PredCNAs3, "Amp")
 
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

