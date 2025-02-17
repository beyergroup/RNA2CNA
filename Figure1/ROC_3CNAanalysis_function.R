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

ROC_Del_Amp_3Analysis <- function(CNAsDNA, PredCNAs1, PredCNAs2, PredCNAs3, NamePredCNAs1, NamePredCNAs2, NamePredCNAs3,colors){
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
  PredDel1 <- EstimateROC(CNAsDNA_Del, PredCNAs1, "Del")
  PredAmp1 <- EstimateROC(CNAsDNA_Amp, PredCNAs1, "Amp")
  PredDel2 <- EstimateROC(CNAsDNA_Del, PredCNAs2, "Del")
  PredAmp2 <- EstimateROC(CNAsDNA_Amp, PredCNAs2, "Amp")
  PredDel3 <- EstimateROC(CNAsDNA_Del, PredCNAs3, "Del")
  PredAmp3 <- EstimateROC(CNAsDNA_Amp, PredCNAs3, "Amp")
  cat("Starting with plotting.\n")
  DelPlot <- GetPlot(PredDel1, PredDel2, PredDel3, "Deletions", NamePredCNAs1, NamePredCNAs2,NamePredCNAs3,colors)
  AmpPlot <- GetPlot(PredAmp1, PredAmp2, PredAmp3, "Amplifications", NamePredCNAs1, NamePredCNAs2,NamePredCNAs3,colors)
  legend <- GetLegend(NamePredCNAs1, NamePredCNAs2,NamePredCNAs3,colors)
  cat("Plotting...\n")

  return(list(AmpPlot,DelPlot, legend))
  cat("Done!\n")
}

EstimateROC <- function(CNAs, PredCNAs, type){
  CNAs <- as.vector(as.matrix(CNAs[row.names(PredCNAs), colnames(PredCNAs)]))
  if(type == "Del"){
    PredCNAs <- as.vector(as.matrix(PredCNAs *(-1)))
  }
  if(type == "Amp"){
    PredCNAs <- as.vector(as.matrix(PredCNAs))
  }
  df <- data.frame(ref = CNAs,
                   Pred = PredCNAs)
  p <- ggplot(df, aes(d = ref, m = Pred)) +
    geom_roc()
  auc <- calc_auc(p)
  return(list(CNAs, PredCNAs, auc))
}

GetPlot <- function(Pred1, Pred2, Pred3, title, NamePred1, NamePred2, NamePred3,colors){
  df <- data.frame(ref = c(Pred1[[1]],
                           Pred2[[1]],
                           Pred3[[1]]),
                   Pred = c(Pred1[[2]],
                            Pred2[[2]],
                            Pred3[[2]]),
                   Method = c(rep(NamePred1, length(Pred1[[1]])),
                              rep(NamePred2, length(Pred2[[1]])),
                              rep(NamePred3, length(Pred3[[1]]))))
  df$Method <- factor(df$Method, levels = c(NamePred1,NamePred2,NamePred3))

  Plot <- ggplot(df, aes(d = ref, m = Pred, col = Method)) +
    geom_roc(labels = FALSE, n.cuts = 0)+
    scale_color_manual(breaks=c(NamePred1, NamePred2, NamePred3),
                       values = colors)+
    labs(title = title, x = "FPR", y = "TPR") +
    annotate("text", x = 0.5, y = 0.25,
             label = paste0("AUC\n", NamePred1, " = ", round(Pred1[[3]][[3]],3),
                            "\n", NamePred2," = ", round(Pred2[[3]][[3]],3),
                            "\n", NamePred3," = ", round(Pred3[[3]][[3]],3)),
             size = 1.5, hjust = 0) +
    theme(plot.title = element_text(size = 8, face = "bold"),
          legend.position ="none",
          axis.title = element_text(size = 8), axis.text = element_text(size = 8),aspect.ratio=1) +
    geom_abline(intercept = 0, slope = 1,
                color = "darkgrey", linetype = "dashed") #+ style_roc()
  return(Plot)
}

GetLegend <- function(NamePred1, NamePred2, NamePred3,colors){
  df <- data.frame(x = c(1,2,3,4,5,6,7,8,9),
                   y = c(1,2,3,4,5,6,7,8,9),
                   Method = c(rep(NamePred1, 3),
                              rep(NamePred2, 3),
                              rep(NamePred3, 3)))
  df$Method <- factor(df$Method, levels = c(NamePred1,NamePred2,NamePred3))

  leg <- ggplot(df, aes(x = x, y = y, col = Method)) +
    geom_line() +
    scale_color_manual(breaks = c(NamePred1, NamePred2, NamePred3),
                       values = colors) +
    theme(legend.text = element_text(size = 8), legend.title = element_text(size = 8))
  leg <- as_ggplot(get_legend(leg))
  return(leg)
}
