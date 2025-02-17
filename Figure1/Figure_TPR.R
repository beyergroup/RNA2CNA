source.dir <- "../../../Analysis/"
source(paste0(source.dir,"packages.R"))
source(paste0(source.dir,"functions.R"))
library(colorRamp2)
source(paste0('scripts/calculate_TPR.R'))

Zacharidias_DNA_CNA <- readRDS("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/inferCNV_copyKat_ROC/output_manuscript/Zacharidias_DNA_CNA.rds")
Zacharidias_RNA_CNA_RNA2CNA <- readRDS("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/inferCNV_copyKat_ROC/output_manuscript/Zacharidias_RNA_CNA_RNA2CNA.rds")
Zacharidias_RNA_CNA_infercnv <- readRDS("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/inferCNV_copyKat_ROC/output_manuscript/Zacharidias_RNA_CNA_infercnv.rds")
Zacharidias_RNA_CNA_copykat <- readRDS("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/inferCNV_copyKat_ROC/output_manuscript/Zacharidias_RNA_CNA_copykat.rds")

Bian_DNA_CNA <- readRDS("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/inferCNV_copyKat_ROC/output_manuscript/Bian_DNA_CNA.rds")
Bian_RNA_CNA_RNA2CNA <- readRDS("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/inferCNV_copyKat_ROC/output_manuscript/Bian_RNA_CNA_RNA2CNA.rds")
Bian_RNA_CNA_infercnv <- readRDS("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/inferCNV_copyKat_ROC/output_manuscript/Bian_RNA_CNA_infercnv.rds")
Bian_RNA_CNA_copykat <- readRDS("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/inferCNV_copyKat_ROC/output_manuscript/Bian_RNA_CNA_copykat.rds")

Bonder_DNA_RNA2CNA <- readRDS("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/inferCNV_copyKat_ROC/output_manuscript/Bonder_DNA_RNA2CNA.rds")
Bonder_DNA_infercnv <- readRDS("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/inferCNV_copyKat_ROC/output_manuscript/Bonder_DNA_infercnv.rds")
Bonder_DNA_copykat <- readRDS("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/inferCNV_copyKat_ROC/output_manuscript/Bonder_DNA_copykat.rds")

Bonder_RNA_CNA_RNA2CNA <- readRDS("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/inferCNV_copyKat_ROC/output_manuscript/Bonder_RNA_RNA2CNA.rds")
Bonder_RNA_CNA_infercnv <- readRDS("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/inferCNV_copyKat_ROC/output_manuscript/Bonder_RNA_infercnv.rds")
Bonder_RNA_CNA_copykat <- readRDS("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/inferCNV_copyKat_ROC/output_manuscript/Bonder_RNA_copykat.rds")

### calculate the TPR at FDR 0.1
library(pROC)

Bian <- calulcate_TPR(CNAsDNA = Bian_DNA_CNA, 
              PredCNAs1 = Bian_RNA_CNA_RNA2CNA,
              PredCNAs2 = Bian_RNA_CNA_infercnv,
              PredCNAs3 = Bian_RNA_CNA_copykat,
              NamePredCNAs1 = "RNA2CNA",
              NamePredCNAs2 = "InferCNV",
              NamePredCNAs3 = "CopyKAT",
              colors=c("#BB5566","#004488","#DDAA33"))
Zacharidias <-calulcate_TPR(CNAsDNA = Zacharidias_DNA_CNA, 
                          PredCNAs1 = Zacharidias_RNA_CNA_RNA2CNA,
                          PredCNAs2 = Zacharidias_RNA_CNA_infercnv,
                          PredCNAs3 = Zacharidias_RNA_CNA_copykat,
                          NamePredCNAs1 = "RNA2CNA",
                          NamePredCNAs2 = "InferCNV",
                          NamePredCNAs3 = "CopyKAT",
                          colors=c("#BB5566","#004488","#DDAA33"))
source(paste0('scripts/calculate_TPR_Bonder.R'))
Bonder<-calculate_TPR_Bonder(CNAsDNA1 = Bonder_DNA_RNA2CNA, 
                          CNAsDNA2 = Bonder_DNA_infercnv, 
                          CNAsDNA3 = Bonder_DNA_copykat,
                          PredCNAs1 = Bonder_RNA_CNA_RNA2CNA,
                          PredCNAs2 = Bonder_RNA_CNA_infercnv,
                          PredCNAs3 = Bonder_RNA_CNA_copykat,
                          NamePredCNAs1 = "RNA2CNA",
                          NamePredCNAs2 = "InferCNV",
                          NamePredCNAs3 = "CopyKAT",
                          colors=c("#BB5566","#004488","#DDAA33"))


Bian$Dataset <- "Bian"
Zacharidias$Dataset <- "Zacharidias"
Bonder$Dataset <- "Bonder"
Bonder <- Bonder[!(Bonder$Dataset=="Bonder" & Bonder$CNA=="Deletion"),]

df <- rbind(Bian,Zacharidias,Bonder)

df$Method <- as.factor(df$Method)
levels(df$Method)
df$Method <- factor(df$Method,levels=c("RNA2CNA","InferCNV","CopyKAT"))
levels(df$Method)
df$Dataset <- as.factor(df$Dataset)
levels(df$Dataset)
df$Dataset <- factor(df$Dataset,levels=c("Bian","Zacharidias","Bonder"))
levels(df$Dataset)

# remove deletions from bonder dataset
df <- droplevels(df)
str(df)


summary(Bian_DNA_CNA[,1])
summary(Zacharidias_DNA_CNA[,1])
summary(Bonder_DNA_RNA2CNA[,1])
length(which(Bonder_DNA_RNA2CNA < 0))

df_label <- data.frame("Bian_Amp"=length(which(Bian_DNA_CNA > 2)), "Bian_Del"=length(which(Bian_DNA_CNA < 2)),
                      "Zacharidias_Amp"=length(which(Zacharidias_DNA_CNA > 2)), "Zacharidias_Del"=length(which(Zacharidias_DNA_CNA < 2)),
                      "Bonder_Amp"=length(which(Bonder_DNA_RNA2CNA > 0)))

annotation <- data.frame(
  TPR = c(0.81,0.81,0.81,0.81,0.81),
  Dataset=c("Bian","Bian","Zacharidias","Zacharidias","Bonder"),
  CNA=c("Amplification","Deletion","Amplification","Deletion","Amplification"),
  label = as.numeric(df_label)
)

pdf('plots/Figure1b_TPR.pdf',width=2.5,height =2.5)
ggplot(df)+geom_bar(aes(y=TPR,x=Dataset,fill=Method),stat="identity",position = "dodge")+facet_grid(.~CNA, scales = "free", space = "free") + theme_manuscript() + theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top") + 
  scale_fill_manual(breaks = c("RNA2CNA","InferCNV","CopyKAT"), values=c("#BB5566","#004488","#DDAA33")) + ylim(0,0.85) # geom_text(data=annotation, aes( x=Dataset, y=TPR, label=label), size=3) +
dev.off()

rm();gc();malloc.trim()
