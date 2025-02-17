source.dir <- "../../../Analysis/"
source(paste0(source.dir,"packages.R"))
source(paste0(source.dir,"functions.R"))
library(colorRamp2)
source(paste0('scripts/ROC_3CNAanalysis_function.R'))

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


p1<-ROC_Del_Amp_3Analysis(CNAsDNA = Bian_DNA_CNA, 
                          PredCNAs1 = Bian_RNA_CNA_RNA2CNA,
                          PredCNAs2 = Bian_RNA_CNA_infercnv,
                          PredCNAs3 = Bian_RNA_CNA_copykat,
                          NamePredCNAs1 = "RNA2CNA",
                          NamePredCNAs2 = "InferCNV",
                          NamePredCNAs3 = "CopyKAT",
                          colors=c("#BB5566","#004488","#DDAA33"))


p2<-ROC_Del_Amp_3Analysis(CNAsDNA = Zacharidias_DNA_CNA, 
                          PredCNAs1 = Zacharidias_RNA_CNA_RNA2CNA,
                          PredCNAs2 = Zacharidias_RNA_CNA_infercnv,
                          PredCNAs3 = Zacharidias_RNA_CNA_copykat,
                          NamePredCNAs1 = "RNA2CNA",
                          NamePredCNAs2 = "InferCNV",
                          NamePredCNAs3 = "CopyKAT",
                          colors=c("#BB5566","#004488","#DDAA33"))

source(paste0('scripts//ROC_3CNAanalysis_function_diff_datasets_GINKGO.R'))
p3<-ROC_Del_Amp_3Analysis(CNAsDNA1 = Bonder_DNA_RNA2CNA, 
                          CNAsDNA2 = Bonder_DNA_infercnv, 
                          CNAsDNA3 = Bonder_DNA_copykat,
                          PredCNAs1 = Bonder_RNA_CNA_RNA2CNA,
                          PredCNAs2 = Bonder_RNA_CNA_infercnv,
                          PredCNAs3 = Bonder_RNA_CNA_copykat,
                          NamePredCNAs1 = "RNA2CNA",
                          NamePredCNAs2 = "InferCNV",
                          NamePredCNAs3 = "CopyKAT",
                          colors=c("#BB5566","#004488","#DDAA33"))



pdf('plots/Figure1b_Legende.pdf',width=3.5)
p1[[3]]
dev.off()

p1 <- ggarrange(p1[[1]],p1[[2]],ncol=2, widths = c(2,2))
p2 <- ggarrange(p2[[1]],p2[[2]], ncol=2, widths = c(2,2))
p3 <- ggarrange(p3[[1]],p3[[2]], ncol=2, widths = c(2,2))

pdf('plots/Figure1b_Bian.pdf',width=3.5)
p1
dev.off()
pdf('plots/Figure1b_Zacharidias.pdf',width=3.5)
p2
dev.off()
pdf('plots/Figure1b_Bonder.pdf',width=3.5)
p3
dev.off()
