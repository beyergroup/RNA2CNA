source.dir <- "../../../../"
source(paste0(source.dir,"packages.R"))
source(paste0(source.dir,"functions.R"))

###### ROC curves
source(paste0(source.dir,'Functions/ROC_3CNAanalysis_function.R'))
cellnet.dir <- "../../../../../../../../../cellnet/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_CBS/"


input_ROC <- function(CNAs_Matrix, batch){
  # DNA
  DNA_gene_matrix <- readRDS(file=paste0("",path,"/DNA_gene_matrix_",batch,"_woX.rds"))
  DNA_gene_matrix <- na.omit(DNA_gene_matrix)
  metadata_DNA <- readRDS(file=paste0("",path,"/metadata_DNA_",batch,".rds"))
  
  SegStats <- read.delim(file=paste0("/cellfile/cellnet/Ronja/Tools/ginkgo-master/uploads/global/bin_10Mb/Bonder_",batch,"/SegStats"))
  summary(SegStats$Reads)
  
  keep <- rownames(SegStats[which(SegStats$Reads >= threshold),])
  DNA_gene_matrix <- DNA_gene_matrix[,keep]
  
  ########## 
  CNAs_Matrix
  
  #colnames(DNA_gene_matrix)
  #colnames(CNAs_Matrix)
  #head(metadata_DNA)
  # intersect DNA and RNA cell names
  #colnames(CNAs_Matrix)
  DNA_cells<- intersect(rownames(metadata_DNA),colnames(DNA_gene_matrix))
  RNA_cells<- intersect(metadata_DNA$Cell_ID,colnames(CNAs_Matrix))
  
  metadata_DNA <- na.omit(metadata_DNA[DNA_cells,])
  rownames(metadata_DNA) <- metadata_DNA$Cell_ID
  metadata_DNA <- na.omit(metadata_DNA[RNA_cells,])
  
  CNAs_Matrix <- CNAs_Matrix[,metadata_DNA$Cell_ID]
  DNA_gene_matrix <- DNA_gene_matrix[,colnames(DNA_gene_matrix)]
  
  all.equal(rownames(metadata_DNA),colnames(CNAs_Matrix))
  colnames(CNAs_Matrix) <- metadata_DNA$DNA_name
  DNA_gene_matrix <- DNA_gene_matrix[,colnames(CNAs_Matrix)]
  dim(CNAs_Matrix) # batch1 338 cells ; batch2 129 cells
  dim(DNA_gene_matrix) 
  
  # center the DNA data
  DNA_gene_matrix_centered <- apply(DNA_gene_matrix, 2, function(cell) cell-median(cell,na.rm = T))
  
  
  return(list(DNA=DNA_gene_matrix_centered,RNA=CNAs_Matrix))
}

########## BONDER #############
library(ggpubr)

# source scripts
files.sources = list.files("/cellfile/datapublic/rjohnen/Ronja/Tumor_Evolution/Analysis/CNAcalling/R/")
sapply(files.sources, function(x) source(paste0("/cellfile/datapublic/rjohnen/Ronja/Tumor_Evolution/Analysis/CNAcalling/R/",x)))
library(colorRamp2)
source(paste0(source.dir,'Functions/ROC_2CNAanalysis_function_diff_datasets_GINKGO.R'))
source(paste0(source.dir,'Functions/ROC_3CNAanalysis_function_diff_datasets_GINKGO.R'))
source(paste0(source.dir,'Functions/ROC_2CNAanalysis_function_diff_datasets_GINKGO_PVAL.R'))
# source(paste0(source.dir,'Functions/Create_Heatmap_pval.R'))
source(paste0(source.dir,'Functions/CreatePVALMatrix_fromLargeDNAcopy.R'))
source(paste0(source.dir,'Functions/Imputation_CNA_PVAL.R'))

######## gencode##########
gencode <- readRDS(file='/cellfile/datapublic/rjohnen/Ronja/Tumor_Evolution/Analysis/Datasets/Gencode/Mouse/gencode.vM23.annotation.rds')
gencode_gene_types <- readRDS(file='/cellfile/datapublic/rjohnen/Ronja/Tumor_Evolution/Analysis/Datasets/Gencode/Mouse/gencode_gene_types.rds')
###### RNA2CNA data #######
CBS_pvalues <- readRDS(file=paste0(cellnet.dir,'/Blood/Bonder/output/scaledata/genewise/FixedMovingWindow/CBS_pvalues.rds'))
CNAs_per_gene <- readRDS(file=paste0(cellnet.dir,'/Blood/Bonder/output/scaledata/genewise/FixedMovingWindow/CNAs_per_gene.rds'))
data <- readRDS(file=paste0(cellnet.dir,'/Blood/Bonder/output/scaledata/genewise/FixedMovingWindow/data.rds'))
dim(data)
metadata <- data@meta.data

CBS_pvalues[which(is.na(CBS_pvalues$num.mark)),"num.mark"] <- 0
CNA_Matrix<-CreateMatrix_fromLargeDNAcopy(DNAcopy = list(output=CBS_pvalues), CNAProf = CNAs_per_gene,gencode = gencode)

CBS_Pval<-CBS_pvalues
CBS_Pval$seg.mean <- CBS_Pval$FDR
Pval_M<-CreatePVALMatrix_fromLargeDNAcopy(DNAcopy = list(output=CBS_Pval), CNAProf = CNAs_per_gene,gencode = gencode)

colnames(CBS_pvalues)[2]<-"Chr"
colnames(CBS_pvalues)[6]<-"CNA_value"
colnames(CBS_Pval)[2]<-"Chr"
colnames(CBS_Pval)[6]<-"CNA_value"

#rm CNAs with start smaller than end location 
CBS_pvalues<-CBS_pvalues[CBS_pvalues$loc.start < CBS_pvalues$loc.end,]
CBS_Pval<-CBS_Pval[CBS_Pval$loc.start < CBS_Pval$loc.end,]

CNA_Matrix_Imputed <- Imputation_CNA(CNA_Matrix,CBS_pvalues,gencode_gene_types) # no threshold 
Pval_M_Imputed <- Imputation_CNA_PVAL(Pval_M,CBS_Pval,gencode_gene_types) # no threshold, but use the p-values as classifier

path <- "../Validation_Healthy/output/global/bin_10Mb"# #"independent/bin_10Mb/read150/"
threshold<-2e+5

RNA2CNA_batch1_pval<-input_ROC(Pval_M_Imputed, batch='batch1')
RNA2CNA_batch2_pval<-input_ROC(Pval_M_Imputed, batch='batch2')
RNA2CNA_batch1_not_imp<-input_ROC(CNA_Matrix, batch='batch1')
RNA2CNA_batch2_not_imp<-input_ROC(CNA_Matrix, batch='batch2')
RNA2CNA_batch1_raw<-input_ROC(CNAs_per_gene, batch='batch1')
RNA2CNA_batch2_raw<-input_ROC(CNAs_per_gene, batch='batch2')

#### RNA data #### 

#copykat
Bonder_copykat <- read.delim("/cellfile/datapublic/rjohnen/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/copyKat/Bonder//Bonder_copykat_CNA_raw_results_gene_by_cell.txt")
rownames(Bonder_copykat) <- Bonder_copykat$mgi_symbol
Bonder_copykat <- Bonder_copykat[,-c(1:7)]
# match cell names
copykat_batch1<-input_ROC(Bonder_copykat, batch='batch1')
copykat_batch2<-input_ROC(Bonder_copykat, batch='batch2')

#inferCNV
infercnv.preliminary.observations <- read.csv("/cellfile/cellnet/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/inferCNV/Bonder/output/infercnv.preliminary.observations.txt", sep="")

# match cell names
infercnv_batch1<-input_ROC(infercnv.preliminary.observations, batch='batch1')
infercnv_batch2<-input_ROC(infercnv.preliminary.observations, batch='batch2')


saveRDS(  cbind(RNA2CNA_batch1_raw$DNA,RNA2CNA_batch2_raw$DNA), "output_manuscript/Bonder_DNA_RNA2CNA.rds")
saveRDS( cbind(infercnv_batch1$DNA,infercnv_batch2$DNA), "output_manuscript/Bonder_DNA_infercnv.rds")
saveRDS( cbind(copykat_batch1$DNA,copykat_batch2$DNA), "output_manuscript/Bonder_DNA_copykat.rds")

saveRDS(  cbind(RNA2CNA_batch1_raw$RNA,RNA2CNA_batch2_raw$RNA), "output_manuscript/Bonder_RNA_RNA2CNA.rds")
saveRDS( cbind(infercnv_batch1$RNA,infercnv_batch2$RNA), "output_manuscript/Bonder_RNA_infercnv.rds")
saveRDS( cbind(copykat_batch1$RNA,copykat_batch2$RNA), "output_manuscript/Bonder_RNA_copykat.rds")

###
source(paste0(source.dir,'Functions/ROC_3CNAanalysis_function_diff_datasets_GINKGO.R'))
pdf(paste0("plots_manuscript/Bonder_inferCNV_copyKat_RNA2CNA_imputed_QC_",threshold,"_centered.pdf"),width=10)
ROC_Del_Amp_3Analysis(CNAsDNA1 = cbind(RNA2CNA_batch1_not_imp$DNA,RNA2CNA_batch2_not_imp$DNA),
                      CNAsDNA2 = cbind(copykat_batch1$DNA,copykat_batch2$DNA),
                      CNAsDNA3 = cbind(infercnv_batch1$DNA,infercnv_batch2$DNA),
                      PredCNAs1 = cbind(RNA2CNA_batch1_not_imp$RNA,RNA2CNA_batch2_not_imp$RNA),
                      PredCNAs2 = cbind(copykat_batch1$RNA,copykat_batch2$RNA),
                      PredCNAs3 = cbind(infercnv_batch1$RNA,infercnv_batch2$RNA),
                      NamePredCNAs1 = "RNA2CNA", NamePredCNAs2 = "CopyKat",NamePredCNAs3 = "InferCNV",
                      colors=c("#BB5566","#004488","#DDAA33"))
dev.off()
pdf(paste0("plots_manuscript/Bonder_inferCNV_copyKat_RNA2CNA_imputed_QC_",threshold,"_centered2.pdf"),width=10)
ROC_Del_Amp_3Analysis(CNAsDNA1 = cbind(RNA2CNA_batch1_raw$DNA,RNA2CNA_batch2_raw$DNA),
                      CNAsDNA2 = cbind(infercnv_batch1$DNA,infercnv_batch2$DNA),
                      CNAsDNA3 = cbind(copykat_batch1$DNA,copykat_batch2$DNA),
                      PredCNAs1 = cbind(RNA2CNA_batch1_raw$RNA,RNA2CNA_batch2_raw$RNA),
                      PredCNAs2 = cbind(infercnv_batch1$RNA,infercnv_batch2$RNA),
                      PredCNAs3 = cbind(copykat_batch1$RNA,copykat_batch2$RNA),
                      NamePredCNAs1 = "RNA2CNA", NamePredCNAs2 = "InferCNV",NamePredCNAs3 = "CopyKat",
                      colors=c("#BB5566","#004488","#DDAA33"))
dev.off()

source("Precred_Del_Amp.R")
library(precrec)
p<-Precrec_Del_Amp(CNAsDNA1 = cbind(RNA2CNA_batch1_raw$DNA,RNA2CNA_batch2_raw$DNA),
                   CNAsDNA2 = cbind(infercnv_batch1$DNA,infercnv_batch2$DNA),
                   CNAsDNA3 = cbind(copykat_batch1$DNA,copykat_batch2$DNA),
                scRNA_df1 = cbind(RNA2CNA_batch1_raw$RNA,RNA2CNA_batch2_raw$RNA),
                scRNA_df2 = cbind(infercnv_batch1$RNA,infercnv_batch2$RNA),
                scRNA_df3 = cbind(copykat_batch1$RNA,copykat_batch2$RNA),
                NamePred1 = "RNA2CNA", NamePred2 = "InferCNV",NamePred3 = "CopyKat",
                colors=c("#BB5566","#004488","#DDAA33"))
png(paste0("plots_manuscript/PRECREC_Bonder_inferCNV_copyKat_RNA2CNA_imputed_QC_",threshold,"_centered.png"),width=1000)
print(p[[1]])
dev.off()
