source.dir <- "../../../../"
source(paste0(source.dir,"packages.R"))
source(paste0(source.dir,"functions.R"))
# source scripts
files.sources = list.files("/cellfile/datapublic/rjohnen/Ronja/Tumor_Evolution/Analysis/CNAcalling/R/")
sapply(files.sources, function(x) source(paste0("/cellfile/datapublic/rjohnen/Ronja/Tumor_Evolution/Analysis/CNAcalling/R/",x)))
library(colorRamp2)

source(paste0(source.dir,'Functions/ROC_2CNAanalysis_function_diff_datasets_PVAL.R'))
source(paste0(source.dir,'Functions/CreatePVALMatrix_fromLargeDNAcopy.R'))
source(paste0(source.dir,'Functions/Imputation_CNA_PVAL.R'))

gencode <- readRDS(file='/cellfile/datapublic/rjohnen/Ronja/Tumor_Evolution/Analysis/Datasets/Gencode/gencode_v28.rds')
gencode_gene_types <- readRDS(file='/cellfile/datapublic/rjohnen/Ronja/Tumor_Evolution/Analysis/Datasets/Gencode/gencode_gene_types.rds')

###### ROC curves
source(paste0(source.dir,'Functions/ROC_3CNAanalysis_function.R'))
cellnet.dir <- "../../../../../../../../../cellnet/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_CBS/"

########## Bian #############
## Reference data
mapped_CNA_ratios <- readRDS(paste0(source.dir,"Datasets/Colon_Cancer/Bian2018/data/scDNA/mapped_CNA_ratios.rds"))
CNA_ratios_diploid <- round(na.omit(mapped_CNA_ratios*2))
table(CNA_ratios_diploid[,1])


#### RNA data #### 
##non UMI
### scores imputed
nonUMI_CBS_pvalues <- readRDS(file=paste0(cellnet.dir,'Colon_Cancer/Bian/output/scaledata/CRC01_02/v2/without_covariates/genewise/FixedMovingWindow/CBS_pvalues_Matrix_Imputed.rds'))
nonUMI_data <- readRDS(file=paste0(cellnet.dir,'Colon_Cancer/Bian/output/scaledata/CRC01_02/v2/without_covariates/genewise/FixedMovingWindow/data.rds'))
nonUMI_data_metadata <- nonUMI_data@meta.data

#copykat
Bian_CRC01_CRC02_copykat <- read.delim("/cellfile/datapublic/rjohnen/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/copyKat/Colon_Cancer/output_Bian_CRC01_CRC02/Bian_CRC01_CRC02_copykat_CNA_raw_results_gene_by_cell.txt")
rownames(Bian_CRC01_CRC02_copykat) <- Bian_CRC01_CRC02_copykat$hgnc_symbol
Bian_CRC01_CRC02_copykat <- Bian_CRC01_CRC02_copykat[,-c(1:7)]
#inferCNV
infercnv.preliminary.observations <- read.csv("/cellfile/datapublic/rjohnen/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/inferCNV/Colon_Cancer/output/infercnv.preliminary.observations.txt", sep="")

###### remove chromosomes

gencode <- readRDS(file='/cellfile/cellnet/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_ZINBWaVE_NetCNA_CBS/Data/gencode_v28.rds')
genes_keep_not_clonal<-subset(gencode,gencode$seqid %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr9","chr10","chr11","chr12","chr14","chr16","chr19","chr21","chrX"))$gene_name
genes_keep_clonal<-subset(gencode,gencode$seqid %in% c("chr8","chr13","chr15","chr17","chr18","chr20","chr22"))$gene_name

intersect(genes_keep_clonal,rownames(nonUMI_CBS_pvalues))
intersect(genes_keep_not_clonal,rownames(nonUMI_CBS_pvalues))

cells_keep<-rownames(nonUMI_data_metadata[which(nonUMI_data_metadata$Donor == "Donor_1"),])

saveRDS(CNA_ratios_diploid[,cells_keep],"output_manuscript/Bian_DNA_CNA.rds")
saveRDS(nonUMI_CBS_pvalues[intersect(genes_keep_not_clonal,rownames(nonUMI_CBS_pvalues)),cells_keep],"output_manuscript/Bian_RNA_CNA_RNA2CNA.rds")
saveRDS(infercnv.preliminary.observations[intersect(genes_keep_not_clonal,rownames(infercnv.preliminary.observations)),cells_keep],"output_manuscript/Bian_RNA_CNA_infercnv.rds")
saveRDS(Bian_CRC01_CRC02_copykat[intersect(genes_keep_not_clonal,rownames(Bian_CRC01_CRC02_copykat)),cells_keep],"output_manuscript/Bian_RNA_CNA_copykat.rds")

p1<-ROC_Del_Amp_3Analysis(CNAsDNA = CNA_ratios_diploid[,cells_keep], 
                          PredCNAs1 = nonUMI_CBS_pvalues[intersect(genes_keep_not_clonal,rownames(nonUMI_CBS_pvalues)),cells_keep],
                          PredCNAs2 = infercnv.preliminary.observations[intersect(genes_keep_not_clonal,rownames(infercnv.preliminary.observations)),cells_keep],
                          PredCNAs3 = Bian_CRC01_CRC02_copykat[intersect(genes_keep_not_clonal,rownames(Bian_CRC01_CRC02_copykat)),cells_keep],
                          NamePredCNAs1 = "RNA2CNA",
                          NamePredCNAs2 = "InferCNV",
                          NamePredCNAs3 = "CopyKAT",
                          colors=c("#BB5566","#004488","#DDAA33"))

p2<-ROC_Del_Amp_3Analysis(CNAsDNA = CNA_ratios_diploid[,cells_keep], 
                          PredCNAs1 = nonUMI_CBS_pvalues[intersect(genes_keep_clonal,rownames(nonUMI_CBS_pvalues)),cells_keep],
                          PredCNAs2 = infercnv.preliminary.observations[intersect(genes_keep_clonal,rownames(infercnv.preliminary.observations)),cells_keep],
                          PredCNAs3 = Bian_CRC01_CRC02_copykat[intersect(genes_keep_clonal,rownames(Bian_CRC01_CRC02_copykat)),cells_keep],
                          NamePredCNAs1 = "RNA2CNA",
                          NamePredCNAs2 = "InferCNV",
                          NamePredCNAs3 = "CopyKAT",
                          colors=c("#BB5566","#004488","#DDAA33"))
source("Precred_Del_Amp.R")
library(precrec)
p1_prec<-Precrec_Del_Amp(CNAsDNA1 = CNA_ratios_diploid[,cells_keep]-2,
                         CNAsDNA2 = CNA_ratios_diploid[,cells_keep]-2,
                         CNAsDNA3 = CNA_ratios_diploid[,cells_keep]-2,
                         scRNA_df1 = nonUMI_CBS_pvalues[intersect(genes_keep_not_clonal,rownames(nonUMI_CBS_pvalues)),cells_keep],
                         scRNA_df2 = infercnv.preliminary.observations[intersect(genes_keep_not_clonal,rownames(infercnv.preliminary.observations)),cells_keep],
                         scRNA_df3 = Bian_CRC01_CRC02_copykat[intersect(genes_keep_not_clonal,rownames(Bian_CRC01_CRC02_copykat)),cells_keep],
                         NamePred1 = "RNA2CNA", NamePred2 = "InferCNV",NamePred3 = "CopyKat",
                         colors=c("#BB5566","#004488","#DDAA33"))

p2_prec<-Precrec_Del_Amp(CNAsDNA1 = CNA_ratios_diploid[,cells_keep]-2,
                        CNAsDNA2 = CNA_ratios_diploid[,cells_keep]-2,
                        CNAsDNA3 = CNA_ratios_diploid[,cells_keep]-2,
                        scRNA_df1 = nonUMI_CBS_pvalues[intersect(genes_keep_clonal,rownames(nonUMI_CBS_pvalues)),cells_keep],
                        scRNA_df2 = infercnv.preliminary.observations[intersect(genes_keep_clonal,rownames(infercnv.preliminary.observations)),cells_keep],
                        scRNA_df3 = Bian_CRC01_CRC02_copykat[intersect(genes_keep_clonal,rownames(Bian_CRC01_CRC02_copykat)),cells_keep],
                        NamePred1 = "RNA2CNA", NamePred2 = "InferCNV",NamePred3 = "CopyKat",
                        colors=c("#BB5566","#004488","#DDAA33"))


pdf('plots_manuscript/ROC_Curves_Bian2018_RNA2CNA_copykat_inferCNV_clonal_vs_notclonal.pdf',width=10)
ggarrange(NA,p2[[1]],p2[[2]],p2[[3]],
          NA,p1[[1]],p1[[2]],NA,
          ncol=4,nrow=2,
          labels = c("A","","","",
                     "B","","",""), 
          font.label = list(size = 30),widths = c(0.05,0.4,0.4,0.2))
dev.off()
pdf('plots_manuscript/PREC_Curves_Bian2018_RNA2CNA_copykat_inferCNV_clonal_vs_notclonal.pdf',width=20,height = 13)
ggarrange(NA,p2_prec[[1]],
          NA,p1_prec[[1]],
          ncol=2,nrow=2,
          labels = c("A","",
                     "B",""), 
          font.label = list(size = 30),widths = c(0.05,0.4))
dev.off()

######## DNTR #########
#### DNA Data ####
DNA_CNV_Genematrix_all_cells <- readRDS("/cellfile/datapublic/rjohnen/Analysis/ROC_Curves/DNTRseq/CNV_Genematrix_all_cells.rds")
table(DNA_CNV_Genematrix_all_cells[,1])
SRA_runtable <- readRDS(file='/cellfile/datapublic/rjohnen/Ronja/Tumor_Evolution/Analysis/Datasets/DNTR-seq_(SmartSeq2)/output/SraRunTable_DNA.rds')
DNA_metadata <- SRA_runtable[,c('Run','cell_ID','Sample.Name','Time_point','Treatment')]

DNA_metadata$Cell <- sapply(1:nrow(DNA_metadata), function(cell) strsplit(as.character(DNA_metadata$cell_ID)[cell],'[.]')[[1]][1])
DNA_metadata$Well <- sapply(1:nrow(DNA_metadata), function(cell) strsplit(as.character(DNA_metadata$cell_ID)[cell],'[.]')[[1]][2])
DNA_metadata$Sample <- paste0(DNA_metadata$Cell,'_',DNA_metadata$Well)

all.equal(as.character(DNA_metadata$Run), colnames(DNA_CNV_Genematrix_all_cells))
DNA_metadata$Run <- as.character(DNA_metadata$Run)

DNA_metadata <- DNA_metadata[which( DNA_metadata$Run %in% colnames(DNA_CNV_Genematrix_all_cells)  ),]
DNA_CNV_Genematrix_all_cells <- DNA_CNV_Genematrix_all_cells[,DNA_metadata$Run]
all.equal(as.character(DNA_metadata$Run), colnames(DNA_CNV_Genematrix_all_cells))
colnames(DNA_CNV_Genematrix_all_cells) <- DNA_metadata$Sample

#### RNA data #### 
### RNA2CNA ###
RNA_CNV_Genematrix_all_cells <- readRDS(file=paste0(cellnet.dir,'DNTR-seq/output/scaledata/v2/genewise/FixedMovingWindow/CBS_pvalues_Matrix_Imputed.rds'))
UMI_data <- readRDS(file=paste0(cellnet.dir,'DNTR-seq/output/scaledata/v2/genewise/FixedMovingWindow/data.rds'))
metadata <- UMI_data@meta.data
metadata <- droplevels(metadata)
head(metadata)
rm(UMI_data)

RNA_CNV_Genematrix_all_cells <- RNA_CNV_Genematrix_all_cells[,rownames(metadata)]
all.equal(rownames(metadata), colnames(RNA_CNV_Genematrix_all_cells))
colnames(RNA_CNV_Genematrix_all_cells) <- metadata$Sample_name

#copykat
DNTR_copykat <- read.delim("/cellfile/datapublic/rjohnen/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/copyKat/DNTR/output/DNTR_copykat_CNA_raw_results_gene_by_cell.txt")
rownames(DNTR_copykat) <- DNTR_copykat$hgnc_symbol
DNTR_copykat <- DNTR_copykat[,-c(1:7)]
rownames(metadata) <- gsub("\\-", ".", rownames(metadata))
cells <-intersect(rownames(metadata),colnames(DNTR_copykat))

DNTR_copykat <- DNTR_copykat[,cells]
all.equal(rownames(metadata[cells,]), colnames(DNTR_copykat))
colnames(DNTR_copykat) <- metadata[cells,]$Sample_name

#inferCNV
infercnv.preliminary.observations <- read.csv("/cellfile/datapublic/rjohnen/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/inferCNV/DNTR/output/infercnv.preliminary.observations.txt", sep="")
rownames(metadata) <- gsub("\\.", "_", rownames(metadata))
rownames(metadata) <- gsub("\\-", "_", rownames(metadata))
infercnv.preliminary.observations <- infercnv.preliminary.observations[,rownames(metadata)]
all.equal(rownames(metadata), colnames(infercnv.preliminary.observations))
colnames(infercnv.preliminary.observations) <- metadata$Sample_name

###### remove chromosomes

gencode <- readRDS(file='/cellfile/cellnet/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_ZINBWaVE_NetCNA_CBS/Data/gencode_v28.rds')
genes_keep_not_clonal<-subset(gencode,gencode$seqid %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr9","chr11","chr12","chr13","chr14","chr15","chr16","chr18","chr19","chr20","chr21","chr22"))$gene_name
genes_keep_clonal<-subset(gencode,gencode$seqid %in% c("chr8","chr10","chr17","chrX"))$gene_name

saveRDS(DNA_CNV_Genematrix_all_cells,"output_manuscript/Zacharidias_DNA_CNA.rds")
saveRDS(RNA_CNV_Genematrix_all_cells[intersect(genes_keep_not_clonal,rownames(RNA_CNV_Genematrix_all_cells)),],"output_manuscript/Zacharidias_RNA_CNA_RNA2CNA.rds")
saveRDS(infercnv.preliminary.observations[intersect(genes_keep_not_clonal,rownames(infercnv.preliminary.observations)),],"output_manuscript/Zacharidias_RNA_CNA_infercnv.rds")
saveRDS(DNTR_copykat[intersect(genes_keep_not_clonal,rownames(DNTR_copykat)),],"output_manuscript/Zacharidias_RNA_CNA_copykat.rds")

p3<-ROC_Del_Amp_3Analysis(CNAsDNA = DNA_CNV_Genematrix_all_cells, 
                      PredCNAs1 = RNA_CNV_Genematrix_all_cells[intersect(genes_keep_not_clonal,rownames(RNA_CNV_Genematrix_all_cells)),],
                      PredCNAs2 = infercnv.preliminary.observations[intersect(genes_keep_not_clonal,rownames(infercnv.preliminary.observations)),],
                      PredCNAs3 = DNTR_copykat[intersect(genes_keep_not_clonal,rownames(DNTR_copykat)),],
                      NamePredCNAs1 = "RNA2CNA",
                      NamePredCNAs2 = "InferCNV",
                      NamePredCNAs3 = "CopyKAT",
                      colors=c("#BB5566","#004488","#DDAA33"))

p4 <- ROC_Del_Amp_3Analysis(CNAsDNA = DNA_CNV_Genematrix_all_cells, 
                      PredCNAs1 = RNA_CNV_Genematrix_all_cells[intersect(genes_keep_clonal,rownames(RNA_CNV_Genematrix_all_cells)),],
                      PredCNAs2 = infercnv.preliminary.observations[intersect(genes_keep_clonal,rownames(infercnv.preliminary.observations)),],
                      PredCNAs3 = DNTR_copykat[intersect(genes_keep_clonal,rownames(DNTR_copykat)),],
                      NamePredCNAs1 = "RNA2CNA",
                      NamePredCNAs2 = "InferCNV",
                      NamePredCNAs3 = "CopyKAT",
                      colors=c("#BB5566","#004488","#DDAA33"))

p3_prec<-Precrec_Del_Amp(CNAsDNA1 = DNA_CNV_Genematrix_all_cells-2,
                         CNAsDNA2 = DNA_CNV_Genematrix_all_cells-2,
                         CNAsDNA3 = DNA_CNV_Genematrix_all_cells-2,
                         scRNA_df1 = RNA_CNV_Genematrix_all_cells[intersect(genes_keep_not_clonal,rownames(RNA_CNV_Genematrix_all_cells)),],
                         scRNA_df2 = infercnv.preliminary.observations[intersect(genes_keep_not_clonal,rownames(infercnv.preliminary.observations)),],
                         scRNA_df3 = DNTR_copykat[intersect(genes_keep_not_clonal,rownames(DNTR_copykat)),],
                         NamePred1 = "RNA2CNA", NamePred2 = "InferCNV",NamePred3 = "CopyKat",
                         colors=c("#BB5566","#004488","#DDAA33"))

p4_prec <- Precrec_Del_Amp(CNAsDNA1 = DNA_CNV_Genematrix_all_cells-2,
                           CNAsDNA2 = DNA_CNV_Genematrix_all_cells-2,
                           CNAsDNA3 = DNA_CNV_Genematrix_all_cells-2,
                           scRNA_df1 = RNA_CNV_Genematrix_all_cells[intersect(genes_keep_clonal,rownames(RNA_CNV_Genematrix_all_cells)),],
                           scRNA_df2 = infercnv.preliminary.observations[intersect(genes_keep_clonal,rownames(infercnv.preliminary.observations)),],
                           scRNA_df3 = DNTR_copykat[intersect(genes_keep_clonal,rownames(DNTR_copykat)),],
                           NamePred1 = "RNA2CNA", NamePred2 = "InferCNV",NamePred3 = "CopyKat",
                           colors=c("#BB5566","#004488","#DDAA33"))


    
pdf('plots_manuscript/ROC_Curves_DNTR_RNA2CNA_copykat_inferCNV_clonal_vs_notclonal.pdf',width=10)
ggarrange(NA,p4[[1]],p4[[2]],p4[[3]],
          NA,p3[[1]],p3[[2]],NA,
          ncol=4,nrow=2,
          labels = c("A","","","",
                     "B","","",""), 
          font.label = list(size = 30),widths = c(0.05,0.4,0.4,0.2))
dev.off()



pdf('plots_manuscript/PRECREC_Curves_DNTR_RNA2CNA_copykat_inferCNV_clonal_vs_notclonal.pdf',width=20,height = 13)
ggarrange(NA,p4_prec[[1]],
          NA,p3_prec[[1]],
          ncol=2,nrow=2,
          labels = c("A","",
                     "B",""), 
          font.label = list(size = 30),widths = c(0.05,0.4))
dev.off()

