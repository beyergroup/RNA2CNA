source.dir <- "../../../../"
source(paste0(source.dir,"packages.R"))
source(paste0(source.dir,"functions.R"))
source(paste0(source.dir,"Functions/Mutation_Frequency_per_Gene.R"))
source(paste0(source.dir,"Functions/Lineplot.R"))

### Load all Pancreas datasets
CNA_Matrix_comb <- readRDS(file='output/FDR01/CNA_Matrix_comb.rds')
metadata_comb <- readRDS(file='output/FDR01/metadata_comb.rds')
CNA_Matrix_comb <- CNA_Matrix_comb[,-1]
dim(CNA_Matrix_comb)
dim(metadata_comb)

# add TS Pancreas
TS_Pancreas <- readRDS(file=paste0('output/TS/FDR01/Pancreas_CBS_pvalues_FDR01_Matrix_Imputed.rds'))
data <- readRDS(file = paste0('/cellfile/cellnet/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_CBS/Tabula_Sapiens/Pancreas/output/scaledata/v2/genewise/FixedMovingWindow/data.rds'))
metadata <- data@meta.data
colnames(metadata)[6] <- "Donor"
colnames(metadata)[10] <- "Celltype"
CNA_Matrix_comb <- cbind(CNA_Matrix_comb,TS_Pancreas)
metadata_comb <- rbind(metadata_comb,metadata[,c('orig.ident', 'nCount_RNA', 'nFeature_RNA',  'Donor' , 'Celltype', 'n_counts_UMIs', 'n_genes', 'nCount_SCT', 'nFeature_SCT')])
##
dim(CNA_Matrix_comb)
dim(metadata_comb)

metadata_comb <- droplevels(metadata_comb)
levels(metadata_comb$Celltype)
metadata_comb$Celltype <- revalue(metadata_comb$Celltype, 
                                  c('alpha'='Alpha',"acinar"='Acinar','unsure' = 'Unclassified','delta'='Delta','beta'='Beta',
                                    'ductal'='Ductal','Ductal cell'='Ductal',
                                    'mesenchymal' = 'Mesenchymal','Alpha'='Alpha','PP'='PP','Beta'='Beta','Delta'='Delta',
                                    'Acinar'='Acinar','Duct'='Ductal', 'Mesenchym' = 'Mesenchymal','Epsilon' = 'Epsilon',
                                    "Endo" ='Endothelial', 'activated_stellate'='Activated stellate', 'endothelial'='Endothelial',
                                    'epsilon'='Epsilon', "gamma"="Gamma","macrophage"="Macrophage","mast"="Mast",
                                    "quiescent_stellate"="Quiescent stellate","schwann"="Schwann",
                                    "Acinar cell"="Acinar", "B cell"= "B cell", "Ductal cell type 1"='Ductal',
                                    "Endocrine cell"="Endocrine","Endothelial cell"="Endothelial","Fibroblast cell"="Fibroblast",
                                    "Macrophage cell"= "Macrophage","Stellate cell"="Stellate","T cell"="T cell",
                                    "co-expression"="Co-expression","MHC class II cell"="MHC class II cell","not applicable"='Unclassified',
                                    "PSC"="PSC","unclassified cell"='Unclassified',"unclassified endocrine cell"="Unclassified endocrine",
                                    "None/Other"='Unclassified',"Stellate"="Stellate","Ductal"="Ductal","Gamma/PP"="Gamma/PP",                                              "endothelial cell"="Endothelial", "myeloid cell"="Myeloid","pancreatic acinar cell"="Acinar",
                                    "pancreatic ductal cell"="Ductal", "pancreatic stellate cell"="Stellate", "t cell"="T cell"))

#metadata_comb <- subset(metadata_comb, metadata_comb$Celltype %in% c("Alpha","Beta","Ductal","Acinar","Delta","Endothelial","Macrophage")) 
levels(metadata_comb$Celltype)

metadata_comb$Age <- metadata_comb$Donor
levels(metadata_comb$Age)
metadata_comb$Age <- revalue(metadata_comb$Age, c('Donor1'=17,'Donor2'=51,'Donor3'=38,'Donor4'=59,
                                                  'Donor_1'=23,'Donor_2'=48,'Donor_3'=54,'Donor_4'=59,
                                                  "N1"="30", "N2"="31","N3"="34", "N4"="41","N5"="42", "N6"="50","N7"="52", "N8"="53","N9"="55", "N10"="64","N11"="65", 
                                                  "T1"="36","T2"="44", "T3"="51","T4"="52","T5"="54", "T6"="54","T7"="54","T8"="56", "T9"="58","T10"="58","T11"="59", "T12"="59",
                                                  "T13"="59","T14"="61", "T15"="64","T16"="64","T17"="65", "T18"="66","T19"="67","T20"="67", "T21"="68","T22"="70","T23"="71", 
                                                  "T24"="72",
                                                  'TSP1'=59,'TSP9'=37))



metadata_comb$Age <- as.numeric(as.character(metadata_comb$Age))
summary(metadata_comb$Age)
metadata_comb$AgeGroup <- NA
metadata_comb[which(metadata_comb$Age < median(metadata_comb$Age)),'AgeGroup'] <- 'Young'
metadata_comb[which(metadata_comb$Age >= median(metadata_comb$Age)),'AgeGroup'] <- 'Old'
table(metadata_comb$AgeGroup)

# per donor 
levels(metadata_comb$orig.ident)
levels(metadata_comb$Donor)
table(metadata_comb$orig.ident,metadata_comb$Donor)
metadata_comb$Dataset <- metadata_comb$Donor
metadata_comb$Dataset <- revalue(metadata_comb$Dataset, c('Donor1'='Baron','Donor2'='Baron','Donor3'='Baron','Donor4'='Baron',
                                                          'Donor_1'='Muraro','Donor_2'='Muraro','Donor_3'='Muraro','Donor_4'='Muraro',
                                                          'N1'='Peng','N2'='Peng','N3'='Peng','N4'='Peng','N5'='Peng','N6'='Peng',
                                                          'N7'='Peng','N8'='Peng','N9'='Peng','N10'='Peng','N11'='Peng','TSP1'='TS','TSP9'='TS'))
donors <- levels(metadata_comb$Donor)
table(metadata_comb$Celltype,metadata_comb$Dataset)

celltypes <- table(metadata_comb$Celltype,metadata_comb$Dataset)
celltypes <- rownames(celltypes)[apply(celltypes,1, function(row) sum(row) > 100 & length(which(row !=0)) >= 2)]
celltypes <- celltypes[ !(celltypes %in% c('Macrophage','Stellate','Delta')) ] # rm because after filter for #CNAs only one dataset

table(metadata_comb$Celltype,metadata_comb$Donor)

metadata_comb <- subset(metadata_comb,metadata_comb$Celltype %in% celltypes)
metadata_comb <- droplevels(metadata_comb)
CNA_Matrix_comb <- CNA_Matrix_comb[,rownames(metadata_comb)]
dim(CNA_Matrix_comb)



###### FREQ
gencode_gene_types <- readRDS(file='/cellfile/cellnet/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_ZINBWaVE_NetCNA_CBS/Data/gencode_gene_types.rds')
metadata_comb <- subset(metadata_comb,metadata_comb$AgeGroup == "Old")
CNA_Matrix_comb <- CNA_Matrix_comb[,rownames(metadata_comb)]
Freq<-Mutation_Frequency_per_Gene(CNA_Matrix_comb,metadata_comb) # FÜr jedes Gen wird berechnet wie oft das amp oder del wurde
Freq_list <- list("Pancreas"=NULL)
  if(all.equal(rownames(Freq),gencode_gene_types$gene_name)){
    Freq_list$Pancreas <- gencode_gene_types
    Freq_list$Pancreas$Amp <- Freq[,"Amp"]
    Freq_list$Pancreas$Del <- Freq[,"Del"]
    Freq_list$Pancreas$Percentage_NAs <-  apply(CNA_Matrix_comb,1, function(row) (length(which(is.na(row)))/ncol(CNA_Matrix_comb))*100)
    Freq_list$Pancreas$tissue <- "Pancreas"
    Freq_list$Pancreas[which(Freq_list$Pancreas$Percentage_NAs > 50),"Amp"] <- NA
    Freq_list$Pancreas[which(Freq_list$Pancreas$Percentage_NAs > 50),"Del"] <- NA
    
    # rm p arm of chr 6
    Freq_list$Pancreas[which(Freq_list$Pancreas$chr == 6 & Freq_list$Pancreas$start <  59.8*10^6),"Amp"] <- NA
    Freq_list$Pancreas[which(Freq_list$Pancreas$chr == 6 & Freq_list$Pancreas$start <  59.8*10^6),"Del"] <- NA
    
    Freq_list$Pancreas<-Freq_list$Pancreas[which(Freq_list$Pancreas$chr %in% c(1:21,22) ),] # 21 hat nur NAs daher raus
    Freq_list$Pancreas<-droplevels(Freq_list$Pancreas)
  }
saveRDS(Freq_list,'output/Freq_list_TS_rm6p_old.rds')



# per celltype
Freq_list <- list()
for (celltype in celltypes){
  metadata_celltype <- subset(metadata_comb,metadata_comb$Celltype == celltype)
  CNA_Matrix_celltype <- CNA_Matrix_comb[,rownames(metadata_celltype)]
  Freq<-Mutation_Frequency_per_Gene(CNA_Matrix_celltype,metadata_celltype) # FÜr jedes Gen wird berechnet wie oft das amp oder del wurde
  
  if(all.equal(rownames(Freq),gencode_gene_types$gene_name)){
    Freq_list[[celltype]] <- gencode_gene_types
    Freq_list[[celltype]]$Amp <- Freq[,"Amp"]
    Freq_list[[celltype]]$Del <- Freq[,"Del"]
    Freq_list[[celltype]]$Percentage_NAs <-  apply(CNA_Matrix_celltype,1, function(row) (length(which(is.na(row)))/ncol(CNA_Matrix_celltype))*100)
    Freq_list[[celltype]]$tissue <- celltype
    Freq_list[[celltype]][which(Freq_list[[celltype]]$Percentage_NAs > 50),"Amp"] <- NA
    Freq_list[[celltype]][which(Freq_list[[celltype]]$Percentage_NAs > 50),"Del"] <- NA
    
    # rm p arm of chr 6
    Freq_list[[celltype]][which(Freq_list[[celltype]]$chr == 6 & Freq_list[[celltype]]$start <  59.8*10^6),"Amp"] <- NA
    Freq_list[[celltype]][which(Freq_list[[celltype]]$chr == 6 & Freq_list[[celltype]]$start <  59.8*10^6),"Del"] <- NA
    
    Freq_list[[celltype]]<-Freq_list[[celltype]][which(Freq_list[[celltype]]$chr %in% c(1:21,22) ),] # 21 hat nur NAs daher raus
    Freq_list[[celltype]]<-droplevels(Freq_list[[celltype]])
  }
}
names(Freq_list)<- celltypes
saveRDS(Freq_list,'output/Freq_list_TS_rm6p_per_celltype.rds')



Freq_list <- readRDS("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Onkogenic_CNAs/output/Freq_list_TS_rm6p_per_celltype.rds")

lineplot_CNAs(Del_Amp_profile = Freq_list,
              file_name = "per_celltype", 
              ylimit = NULL, 
              per_celltype = T,
              pdf=TRUE)

Freq_list <- readRDS("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Onkogenic_CNAs/output/Freq_list_TS_rm6p.rds")
Freq_list <- Freq_list$Pancreas

lineplot_CNAs(Del_Amp_profile = Freq_list,
              file_name = "Combined_Pancreas", 
              ylimit = NULL, 
              per_celltype = F,
              pdf=TRUE)

Freq_list <- readRDS("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Onkogenic_CNAs/output/Freq_list_TS_rm6p_young.rds")
Freq_list <- Freq_list$Pancreas

lineplot_CNAs(Del_Amp_profile = Freq_list,
              file_name = "Combined_young", 
              ylimit = NULL, 
              per_celltype = F,
              pdf=TRUE)


#### all tissue from TS

Freq_list<-readRDS('../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Overview_Plots/output/Freq_list_TS_rm6p_all_tissues.rds')

lineplot_CNAs(Del_Amp_profile = Freq_list,
              file_name = "Combined", 
              ylimit = NULL, 
              per_celltype = F,
              pdf=TRUE)

