source.dir <- "../../../../../"
source(paste0(source.dir,"packages.R"))
source(paste0(source.dir,"functions.R"))
source(paste0(source.dir,"Functions/GO_enrichment.R"))
output_dir <- c("output/FDR01/")
library(topGO)
library(biomaRt)
library(viridis)

gencode <- readRDS(file='/cellfile/cellnet/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_ZINBWaVE_NetCNA_CBS/Data/gencode_v28.rds')
gencode_gene_types<- readRDS(file='/cellfile/cellnet/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_ZINBWaVE_NetCNA_CBS/Data/gencode_gene_types.rds')

# correlate CNAs with genes → ranking of genes
#add up the number of CNAs per cell (segments)
#possibly per tissue, or per cell type if there are cell types with more # CNAs → Attention: here I have to take the normalized residuals

# all endothelial cells together, so that you have more signal (because each cell has many zeros) 
# → fit glm poission, because n CNAs is discrete (nCNA ~ Expr) 
# → sort genes by p-value (note: slope is not interpretable, because it can be easily scaled, low expr. genes have a higher slope than high expr. genes) 
# → y= beta * x (if x is large, beta becomes small and vice versa)

dir_data <- "../../Aging_Analysis_Segments/Pancreas/Combine_loop/output/FDR01/"
#dir_plot <- "UMI/FDR01/"


#### Load Gene expression data (after normalization)
### TS
Tabula_Sapiens <- readRDS(file = paste0('/cellfile/cellnet/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_CBS/Tabula_Sapiens/Pancreas/output/scaledata/v2/genewise/FixedMovingWindow/data.rds'))
dim(Tabula_Sapiens@assays$SCT$scale.data)
### Muraro
Muraro2016 <- readRDS(file = paste0('/cellfile/cellnet/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_CBS/Pancreas/Muraro2016/output/scaledata/genewise/FixedMovingWindow/data.rds'))
dim(Muraro2016@assays$SCT$scale.data)
### Baron
Baron2016 <- readRDS(file = paste0('/cellfile/cellnet/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_CBS/Pancreas/Baron2016/output/scaledata/genewise/FixedMovingWindow/data.rds'))
dim(Baron2016@assays$SCT$scale.data)
### Peng
Peng2019 <- readRDS(file = paste0('/cellfile/cellnet/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_CBS/Pancreas/Peng2019/output/scaledata/v2/genewise/healthy/FixedMovingWindow/data.rds'))
dim(Peng2019@assays$SCT$scale.data)

genes.intersect <- intersect(intersect(intersect(rownames(Tabula_Sapiens),rownames(Muraro2016)),rownames(Baron2016)),rownames(Peng2019))
length(genes.intersect)  
genes.keep <-unique(c(rownames(Tabula_Sapiens),rownames(Muraro2016),rownames(Baron2016),rownames(Peng2019)))
length(genes.keep)

Expr <- cbind(Tabula_Sapiens@assays$SCT$scale.data[genes.intersect,],Muraro2016@assays$SCT$scale.data[genes.intersect,], Baron2016@assays$SCT$scale.data[genes.intersect,],Peng2019@assays$SCT$scale.data[genes.intersect,])


####### Load Number CNAs per Gene 
Baron2016_Number_CNA_metadata <- readRDS(paste0(dir_data,'Baron2016_Number_CNA_metadata.rds'))
Baron2016_Number_CNA_metadata$Dataset <- 'Baron'
Baron2016_Number_CNA_metadata$Disease <- "no"
Baron2016_Number_CNA_metadata$Disease[which(Baron2016_Number_CNA_metadata$Age == 59)] <- "yes"
Peng2019_Number_CNA_metadata <- readRDS(paste0(dir_data,'Peng2019_Number_CNA_metadata.rds'))
Peng2019_Number_CNA_metadata$Dataset <- 'Peng'
Peng2019_Number_CNA_metadata_Normal <- subset(Peng2019_Number_CNA_metadata, Disease == 'Normal')
Peng2019_Number_CNA_metadata_Normal <- Peng2019_Number_CNA_metadata_Normal[,colnames(Baron2016_Number_CNA_metadata)]
Peng2019_Number_CNA_metadata_Normal$Disease <- 'no'
Muraro2016_Number_CNA_metadata <- readRDS(paste0(dir_data,'Muraro2016_Number_CNA_metadata.rds'))
Muraro2016_Number_CNA_metadata$Dataset <- 'Muraro'
Muraro2016_Number_CNA_metadata$Disease <- "no"
Muraro2016_Number_CNA_metadata <- Muraro2016_Number_CNA_metadata[,colnames(Baron2016_Number_CNA_metadata)]

TS_Pancreas_Number_CNA_metadata <- readRDS(paste0('../../Aging_Analysis_Segments/Pancreas/Combine_loop/output/TS/Pancreas/',strsplit(dir_data,"/")[[1]][7],'/Pancreas_Number_CNA_metadata.rds'))
TS_Pancreas_Number_CNA_metadata$Dataset <- 'TS_Pancreas'
TS_Pancreas_Number_CNA_metadata$Disease <- "no"
colnames(TS_Pancreas_Number_CNA_metadata)[8] <- "Donor"
colnames(TS_Pancreas_Number_CNA_metadata)[12] <- "Celltype"
TS_Pancreas_Number_CNA_metadata <- TS_Pancreas_Number_CNA_metadata[,colnames(Baron2016_Number_CNA_metadata)]

Number_CNA_metadata <- rbind(Muraro2016_Number_CNA_metadata,Baron2016_Number_CNA_metadata,Peng2019_Number_CNA_metadata_Normal, TS_Pancreas_Number_CNA_metadata) # Peng2019_Number_CNA_metadata_Normal
Number_CNA_metadata$Age <-as.numeric(as.character(Number_CNA_metadata$Age))

Number_CNA_metadata <- droplevels(Number_CNA_metadata)
levels(Number_CNA_metadata$Celltype)
Number_CNA_metadata$Celltype <- revalue(Number_CNA_metadata$Celltype, 
                                        c('alpha'='Alpha',"acinar"='Acinar','unsure' = 'Unclassified','delta'='Delta','beta'='Beta',
                                          'ductal'='Ductal',"Ductal cell" = 'Ductal', 
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
                                          "None/Other"='Unclassified',"Stellate"="Stellate","Ductal"="Ductal","Gamma/PP"="Gamma/PP",
                                          "endothelial cell"="Endothelial", "myeloid cell"="Myeloid","pancreatic acinar cell"="Acinar",
                                          "pancreatic ductal cell"="Ductal", "pancreatic stellate cell"="Stellate", "t cell"="T cell"))
levels(Number_CNA_metadata$Celltype)
levels(Number_CNA_metadata$Sex)
Number_CNA_metadata$Sex <- revalue(Number_CNA_metadata$Sex, c('m'='m','f'='f','Male'='m','Female'='f',
                                                              'F'='f','M'='m'))
table(Number_CNA_metadata$Celltype,Number_CNA_metadata$Dataset)
#Number_CNA_metadata <- subset(Number_CNA_metadata, Number_CNA_metadata$Celltype %in% c("Alpha","Beta","Ductal","Acinar","Delta","Endothelial")) 
# rm < 20 years
Number_CNA_metadata <- subset(Number_CNA_metadata, Number_CNA_metadata$Age > 15) 

#celltype<-levels(Number_CNA_metadata$Celltype)
celltype<-c("Alpha","Beta","Ductal","Acinar","Delta","Endothelial","Myeloid","Macrophage")

ggplot(Number_CNA_metadata)+geom_violin(aes(x=Celltype,y=Number_CNAs_per_cell))

celltype <- "Stellate"
# subset endothelial cells
Endo <- subset(Number_CNA_metadata,Number_CNA_metadata$Celltype == celltype) #  Acinar             Ductal 
summary(Endo$Number_CNAs_per_cell)
# dim(Endo)
# dim(Number_CNA_metadata)
# dim(Expr)

cells_keep <- intersect(rownames(Endo),colnames(Expr))
#MODEL
Endo_pval <- lapply(1:nrow(Expr), function(gene) rbind(cor.test(Endo[cells_keep,]$Number_CNAs_per_cell, Expr[gene,cells_keep], method = "spearman", exact=FALSE)[c("estimate", "p.value")]) )
#str(Endo_pval)
Endo_pval<-do.call(rbind,Endo_pval)
rownames(Endo_pval) <- rownames(Expr)
#hist(Endo_pval[,2])
Endo_pval <- data.frame('coef'=as.numeric(Endo_pval[,"estimate"]),"pval"=as.numeric(Endo_pval[,"p.value"]),"gene"=rownames(Endo_pval))
Endo_pval$pval_adj<- p.adjust(Endo_pval$pval,method = "fdr")

# TOP 5 HITS
saveRDS(Endo_pval, file=paste0("output/spearman/",celltype,"/Gene_list.rds"))
Endo_pval[order(Endo_pval$pval_adj),][1:5,"gene"] 
Endo_pval[order(Endo_pval$pval_adj),][1:5,]

#endo
#  Neg"FABP5"  "HMGB1"  "IGFBP4" "PRMT1" 
# Pos"KCTD12" 


# Neg "TPST1"   "LPCAT1"  "CYP27A1" "THBS1"   "DEDD"  
# Pos

############# GO TERM ENRICHMENT

Endo_pval_sig_pos <- subset(Endo_pval, pval_adj <= 0.05 & coef > 0 ) #& coef > 0 
Endo_pval_sig_neg <- subset(Endo_pval, pval_adj <= 0.05 & coef < 0 ) #& coef > 0 
Endo_pval_sig_all <- subset(Endo_pval, pval_adj <= 0.05 )

sig_genes <- list()
sig_genes$All <- Endo_pval_sig_all$gene
sig_genes$Neg <- Endo_pval_sig_neg$gene
sig_genes$Pos <- Endo_pval_sig_pos$gene

# go enrichment
allGenes <- genes.intersect
t <- list()
for (i in 1:length (sig_genes)){
  ClusterGenes <- function(x){
    y <- allGenes
    return(as.integer(y %in% x))
  }
  t[[i]]<- ClusterGenes(sig_genes[[i]])
}

tt <- do.call(cbind,t) #makes a table
rownames (tt) <- allGenes
all_genes_clusterwise <- tt

# # make named factor showing which genes are of interest
# geneList=factor(as.integer(allGenes %in% sig_genes$Endothelial))
# names(geneList)= allGenes

EPPId2GO <- readRDS('output/bioMart/EPPId2GO_2024-02-11.rds')

# Remove blank entries
EPPId2GO <- EPPId2GO[EPPId2GO$go_id != '',]
# convert from table format to list format
PPID2GO <- by(EPPId2GO$go_id,EPPId2GO$hgnc_symbol,function(x) as.character(x))

#all_genes_clusterwise<- apply(all_genes_clusterwise, 2, as.factor)
sig_genes$GOdata$BP <- apply(all_genes_clusterwise,2,run_topgoBP)

### Performing statistical tests for the genes in the individual clusters:
# BP
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
sig_genes$ClassicBP <-lapply(seq_along(sig_genes$GOdata$BP), function(x) getSigGroups(sig_genes$GOdata$BP[[x]], test.stat))

test.stat <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher test")
sig_genes$ElimBP <-lapply(seq_along(sig_genes$GOdata$BP), function(x) getSigGroups(sig_genes$GOdata$BP[[x]], test.stat))

sig_genes$allResClustersBP_EF<- lapply(seq_along(sig_genes$GOdata$BP), function(x) GenTable(sig_genes$GOdata$BP[[x]],classicFisher= sig_genes$ClassicBP[[x]],
                                                                                            elimFisher= sig_genes$ElimBP[[x]],orderBy = "elimFisher", topNodes = 300))

names (sig_genes$allResClustersBP_EF)<- c("All",'Neg','Pos')


sig_genes$allResClustersBP_EF <- lapply(1:length(sig_genes$allResClustersBP_EF), function (x) compute_log2SignDivExp(sig_genes$allResClustersBP_EF[[x]]))
#View(sig_genes$allResClustersBP_EF[[1]])

### Adding the GO domains on the enrichment tables
for (i in 1:length(sig_genes$allResClustersBP_EF)){
  sig_genes$allResClustersBP_EF[[i]]$elimFisher <- as.numeric( gsub( "<","", sig_genes$allResClustersBP_EF[[i]]$elimFisher))
  sig_genes$allResClustersBP_EF[[i]]$domain <- rep("BP",nrow(sig_genes$allResClustersBP_EF[[i]]))
}

sig_genes$allResClustersBP_EF <- lapply(1:length(sig_genes$allResClustersBP_EF), function (x) AddingGeneNamesOnGO(sig_genes$allResClustersBP_EF[[x]], sig_genes$GOdata$BP[[x]], sig_genes[[x]]))
str(sig_genes$allResClustersBP_EF[[1]])

sig_genes$allResClustersBP_EF_filtered <- lapply(1:length(sig_genes$allResClustersBP_EF), function(x) filterTopGo(sig_genes$allResClustersBP_EF[[x]], mingenes=10, alpha=0.05, alphaTerm="elimFisher"))

All <- sig_genes$allResClustersBP_EF_filtered[[1]]
Neg <- sig_genes$allResClustersBP_EF_filtered[[2]]
Pos <- sig_genes$allResClustersBP_EF_filtered[[3]]

saveRDS(All, file=paste0("output/spearman/",celltype,"/All.rds"))
saveRDS(Neg, file=paste0("output/spearman/",celltype,"/Neg.rds") )
saveRDS(Pos, file=paste0("output/spearman/",celltype,"/Pos.rds")) 

#Pos <- readRDS(file=paste0("output/spearman/",celltype,"/Pos.rds")) 

#Endothelial<-readRDS( file="output/Endothelial_pos_coeff.rds")

######### Visualization
All <- All[which (All$Significant >= 10 & All$Enrichment >= 1),]
All <- All[order(All$Enrichment),]
All$Term <- as.factor(All$Term)
#levels(Endothelial$Term)
All$Term <- factor(All$Term, levels = All$Term)

# ggplot(Endothelial, aes(x=domain, y=Term, color = elimFisher, size = Enrichment)) + 
#   geom_point()+ scale_color_viridis(direction = -1, limits=c(0,0.1), na.value="#440154FF") # -log10(p.adjust)
png(paste0('plots/spearman/',celltype,'/All_GO_enrichment.png'))
ggplot(All, aes(x=Enrichment, y=Term, color = elimFisher, size = Significant)) +  theme_readable() +
  geom_point()+ scale_color_viridis(direction = -1, limits=c(0,0.05), na.value="#440154FF") # -log10(p.adjust) # Enrichment = log2(Signifcant/Expected)
dev.off()
Neg <- Neg[which (Neg$Significant >= 10 & Neg$Enrichment >= 1),]
Neg <- Neg[order(Neg$Enrichment),]
Neg$Term <- as.factor(Neg$Term)
#levels(Endothelial$Term)
Neg$Term <- factor(Neg$Term, levels = Neg$Term)
png(paste0('plots/spearman/',celltype,'/Neg_GO_enrichment.png'))
ggplot(Neg, aes(x=Enrichment, y=Term, color = elimFisher, size = Significant)) +  theme_readable() +
  geom_point()+ scale_color_viridis(direction = -1, limits=c(0,0.05), na.value="#440154FF") # -log10(p.adjust) # Enrichment = log2(Signifcant/Expected)
dev.off()
Pos <- Pos[which (Pos$Significant >= 10 & Pos$Enrichment >= 1),]
Pos <- Pos[order(Pos$Enrichment),]
Pos$Term <- as.factor(Pos$Term)
#levels(Endothelial$Term)
Pos$Term <- factor(Pos$Term, levels = Pos$Term)
png(paste0('plots/spearman/',celltype,'/Pos_GO_enrichment.png'))
ggplot(Pos, aes(x=Enrichment, y=Term, color = elimFisher, size = Significant)) + theme_readable() +
  geom_point()+ scale_color_viridis(direction = -1, limits=c(0,0.05), na.value="#440154FF") # -log10(p.adjust) # Enrichment = log2(Signifcant/Expected)
dev.off()


pdf(paste0('plots/spearman/',celltype,'/All_GO_enrichment.pdf'),height = 4,width = 7)
ggplot(All, aes(x=Enrichment, y=Term, color = elimFisher, size = Significant)) +  theme_readable() +
  geom_point()+ scale_color_viridis(direction = -1, limits=c(0,0.05), na.value="#440154FF") # -log10(p.adjust) # Enrichment = log2(Signifcant/Expected)
dev.off()
pdf(paste0('plots/spearman/',celltype,'/Neg_GO_enrichment.pdf'),height = 4,width = 7)
ggplot(Neg, aes(x=Enrichment, y=Term, color = elimFisher, size = Significant)) +  theme_readable() +
  geom_point()+ scale_color_viridis(direction = -1, limits=c(0,0.05), na.value="#440154FF") # -log10(p.adjust) # Enrichment = log2(Signifcant/Expected)
dev.off()
pdf(paste0('plots/spearman/',celltype,'/Pos_GO_enrichment.pdf'),height = 4,width = 7)
ggplot(Pos, aes(x=Enrichment, y=Term, color = elimFisher, size = Significant)) +  theme_readable() +
  geom_point()+ scale_color_viridis(direction = -1, limits=c(0,0.05), na.value="#440154FF") # -log10(p.adjust) # Enrichment = log2(Signifcant/Expected)
dev.off()
