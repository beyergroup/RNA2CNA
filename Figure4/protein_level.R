source.dir <- "../../../../../"
source(paste0(source.dir,"packages.R"))
source(paste0(source.dir,"functions.R"))
source(paste0(source.dir,"Functions/GO_enrichment.R"))
library(topGO)
library(biomaRt)
library(viridis)
library(readr)
library("rjson")
library(XML)
library(readxl)

ribosomal_genes <- read.delim2("/cellfile/datapublic/rjohnen/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Correlation/ribosomal_genes.txt")
ribosomal_genes <- subset(ribosomal_genes,ribosomal_genes$Group.name %in% c('L ribosomal proteins', 'S ribosomal proteins'))$Approved.symbol


Model <- read.csv("/cellfile/datapublic/rjohnen/Ronja/Tumor_Evolution/Analysis/Datasets/CCLE/DepMap/Model.csv")
OmicsAbsoluteCNGene <- read.csv("/cellfile/datapublic/rjohnen/Ronja/Tumor_Evolution/Analysis/Datasets/CCLE/DepMap/OmicsAbsoluteCNGene.csv")
summary(as.numeric(OmicsAbsoluteCNGene[1,]))
OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected <- read.csv("/cellfile/datapublic/rjohnen/Ronja/Tumor_Evolution/Analysis/Datasets/CCLE/DepMap/OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv")
OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected[1:5,1:5]

#Proteomics
protein <- read.csv("/cellfile/datapublic/rjohnen/Ronja/Tumor_Evolution/Analysis/Datasets/CCLE/DepMap/protein_quant_current_normalized.csv")
protein[1:5,1:10]
Table_S1_Sample_Information <- read_excel("/cellfile/datapublic/rjohnen/Ronja/Tumor_Evolution/Analysis/Datasets/CCLE/DepMap/Table_S1_Sample_Information.xlsx",  sheet = "Sample_Information")

cell_lines_names <- sapply(colnames(protein[49:ncol(protein)]), function(x) strsplit(x,"\\_")[[1]][1])

# subset 
all.equal(OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected$X,OmicsAbsoluteCNGene)
dim(OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected)
dim(OmicsAbsoluteCNGene)

RNA_gene_names <- sapply(2:ncol(OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected), function(x) strsplit(colnames(OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected)[x],"\\..")[[1]][1])
CNA_gene_names <- sapply(2:ncol(OmicsAbsoluteCNGene), function(x) strsplit(colnames(OmicsAbsoluteCNGene)[x],"\\..")[[1]][1])

RNA_data <- t(OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected[,-1])
colnames(RNA_data) <- OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected$X
rownames(RNA_data) <- RNA_gene_names

CNA_data <- t(OmicsAbsoluteCNGene[,-1])
colnames(CNA_data) <- OmicsAbsoluteCNGene$X
rownames(CNA_data) <- CNA_gene_names

Protein_data <- protein[,names(cell_lines_names)]
colnames(Protein_data) <- cell_lines_names
Protein_genes<- protein[,c( "Protein_Id","Gene_Symbol"  )]
  
Protein_genes$Gene_Symbol[duplicated(Protein_genes$Gene_Symbol)]

# intersect genes and sample names for all levels: RNA, CNA and protein
genes <- intersect(rownames(CNA_data),rownames(RNA_data))
genes <- intersect(Protein_genes$Gene_Symbol,genes)
ribosomal_genes <- intersect(ribosomal_genes,genes)

samples <- intersect(colnames(CNA_data),colnames(RNA_data))
metadata <- Model[Model$ModelID %in% samples,]

cell_lines_names <- intersect(cell_lines_names,metadata$StrippedCellLineName)
metadata <- metadata[metadata$StrippedCellLineName %in% cell_lines_names,]

# subset genes and samples
CNA_data <- CNA_data[genes,metadata$ModelID]
RNA_data <- RNA_data[genes,metadata$ModelID]
colnames(CNA_data) <- metadata$StrippedCellLineName
colnames(RNA_data) <- metadata$StrippedCellLineName

Protein_data <- Protein_data[,metadata$StrippedCellLineName]

# susbet genes protein

Protein_data <- Protein_data[!duplicated(Protein_genes$Gene_Symbol),]
Protein_genes <- Protein_genes[!duplicated(Protein_genes$Gene_Symbol),]
rownames(Protein_data)<-Protein_genes$Gene_Symbol
Protein_data <- Protein_data[genes,]
all.equal(rownames(Protein_data),rownames(CNA_data))
all.equal(colnames(Protein_data),colnames(CNA_data))

samples_without_del_in_ribosomes <- readRDS("output/samples_without_del_in_ribosomes.rds")
intersect(samples_without_del_in_ribosomes,colnames(Protein_data))

##### repeat GO with protein expr
CNA_table <- list()
CNA_table$Amp <- apply(CNA_data, 2, function(x) length( which(x > 2)))
CNA_table$Del <- apply(CNA_data, 2, function(x) length( which(x < 2)))
CNA_table <- data.frame(CNA_table)
head(CNA_table)
all.equal(colnames(Protein_data),rownames(CNA_table))
Protein_data<-as.matrix(Protein_data)


#MODEL
Protein_data <- Protein_data[,samples_without_del_in_ribosomes]
n_NA <- apply(Protein_data,1,  function(x) length(which(is.na(x))) )
table(n_NA)

Protein_data <- Protein_data[n_NA < ncol(Protein_data)/2,]

#PT_pval <- lapply(1:nrow(RNA_data), function(gene) rbind(cor.test(rowSums(CNA_table[,c("Del","Amp")]), Protein_data[gene,], method = "spearman", exact=FALSE)[c("estimate", "p.value")]) )
PT_pval <- lapply(1:nrow(Protein_data), function(gene) rbind(cor.test(CNA_table[samples_without_del_in_ribosomes,"Amp"], Protein_data[gene, samples_without_del_in_ribosomes], method = "spearman", exact=FALSE)[c("estimate", "p.value")]) )
#PT_pval <- lapply(1:nrow(Protein_data), function(gene) rbind(cor.test(rowSums(CNA_table[samples_without_del_in_ribosomes,c("Del","Amp")]), Protein_data[gene, samples_without_del_in_ribosomes], method = "spearman", exact=FALSE)[c("estimate", "p.value")]) )

PT_pval<-do.call(rbind,PT_pval)

# rbind(cor.test(CNA_table[samples_without_del_in_ribosomes,"Amp"], Protein_data[947, samples_without_del_in_ribosomes], method = "spearman", exact=FALSE)[c("estimate", "p.value")])
# 
# test <- cor.test(CNA_table[samples_without_del_in_ribosomes,"Amp"], Protein_data[947, samples_without_del_in_ribosomes], method = "spearman", exact=FALSE)
# test$p.value
# test$estimate
# 
# cor(CNA_table[samples_without_del_in_ribosomes,"Amp"], Protein_data[947, samples_without_del_in_ribosomes], method = "spearman")
# 
# plot(CNA_table[samples_without_del_in_ribosomes,"Amp"], Protein_data["OSMR", samples_without_del_in_ribosomes])

#str(PT_pval)
rownames(PT_pval) <- rownames(Protein_data)
#hist(PT_pval[,2])
PT_pval <- data.frame('coef'=as.numeric(PT_pval[,"estimate"]),"pval"=as.numeric(PT_pval[,"p.value"]),"gene"=rownames(PT_pval))
PT_pval$pval_adj<- p.adjust(PT_pval$pval,method = "fdr")


### subset ribosomal proteins
PT_pval$Ribo <- "No"
PT_pval[PT_pval$gene %in% ribosomal_genes,"Ribo"] <- "Yes"
PT_pval$Sig <- "No"
PT_pval[PT_pval$pval_adj < 0.05 ,"Sig"] <- "Yes"

#saveRDS(PT_pval,'output/Spearman_cor_Number_Amp_vs_protein_abundance_CCLE_samples_without_del_in_ribosomal_proteins.rds')

png(paste0('plots/Amp_Hist_coeff.png'))
ggplot(PT_pval)+geom_histogram(aes(x=coef,fill=Ribo),alpha=0.2, position="identity")+theme_readable()+facet_wrap(.~Ribo, scales = "free")+xlim(c(-0.35,0.3))
dev.off()
png(paste0('plots/Amp_Dens_coeff.png'))
ggplot(PT_pval)+geom_density(aes(x=coef,fill=Ribo),alpha=0.2, position="identity")+theme_readable()+geom_vline(xintercept = -0.13)#+facet_wrap(.~Sig)
dev.off()

subset(PT_pval,PT_pval$Ribo=="Yes")[which(subset(PT_pval,PT_pval$Ribo=="Yes")$coef < -0.13),]

t.test(coef~Ribo,PT_pval) #  p-value < 2.2e-16
t.test(subset(PT_pval,PT_pval$Ribo=="Yes")$coef,subset(PT_pval,PT_pval$Ribo=="No")$coef) #  p-value < 2.2e-16
wilcox.test(subset(PT_pval,PT_pval$Ribo=="Yes")$coef,subset(PT_pval,PT_pval$Ribo=="No")$coef) #  p-value < 2.2e-16

# TOP 5 HITS
saveRDS(PT_pval, file=paste0("output/Protein_list.rds"))
PT_pval[order(PT_pval$pval_adj),][1:5,"gene"] 
PT_pval[order(PT_pval$pval_adj),][1:5,]

png(paste0('plots/Hist_pvalues.png'))
ggplot(PT_pval)+geom_histogram(aes(x=pval),alpha=0.2, position="identity")+theme_readable()
dev.off()
png(paste0('plots/Hist_coeff.png'))
ggplot(PT_pval)+geom_histogram(aes(x=coef),alpha=0.2, position="identity")+theme_readable()
dev.off()
ggplot(PT_pval)+geom_point(aes(x=coef,y= -log10(pval)),alpha=0.2)+theme_readable()

############# GO TERM ENRICHMENT
PT_pval_sig_pos <- subset(PT_pval, pval_adj <= 0.05 & coef > 0 ) #& coef > 0 
PT_pval_sig_neg <- subset(PT_pval, pval_adj <= 0.05 & coef < 0 ) #& coef > 0 
PT_pval_sig_all <- subset(PT_pval, pval_adj <= 0.05 )

sig_genes <- list()
sig_genes$All <- PT_pval_sig_all$gene
sig_genes$Neg <- PT_pval_sig_neg$gene
sig_genes$Pos <- PT_pval_sig_pos$gene

# go enrichment
allGenes <-  rownames(RNA_data)
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

EPPId2GO <- readRDS('../Pancreas/output/bioMart/EPPId2GO_2024-02-11.rds')

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

saveRDS(All, file=paste0("output/All_Protein.rds"))
saveRDS(Neg, file=paste0("output/Neg_Protein.rds") )
saveRDS(Pos, file=paste0("output/Pos_Protein.rds")) 


######### Visualization
All <- All[which (All$Significant >= 10 & All$Enrichment >= 1),]
All <- All[order(All$Enrichment),]
All$Term <- as.factor(All$Term)
#levels(Endothelial$Term)
All$Term <- factor(All$Term, levels = All$Term)

# ggplot(Endothelial, aes(x=domain, y=Term, color = elimFisher, size = Enrichment)) + 
#   geom_point()+ scale_color_viridis(direction = -1, limits=c(0,0.1), na.value="#440154FF") # -log10(p.adjust)
png(paste0('plots/Amp_All_GO_enrichment_Protein.png'))
ggplot(All, aes(x=Enrichment, y=Term, color = elimFisher, size = Significant)) +  theme_readable() +
  geom_point()+ scale_color_viridis(direction = -1, limits=c(0,0.05), na.value="#440154FF") # -log10(p.adjust) # Enrichment = log2(Signifcant/Expected)
dev.off()
Neg <- Neg[which (Neg$Significant >= 10 & Neg$Enrichment >= 1),]
Neg <- Neg[order(Neg$Enrichment),]
dupl<-which(duplicated(Neg$Term))
#Neg[dupl[1],"Term"] <- paste0(Neg[dupl[1],"Term"],".")
#Neg[dupl[2],"Term"] <- paste0(Neg[dupl[2],"Term"],".")
Neg$Term <- as.factor(Neg$Term)
levels(Neg$Term)
Neg$Term <- factor(Neg$Term, levels = Neg$Term)
png(paste0('plots/Amp_Neg_GO_enrichment_Protein.png'),height = 900)
ggplot(Neg, aes(x=Enrichment, y=Term, color = elimFisher, size = Significant)) +  theme_readable() +
  geom_point()+ scale_color_viridis(direction = -1, limits=c(0,0.05), na.value="#440154FF") # -log10(p.adjust) # Enrichment = log2(Signifcant/Expected)
dev.off()
Pos <- Pos[which (Pos$Significant >= 10 & Pos$Enrichment >= 1),]
Pos <- Pos[order(Pos$Enrichment),]
Pos$Term <- as.factor(Pos$Term)
#levels(Endothelial$Term)
Pos$Term <- factor(Pos$Term, levels = Pos$Term)
png(paste0('plots/Amp_Pos_GO_enrichment_Protein.png'),height = 700)
ggplot(Pos, aes(x=Enrichment, y=Term, color = elimFisher, size = Significant)) + theme_readable() +
  geom_point()+ scale_color_viridis(direction = -1, limits=c(0,0.05), na.value="#440154FF") # -log10(p.adjust) # Enrichment = log2(Signifcant/Expected)
dev.off()



#### correlation RNA and protein 

all.equal(rownames(RNA_data),rownames(Protein_data))
all.equal(colnames(RNA_data),colnames(Protein_data))

cor(x=RNA_data[1,],y=Protein_data[1,])
cor(c(RNA_data),c(Protein_data), method = "spearman")

plot(c(RNA_data),c(Protein_data))
plot(RNA_data[2,],Protein_data[2,])

c(RNA_data[ribosomal_genes,]) # concanate columns
plot(c(RNA_data[ribosomal_genes,]),c(Protein_data[ribosomal_genes,]))
cor(c(RNA_data[ribosomal_genes,]),c(Protein_data[ribosomal_genes,]), method = "spearman")

plot(c(RNA_data[ribosomal_genes,1]),c(Protein_data[ribosomal_genes,1]))


# check whether these cell lines have rRNA and ribosomal proteins deleted
# ribosomal proteins
ribosomal_genes <- read.delim2("/cellfile/datapublic/rjohnen/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Correlation/ribosomal_genes.txt")
ribosomal_genes <- subset(ribosomal_genes,ribosomal_genes$Group.name %in% c('L ribosomal proteins', 'S ribosomal proteins'))$Approved.symbol

ribosomal_genes <- intersect(rownames(CNA_data), ribosomal_genes)

boxplot(CNA_data[ribosomal_genes,])

samples_without_del_in_ribosomes<- names(which(colSums(CNA_data[ribosomal_genes,] < 2) == 0 ))
samples_with_amp_in_ribosomes<- names(which(colSums(CNA_data[ribosomal_genes,] > 2) > 50 ))

pot_cell_lines <- intersect(samples_with_amp_in_ribosomes,samples_without_del_in_ribosomes)

boxplot(CNA_data[ribosomal_genes,samples_without_del_in_ribosomes])
boxplot(CNA_data[ribosomal_genes,pot_cell_lines])

CNA_table[samples_without_del_in_ribosomes,]

saveRDS(samples_without_del_in_ribosomes,"output/samples_without_del_in_ribosomes.rds")

saveRDS(pot_cell_lines,"output/pot_cell_lines.rds")
