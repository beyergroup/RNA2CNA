source.dir <- "../../../../"
source(paste0(source.dir,"packages.R"))
source(paste0(source.dir,"functions.R"))
source(paste0(source.dir,"Functions/Compute_Number_CNAs_per_Cell.R"))
library(ggpubr)

output_dir <- c("../Onkogenic_CNAs/output/TS/FDR01/")
Tissues <- c("Blood","Intestine","Marrow","Muscle","Node","Pancreas","Prostate","Skin","Spleen","Thymus","Tongue")
Tissue_anot <- c("Blood","Large\nIntestine","Bone\nMarrow","Muscle","Lymph\nNode","Pancreas","Prostate","Skin","Spleen","Thymus","Tongue")
names(Tissue_anot) <- Tissues

FDR<-0.1

Number_CNAs<-list()
p<-list()
for (tissue in Tissues){
  print(tissue)
  metadata<-readRDS(file=paste0("/cellfile/cellnet/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_CBS/Tabula_Sapiens/",tissue,"/output/scaledata/v2/genewise/FixedMovingWindow/data.rds"))
  metadata<-metadata@meta.data
  Segments <- readRDS(file=paste0("/cellfile/cellnet/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_CBS/Tabula_Sapiens/",tissue,"/output/scaledata/v2/genewise/FixedMovingWindow/CBS_pvalues.rds"))
  Segments$CNA<-NA
  Segments[which(Segments$seg.mean > 0),"CNA"] <- "Amp"
  Segments[which(Segments$seg.mean < 0),"CNA"] <- "Del"
  #rm non sig calls
  Segments[which(Segments$FDR > FDR),"CNA"] <- "Zero"
  
  Segments_split <- split(Segments, Segments$Cell)
  Cells <- unique(Segments$Cell)
  N_CNAs <- Compute_Number_CNAs_per_Cell(Segments_split,metadata,Cells = Cells)#,Sep=F
  N_CNAs$tissue <- tissue
  
  Number_CNAs[[tissue]] <- N_CNAs
  # labs<- c("Amplification", "Deletion")
  # names(labs) <- c("Amp", "Del")
  # p[[tissue]]<-ggplot(Number_CNAs[[tissue]])+geom_violin(aes(x=cell_ontology_class,y=Number_CNAs_per_cell,col=CNA), show.legend = FALSE)+ coord_flip()+ 
  #   facet_grid(.~CNA,labeller = as_labeller(labs))+
  #   xlab("Cell type")+ylab('Number of CNAs per cell') + theme_readable() +  
  #   scale_color_manual(values = c("#D16103","#293352"))
}
names(Number_CNAs)<- Tissue_anot

saveRDS(Number_CNAs,paste0("output/TS_Number_CNAs.rds"))
saveRDS(p,paste0("output/TS_p.rds"))


### select cells
Number_CNAs<-readRDS(paste0("output/TS_Number_CNAs_no_Sep.rds"))
Number_CNAs_comb<-do.call(rbind,Number_CNAs)

table(Number_CNAs_comb$tissue)
pdf("plots_manuscript/Number_detected_Genes_per_tissue.pdf", width = 10)
ggplot(Number_CNAs_comb)+geom_boxplot(aes(y=nFeature_RNA, x=tissue)) + xlab("Tissue")+ylab('Number of detected genes per cell') + theme_readable()+geom_hline(yintercept= quantile(Number_CNAs_comb$nFeature_RNA,0.4),linetype="dashed" )+
  geom_hline(yintercept= quantile(Number_CNAs_comb$nFeature_RNA,0.6) ,linetype="dashed" )
dev.off()

quantile(Number_CNAs_comb$nFeature_RNA)[3]
quantile(Number_CNAs_comb$nFeature_RNA)[4]

Number_CNAs_comb$Quantile <- NA
# Number_CNAs_comb[Number_CNAs_comb$nFeature_RNA <= quantile(Number_CNAs_comb$nFeature_RNA)[2],"Quantile"] <- "Q1"
# Number_CNAs_comb[(Number_CNAs_comb$nFeature_RNA > quantile(Number_CNAs_comb$nFeature_RNA)[2] & Number_CNAs_comb$nFeature_RNA <= quantile(Number_CNAs_comb$nFeature_RNA)[3]),"Quantile"] <- "Q2"
Number_CNAs_comb[(Number_CNAs_comb$nFeature_RNA > quantile(Number_CNAs_comb$nFeature_RNA)[3] & Number_CNAs_comb$nFeature_RNA <= quantile(Number_CNAs_comb$nFeature_RNA)[4]),"Quantile"] <- "Q3"
# Number_CNAs_comb[(Number_CNAs_comb$nFeature_RNA > quantile(Number_CNAs_comb$nFeature_RNA)[4] & Number_CNAs_comb$nFeature_RNA <= quantile(Number_CNAs_comb$nFeature_RNA)[5]),"Quantile"] <- "Q4"
# Number_CNAs_comb[(Number_CNAs_comb$nFeature_RNA > quantile(Number_CNAs_comb$nFeature_RNA,0.4) & Number_CNAs_comb$nFeature_RNA <= quantile(Number_CNAs_comb$nFeature_RNA,0.6)),"Quantile"] <- "Q3"
cells_keep_50_75<- subset(Number_CNAs_comb,Quantile=="Q3")$Cell
length(cells_keep_50_75)
Number_CNAs_comb <- subset(Number_CNAs_comb,Number_CNAs_comb$Cell %in%  cells_keep_50_75)
table(Number_CNAs_comb$tissue)

saveRDS(cells_keep_50_75,'output/cells_keep_50_75.rds')

#compute % of cells with at least one CNA

perc_cells<-c()
for (t in Tissues){
  df<-subset(Number_CNAs_comb,Number_CNAs_comb$tissue == t)
  perc_cells <- append(perc_cells,(length(which(df$Number_CNAs_per_cell != 0))/length(which(df$Number_CNAs_per_cell == 0)))*100)
}
names(perc_cells) <- Tissues
perc_cells <- data.frame(perc_cells,'tissue'=names(perc_cells))

ggplot(perc_cells) + geom_bar(aes(y=perc_cells,x=tissue), stats='identitiy')

###

cells_keep_50_75<-readRDS('output/cells_keep_50_75.rds')
Number_CNAs<-readRDS(paste0("output/TS_Number_CNAs.rds"))
Number_CNAs_comb<-do.call(rbind,Number_CNAs)
Number_CNAs_comb <- subset(Number_CNAs_comb,Number_CNAs_comb$Cell %in%  cells_keep_50_75)

Number_CNAs_comb$CNA<-as.factor(Number_CNAs_comb$CNA)
levels(Number_CNAs_comb$CNA)<-c("Amplification",'Deletion')
Number_CNAs_comb$tissue<-as.factor(Number_CNAs_comb$tissue)
levels(Number_CNAs_comb$tissue)<-Tissue_anot

png("plots/Number_CNAs_Combined.png", width=600, height = 1500)
ggplot(Number_CNAs_comb)+geom_violin(aes(x=cell_ontology_class,y=Number_CNAs_per_cell,col=CNA), show.legend = FALSE)+ coord_flip()+ 
  facet_grid(tissue~CNA,scales = "free", space = "free")+
  xlab("Cell type")+ylab('Number of CNAs per cell') + theme_readable() +  
  scale_color_manual(values = c("#D16103","#293352"))
dev.off()

p<-ggplot(Number_CNAs_comb)+geom_violin(aes(x=cell_ontology_class,y=Number_CNAs_per_cell,col=CNA), show.legend = FALSE)+ coord_flip()+ 
  facet_grid(tissue~CNA,scales = "free", space = "free")+
  xlab("Cell type")+ylab('Number of CNAs per cell') + theme_readable() +  
  scale_color_manual(values = c("#D16103","#293352")) + theme(axis.text.y = element_text(size = 10),
                                                              strip.text = element_text(size = 12))

pdf("plots_manuscript/Number_CNAs_Combined.pdf", width=8.5, height = 20)
print(p)
dev.off()

# per tissue
p<-ggplot(Number_CNAs_comb)+geom_violin(aes(x=tissue,y=Number_CNAs_per_cell,col=CNA), show.legend = FALSE)+ coord_flip()+ 
  facet_grid(.~CNA, space = "free")+
  xlab("Tissue")+ylab('Number of CNAs per cell') + theme_readable() +  
  scale_color_manual(values = c("#D16103","#293352")) + theme(axis.text.y = element_text(size = 10),
                                                              strip.text = element_text(size = 12))

pdf("plots_manuscript/Number_CNAs_Combined2.pdf")
print(p)
dev.off()


p<-ggplot(Number_CNAs_comb)+geom_violin(aes(x=CNA,y=Number_CNAs_per_cell,col=CNA), show.legend = FALSE)+
  ylab('Number of CNAs per cell') + theme_readable() +  
  scale_color_manual(values = c("#D16103","#293352")) + theme(axis.text.y = element_text(size = 10),
                                                              strip.text = element_text(size = 12))

pdf("plots_manuscript/Number_CNAs_Combined3.pdf")
print(p)
dev.off()



# pdf("/data/public/rjohnen/Ronja/Tumor_Evolution/Manuscript/Figures/Figure_2/Number_CNAs_Combined.pdf", width=8.5, height = 20)
# print(p)
# dev.off()

#### downsample so that I can compare between tissues and cell type


#group based on nFeatureRNA
cells_keep_30_60<-readRDS('output/cells_keep_30_60.rds')
Number_CNAs_comb<- subset(Number_CNAs_comb,Number_CNAs_comb$Cell %in% cells_keep_30_60)

pdf("plots_manuscript/Number_detected_Genes_per_tissue_subset.pdf", width = 10)
ggplot(Number_CNAs_comb)+geom_boxplot(aes(y=nFeature_RNA, x=tissue)) + xlab("Tissue")+ylab('Number of detected genes per cell') + theme_readable()
dev.off()

pdf("plots_manuscript/Number_detected_Genes_quantiles.pdf", width = 10)
ggplot(Number_CNAs_comb)+geom_boxplot(aes(y=nFeature_RNA, x=Quantile)) + xlab("Quantiles")+ylab('Number of detected genes per cell') + theme_readable()
dev.off()


p<-ggplot(Number_CNAs_comb)+geom_violin(aes(x=tissue,y=Number_CNAs_per_cell,col=CNA), show.legend = FALSE)+ coord_flip()+ 
  facet_grid(Quantile~CNA, space = "free")+
  xlab("Tissue")+ylab('Number of CNAs per cell') + theme_readable() +  
  scale_color_manual(values = c("#D16103","#293352")) + theme(axis.text.y = element_text(size = 10),
                                                              strip.text = element_text(size = 12))

pdf("plots_manuscript/Number_CNAs_Combined_quantiles_per_tissue.pdf",height = 15)
print(p)
dev.off()


p<-ggplot(Number_CNAs_comb)+geom_violin(aes(x=Quantile,y=Number_CNAs_per_cell,col=CNA), show.legend = TRUE)+
  ylab('Number of CNAs per cell') + theme_readable() +  
  scale_color_manual(values = c("#D16103","#293352")) + theme(axis.text.y = element_text(size = 10),
                                                              strip.text = element_text(size = 12))

pdf("plots_manuscript/Number_CNAs_Combined_quantiles.pdf")
print(p)
dev.off()

ggplot(Number_CNAs_comb)+ geom_bar(aes(x=Quantile,fill=tissue))


pdf("plots_manuscript/Number_detected_Genes_per_tissue_Q3.pdf", width = 10)
ggplot(subset(Number_CNAs_comb,Quantile=="Q3"))+geom_boxplot(aes(y=nFeature_RNA, x=tissue)) + xlab("Tissue")+ylab('Number of detected genes per cell') 
dev.off()

p<-ggplot(subset(Number_CNAs_comb,Quantile=="Q3"))+geom_violin(aes(x=tissue,y=Number_CNAs_per_cell,col=CNA), show.legend = FALSE)+ coord_flip()+ 
  facet_grid(.~CNA, space = "free")+
  xlab("Tissue")+ylab('Number of CNAs per cell') + theme_readable() +  
  scale_color_manual(values = c("#D16103","#293352")) + theme(axis.text.y = element_text(size = 10),
                                                              strip.text = element_text(size = 12))


pdf("plots_manuscript/Number_CNAs_Combined_Q3.pdf")
print(p)
dev.off()

ggplot(subset(Number_CNAs_comb,Quantile=="Q3"))+ geom_bar(aes(x=Quantile,fill=tissue))

p<-ggplot(Number_CNAs_comb)+geom_violin(aes(x=CNA,y=Number_CNAs_per_cell,col=CNA), show.legend = FALSE)+
  ylab('Number of CNAs per cell') + theme_readable() +  
  scale_color_manual(values = c("#D16103","#293352")) + theme(axis.text.y = element_text(size = 10),
                                                              strip.text = element_text(size = 12))

pdf("plots_manuscript/Number_CNAs_Combined_Q3_Amp_vs_Del.pdf")
print(p)
dev.off()


wilcox.test(Number_CNAs_per_cell~CNA,Number_CNAs_comb) # 7.809e-12 Amp vs del 
summary(subset(Number_CNAs_comb,CNA=="Amplification")$Number_CNAs_per_cell)
summary(subset(Number_CNAs_comb,CNA=="Deletion")$Number_CNAs_per_cell)

#### pairwise tests
# test if we observe differences in number cnas per tissue 
kruskal.test(Number_CNAs_per_cell~tissue,Number_CNAs_comb)
pvalue_nCNA_tissue<-pairwise.wilcox.test(Number_CNAs_comb[,"Number_CNAs_per_cell"],Number_CNAs_comb[,"tissue"], p.adjust.method = 'fdr')

library(ComplexHeatmap)
Heatmap(pvalue_nCNA_tissue$p.value,  cluster_rows = FALSE, cluster_columns =  FALSE, name = "FDR",
        show_row_names = TRUE, col=rev(viridis(3)))


#### check if i have the same pvalues when I do the test by myself
Tissues <- levels(Number_CNAs_comb[,"tissue"])
pvalue<-matrix(nrow=length(Tissues), ncol=length(Tissues))
colnames(pvalue)<-Tissues
rownames(pvalue)<-Tissues
for (tissue1 in Tissues){
  for(tissue2 in Tissues){
    if(tissue1 != tissue2){
      print(tissue1)
      print(tissue2)
      #if(is.na(pvalue[tissue1,tissue2])){ 
        pvalue[tissue2,tissue1] <- wilcox.test(Number_CNAs_per_cell~tissue, subset(Number_CNAs_comb,tissue == tissue1 | tissue == tissue2 ))$p.value
      #}
    }
  }
}

# test<-matrix(c(pvalue),nrow=length(Tissues), ncol=length(Tissues))
# colnames(test)<-Tissues
# rownames(test)<-Tissues
# all.equal(test,pvalue)
padjusted<-matrix(p.adjust(pvalue,method = "fdr"),nrow=length(Tissues), ncol=length(Tissues))
colnames(padjusted)<-Tissues
rownames(padjusted)<-Tissues
#padjusted<-padjusted[-1,-ncol(padjusted)]
pvalue_nCNA_tissue$p.value
all.equal(padjusted,pvalue_nCNA_tissue$p.value)

Heatmap(padjusted,  cluster_rows = FALSE, cluster_columns =  FALSE, name = "FDR",
        show_row_names = TRUE, col=rev(viridis(3)))
# it is the same

###### create heatmap
annot <- matrix(ifelse(padjusted < 0.05, "*", ""), nrow(padjusted))
annot[is.na(annot)] <- ""

Tissues <- levels(Number_CNAs_comb[,"tissue"])
direction<-matrix(nrow=length(Tissues), ncol=length(Tissues))
colnames(direction)<-Tissues
rownames(direction)<-Tissues
for (tissue1 in Tissues){
  for(tissue2 in Tissues){
    if(tissue1 != tissue2){
      print(tissue1)
      print(tissue2)
      #if(is.na(direction[tissue1,tissue2])){ # to not have a symmetrix matrix
        direction[tissue2,tissue1] <-  mean(subset(Number_CNAs_comb, tissue == tissue1)$Number_CNAs_per_cell) - mean(subset(Number_CNAs_comb, tissue == tissue2)$Number_CNAs_per_cell)
      #}
    }
  }
}
#direction<-direction[-1,-ncol(direction)]

saveRDS(direction,file='output/direction_tissue.rds')
saveRDS(annot,file='output/annot_tissue.rds')

pdf("plots_manuscript/Number_CNAs_Combined_Q3_wilcoxon_tests.pdf")
Heatmap(direction,  cluster_rows = T, cluster_columns =  T, name = "Mean \ndifference",
        show_row_names = TRUE, 
        cell_fun=function(j, i, x, y, width, height, fill) {
          grid.text(annot[i, j], x, y, gp = gpar(fontsize = 20))
        })
dev.off()

# Heatmap(ifelse(direction, "less", "more"),  cluster_rows = FALSE, cluster_columns =  FALSE, name = "Direction",
#         show_row_names = TRUE, col=c("more"="red","less"="blue"), 
#         cell_fun=function(j, i, x, y, width, height, fill) {
#           grid.text(annot[i, j], x, y, gp = gpar(fontsize = 20))
#         })



############### BONE MARROW
library(RColorBrewer)

BM<-subset(Number_CNAs_comb,tissue=="Bone\nMarrow")
BM <- droplevels(BM)
levels(BM$cell_ontology_class)
BM$cell_ontology_class<-factor(BM$cell_ontology_class,
                                         levels = c("hematopoietic stem cell", "myeloid progenitor","erythroid progenitor","granulocyte",'cd24 neutrophil',"monocyte","macrophage","nk cell","cd8-positive, alpha-beta t cell","cd4-positive, alpha-beta t cell",'memory b cell',"plasma cell"))
table(BM$cell_ontology_class)
#col_celltype<- c(rocket(6),viridis(6)[2:6])
#col_celltype<- c(rocket(6),viridis(6)[2:6])
col_celltype<- c("black",rev(brewer.pal(6, "Blues")),rev(brewer.pal(7, "Oranges"))[1:6])
names(col_celltype) <- levels(BM$cell_ontology_class)
#show_col(col_celltype)
col_celltype


p<-ggplot(BM)+geom_violin(aes(x=cell_ontology_class,y=Number_CNAs_per_cell,fill=cell_ontology_class), show.legend = T)+
  ylab('\n\n\n\n\nNumber of CNAs per cell') + theme_readable() +  
  scale_fill_manual(values = col_celltype) + theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1),axis.text.y = element_text(size = 12),
                                                              strip.text = element_text(size = 12))+xlab("")+labs(fill="Cell type")
p
pdf("plots_manuscript/Number_CNAs_Combined_Q3_BM_per_celltype.pdf",width =8,height = 5)
print(p)
dev.off()
png("plots_manuscript/Number_CNAs_Combined_Q3_BM_per_celltype.png",width =800,height = 500)
print(p)
dev.off()
