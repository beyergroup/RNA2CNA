source.dir <- "../../../../../"
source(paste0(source.dir,"packages.R"))
source(paste0(source.dir,"functions.R"))
dir_data <- "../../Aging_Analysis_Segments/Pancreas/Combine_loop/output/FDR01/"
library(reshape2)
library(readxl)
TranslationStressGenes <- read_excel("/cellfile/datapublic/rjohnen/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Correlation/TranslationStressGenes.xlsx",
                                     sheet = "RSR", col_names = FALSE)



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


####### Load Number CNAs 

Baron2016_Number_CNA_metadata_Sep <- readRDS(paste0(dir_data,'Baron2016_Number_CNA_metadata_Amp_Del.rds'))
Baron2016_Number_CNA_metadata_Sep$Dataset <- 'Baron'
Baron2016_Number_CNA_metadata_Sep$Disease <- "no"
Baron2016_Number_CNA_metadata_Sep$Disease[which(Baron2016_Number_CNA_metadata_Sep$Age == 59)] <- "yes"
Peng2019_Number_CNA_metadata_Sep <- readRDS(paste0(dir_data,'Peng2019_Number_CNA_metadata_Amp_Del.rds'))
Peng2019_Number_CNA_metadata_Sep$Dataset <- 'Peng'
Peng2019_Number_CNA_metadata_Sep_Normal <- subset(Peng2019_Number_CNA_metadata_Sep, Disease == 'Normal')
Peng2019_Number_CNA_metadata_Sep_Normal <- Peng2019_Number_CNA_metadata_Sep_Normal[,colnames(Baron2016_Number_CNA_metadata_Sep)]
Peng2019_Number_CNA_metadata_Sep_Normal$Disease <- 'no'
Muraro2016_Number_CNA_metadata_Sep <- readRDS(paste0(dir_data,'Muraro2016_Number_CNA_metadata_Amp_Del.rds'))
Muraro2016_Number_CNA_metadata_Sep$Dataset <- 'Muraro'
Muraro2016_Number_CNA_metadata_Sep$Disease <- "no"
Muraro2016_Number_CNA_metadata_Sep <- Muraro2016_Number_CNA_metadata_Sep[,colnames(Baron2016_Number_CNA_metadata_Sep)]

TS_Pancreas_Number_CNA_metadata <- readRDS(paste0('../../Aging_Analysis_Segments/Pancreas/Combine_loop/output/TS/Pancreas/',strsplit(dir_data,"/")[[1]][7],'/Pancreas_Number_CNA_metadata_Amp_Del.rds'))
TS_Pancreas_Number_CNA_metadata$Dataset <- 'TS_Pancreas'
TS_Pancreas_Number_CNA_metadata$Disease <- "no"
colnames(TS_Pancreas_Number_CNA_metadata)[9] <- "Donor"
colnames(TS_Pancreas_Number_CNA_metadata)[13] <- "Celltype"
TS_Pancreas_Number_CNA_metadata <- TS_Pancreas_Number_CNA_metadata[,colnames(Baron2016_Number_CNA_metadata_Sep)]

Number_CNA_metadata_Sep <- rbind(Muraro2016_Number_CNA_metadata_Sep,Baron2016_Number_CNA_metadata_Sep,Peng2019_Number_CNA_metadata_Sep_Normal,TS_Pancreas_Number_CNA_metadata)
Number_CNA_metadata_Sep$Age <-as.numeric(Number_CNA_metadata_Sep$Age)

Number_CNA_metadata_Sep <- droplevels(Number_CNA_metadata_Sep)
levels(Number_CNA_metadata_Sep$Celltype)
Number_CNA_metadata_Sep$Celltype <- revalue(Number_CNA_metadata_Sep$Celltype, 
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
levels(Number_CNA_metadata_Sep$Sex)
Number_CNA_metadata_Sep$Sex <- revalue(Number_CNA_metadata_Sep$Sex, c('m'='m','f'='f','Male'='m','Female'='f',
                                                                      'F'='f','M'='m'))
Number_CNA_metadata_Sep <- subset(Number_CNA_metadata_Sep, Number_CNA_metadata_Sep$Celltype %in% c("Alpha","Beta","Ductal","Acinar","Delta","Endothelial","Macrophage")) 
# rm < 20 years
Number_CNA_metadata_Sep <- subset(Number_CNA_metadata_Sep, Number_CNA_metadata_Sep$Age > 15) 

#Number_CNA_metadata_Sep<-subset(Number_CNA_metadata_Sep,Disease == "no")
table(Number_CNA_metadata_Sep$CNA)
Number_CNA_metadata_Sep$CNA <- as.factor(Number_CNA_metadata_Sep$CNA)
Number_CNA_metadata_Sep$CNA <- revalue(Number_CNA_metadata_Sep$CNA, c("Amp"="Amplification","Del"="Deletion"))

celltype<-c("Alpha","Beta","Ductal","Acinar","Delta","Endothelial")

CNA_type<-"Amplification"
celltype <- "Endothelial"
cells_amp <- subset(Number_CNA_metadata_Sep,Number_CNA_metadata_Sep$CNA == 'Amplification')
cells_amp<- subset(cells_amp,cells_amp$Number_CNAs_per_cell !=0)
cells_del <- subset(Number_CNA_metadata_Sep, Number_CNA_metadata_Sep$CNA == 'Deletion')
cells_del<- subset(cells_del,cells_del$Number_CNAs_per_cell !=0)

cells_keep_amp<- setdiff(cells_amp$Cell, cells_del$Cell)
cells_keep_del<- setdiff(cells_del$Cell, cells_amp$Cell)

# subset endothelial cells
Endo <- subset(Number_CNA_metadata_Sep,CNA == CNA_type)
Endo <- subset(Endo,Endo$Celltype == celltype)
Endo<-droplevels(Endo)

table(Endo$Number_CNAs_per_cell)

# Percentage of cells with at least one CNA 
(sum(table(Endo$Number_CNAs_per_cell)[2:length(table(Endo$Number_CNAs_per_cell))])/table(Endo$Number_CNAs_per_cell)[1])*100
(sum(table(Number_CNA_metadata_Sep$Number_CNAs_per_cell)[2:length(table(Number_CNA_metadata_Sep$Number_CNAs_per_cell))])/table(Number_CNA_metadata_Sep$Number_CNAs_per_cell)[1])*100

#### separate cells based on CNAs
if(CNA_type == 'Amplification'){
  Cells_noCNA<- Endo[!(Endo$Cell %in% cells_keep_amp),]
  Cells_noCNA<- Cells_noCNA[Cells_noCNA$Number_CNAs_per_cell == 0,]$Cell
  Cells_CNA<- Endo[(Endo$Cell %in% cells_keep_amp),]$Cell
}
if(CNA_type == 'Deletion'){
  Cells_noCNA<- Endo[!(Endo$Cell %in% cells_keep_del),]
  Cells_noCNA<- Cells_noCNA[Cells_noCNA$Number_CNAs_per_cell == 0,]$Cell
  Cells_CNA<- Endo[(Endo$Cell %in% cells_keep_del),]$Cell
}

stress_genes <- intersect(rownames(Expr),TranslationStressGenes$...1)
Expr_stress <- Expr[stress_genes,]
#Expr_stress<-t(apply(Expr_stress,1, function(x) x-mean(x))) ### center per gene

Expr_stress_NoCNA <- Expr_stress[,Cells_noCNA]
Expr_stress_CNA <- Expr_stress[,Cells_CNA]
Expr_stress_NoCNA_melt <- melt(Expr_stress_NoCNA)
Expr_stress_CNA_melt <- melt(Expr_stress_CNA)
Expr_stress_CNA_melt$value <- as.numeric(Expr_stress_CNA_melt$value)
Expr_stress_NoCNA_melt$value <- as.numeric(Expr_stress_NoCNA_melt$value)

Expr_stress_NoCNA_melt$group<- 'NoCNA'
Expr_stress_CNA_melt$group<- 'CNA'
Expr_stress_melt<-rbind(Expr_stress_CNA_melt,Expr_stress_NoCNA_melt)
table(is.na(Expr_stress_melt$group))
Expr_stress_melt$value <- as.numeric(Expr_stress_melt$value)

png(paste0('plots/',CNA_type,'/Expresssion_stress_genes_cells_with_',CNA_type,'_vs_without_',CNA_type,'.png'))
ggplot(Expr_stress_melt)+geom_violin(aes(x=group,y=value))+theme_readable()+ylab('sctransformed gene expression')
dev.off()

plot(rowMeans(Expr_stress_NoCNA),rowMeans(Expr_stress_CNA))

saveRDS(Expr_stress_melt,paste0("output/Expr_stress_melt_",CNA_type))

# t.test(value~group, Expr_stress_melt) # p-value = 5.424e-10
# summary(subset(Expr_stress_melt,group=='CNA')$value) # 0.05873
# summary(subset(Expr_stress_melt,group=='NoCNA')$value) # -0.00492

# ks.test(Expr_stress_melt$value, 'pnorm')# keine Normalverteilung
# ks.test(subset(Expr_stress_melt,group=='NoCNA')$value, 'pnorm')# keine Normalverteilung
# ks.test(subset(Expr_stress_melt,group=='CNA')$value, 'pnorm')# keine Normalverteilung

# wilcoxon test da nicht normalverteilt
wilcox.test(value~group, Expr_stress_melt) # Deletion p-value = 1.272e-05 | Amplification  0.001118
# mit Mittelwert -> gepaarter test 
summary(subset(Expr_stress_melt,group=='CNA')$value) # Deletion -0.3369 | Amplification   -0.35142 
summary(subset(Expr_stress_melt,group=='NoCNA')$value) #  Deletion -0.366572  | Amplification  -0.366024

######## repeat per donor
Expr_stress_melt$Donor <- Endo[match(Expr_stress_melt$Var2,Endo$Cell),'Donor']
Expr_stress_melt$Age <- Endo[match(Expr_stress_melt$Var2,Endo$Cell),'Age']
table(is.na(Expr_stress_melt$Donor))
table(Expr_stress_melt$Donor,Expr_stress_melt$Age)
unique(Expr_stress_melt[order(Expr_stress_melt$Age),'Donor'])
levels(Expr_stress_melt$Donor)
Expr_stress_melt$Donor <- factor(Expr_stress_melt$Donor,levels = unique(Expr_stress_melt[order(Expr_stress_melt$Age),'Donor']))


donor<-levels(Expr_stress_melt$Donor)

pval<-c()
CNA_smaller<-c()
for(d in donor){
  df<-subset(Expr_stress_melt,Expr_stress_melt$Donor==d)
  if(length(unique(df$group)) == 2){
    pval<-append(pval, wilcox.test(value~group,df)$p.value) 
    print(d)
    print(mean(subset(df,group=='CNA')$value))
    print(mean(subset(df,group=='NoCNA')$value))
    CNA_smaller <- append(CNA_smaller, mean(subset(df,group=='CNA')$value)< mean(subset(df,group=='NoCNA')$value)) # TRUE if CNA is smaller
  }else{pval<-append(pval,NA);CNA_smaller<-append(CNA_smaller,NA)}
}
names(pval) <-donor
names(CNA_smaller) <-donor
pval < 0.05 
padjust<-p.adjust(pval,'fdr')
padjust < 0.05 # Donor1  N2      N3    TSP9  N5 Donor_2  N7      N8  N9 Donor_4    TSP1 N10 N11

# bei allen gleicher trend außer bei Donor1, Donor_2,  TSP1 (dort ist expr bei CNA höher als bei noCNA)

#color by trend
col<-c('#a6a6f7','#ff8a8a')
col <- ifelse(CNA_smaller, col[1], col[2])
# color by significant
for(d in donor){
  if(is.na(padjust[d])){
    col[d] == 'grey'
  }else if( (padjust < 0.05)[d] == TRUE ){
    if(col[d] == '#a6a6f7'){col[d] <-'blue'}else{col[d] <-'red'}
  }else{print('test')}
}
# Only colour strips in x-direction
strip <- strip_themed(background_x = elem_list_rect(fill =  col))

col_celltype<-c('#a6a6f7','blue','#ff8a8a','red')
names(col_celltype)<-c('down','sig.down','up','sig.up')
#add legend
gplot<-ggplot(data.frame("col"=col_celltype,"celltype"=names(col_celltype),"x"=1:length(col_celltype),"y"=1:length(col_celltype))) + 
  geom_point(aes(x, y, color = celltype)) + 
  scale_color_manual(values=as.character(col_celltype), breaks = names(col_celltype))+ theme(legend.title = element_blank(),legend.text = element_text(size=10))
legend <- get_legend(gplot)


# png('plots/Expresssion_translation_stress_genes_cells_with_CNAs_vs_without_CNAs_per_donor.png')
# ggplot(Expr_stress_melt)+geom_violin(aes(x=group,y=value))+theme_readable()+ylab('sctransformed gene expression')+facet_wrap(~Donor)
# dev.off()
# donor<-levels(Expr_stress_melt$Donor)

png(paste0('plots/',CNA_type,'/Expresssion_translation_stress_genes_cells_with_',CNA_type,'_vs_without_',CNA_type,'_per_donor.png'),width = 600,height = 600)
p<-ggplot(Expr_stress_melt,aes(x=group,y=value))+geom_violin() +theme_readable()+ylab('sctransformed gene expression')+xlab('')+
  facet_wrap2(~Donor, strip.position='top', strip = strip) 
print(ggarrange(p,legend,ncol = 2,widths = c(0.8,0.2)))
dev.off()
