source.dir_plot <- "../../../Analysis/"
source(paste0(source.dir_plot,"packages.R"))
source(paste0(source.dir_plot,"functions.R"))
library(viridis)

dir_data <- "../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Aging_Analysis_Segments/Pancreas/Combine_loop/output/FDR01/"

Baron2016_Number_CNA_metadata <- readRDS(paste0(dir_data,'Baron2016_Number_CNA_metadata.rds'))
Baron2016_Number_CNA_metadata$Dataset <- 'Baron'
Peng2019_Number_CNA_metadata <- readRDS(paste0(dir_data,'Peng2019_Number_CNA_metadata.rds'))
Peng2019_Number_CNA_metadata$Dataset <- 'Peng'
Peng2019_Number_CNA_metadata_Normal <- subset(Peng2019_Number_CNA_metadata, Disease == 'Normal')
Peng2019_Number_CNA_metadata_Normal <- Peng2019_Number_CNA_metadata_Normal[,colnames(Baron2016_Number_CNA_metadata)]
Muraro2016_Number_CNA_metadata <- readRDS(paste0(dir_data,'Muraro2016_Number_CNA_metadata.rds'))
Muraro2016_Number_CNA_metadata$Dataset <- 'Muraro'
Muraro2016_Number_CNA_metadata <- Muraro2016_Number_CNA_metadata[,colnames(Baron2016_Number_CNA_metadata)]

Number_CNA_metadata <- rbind(Muraro2016_Number_CNA_metadata,Baron2016_Number_CNA_metadata,Peng2019_Number_CNA_metadata_Normal)
Number_CNA_metadata$Age <-as.numeric(Number_CNA_metadata$Age)

# I have to give the same celltype names
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
levels(Number_CNA_metadata$Sex)
levels(Number_CNA_metadata$Celltype)
Number_CNA_metadata$Sex <- revalue(Number_CNA_metadata$Sex, c('m'='m','f'='f','Male'='m','Female'='f',
                                                              'F'='f','M'='m'))

# min 100 Zellen insgesamt und in mind.2  Datasets
celltypes<-names(table(Number_CNA_metadata$Celltype)[table(Number_CNA_metadata$Celltype) > 100])
table(Number_CNA_metadata$Celltype,Number_CNA_metadata$Dataset)[celltypes,]

Number_CNA_metadata <- subset(Number_CNA_metadata, Number_CNA_metadata$Celltype %in% c("Alpha","Beta","Ductal","Acinar","Delta","Endothelial")) 
# rm < 15 years
Number_CNA_metadata <- subset(Number_CNA_metadata, Number_CNA_metadata$Age > 15) 


############# Age vs Number detected genes 
####### add graph on the side 
library(colorBlindness)
library(ggpubr)
col_pal<-palette.colors(palette = "Okabe-Ito")
col_pal<-c("Baron"="#DC267F","Muraro"="#648FFF","Peng"="#FFB000") # 785EF0

Number_CNA_metadata$Dataset <- as.factor(Number_CNA_metadata$Dataset)
levels(Number_CNA_metadata$Dataset)
Number_CNA_metadata$Dataset <- factor(Number_CNA_metadata$Dataset, levels=c("Muraro","Peng","Baron"))
levels(Number_CNA_metadata$Dataset)

main <- ggplot(Number_CNA_metadata)+geom_point(aes(y=nFeature_RNA,x=Age,col=Dataset),alpha=0.5, size=0.5)+geom_smooth(aes(x=Age,y=nFeature_RNA),method="lm",color="black")+ylab("Number of detected genes")+
  theme_manuscript()+  theme(legend.position = "none") + scale_color_manual(values=as.character(c(col_pal["Muraro"],col_pal["Peng"],col_pal["Baron"])))+xlim(c(15,67)) + guides(color = FALSE, fill = FALSE) # ylim(c(1900,3500))+
Hist_age <- ggplot(Number_CNA_metadata)+geom_histogram(aes(x=Age,fill=Dataset,col=Dataset), alpha=0.5, position="identity")  + theme_manuscript() + theme(axis.title.x=element_blank(),
                                                                                                                                                       axis.text.x=element_blank(),
                                                                                                                                                       axis.ticks.x=element_blank(),
                                                                                                                                                       legend.position = "none") + ylab('cell count')+
  scale_color_manual(values=as.character(c(col_pal["Muraro"],col_pal["Peng"],col_pal["Baron"]))) + scale_fill_manual(values=as.character(c(col_pal["Muraro"],col_pal["Peng"],col_pal["Baron"])))+xlim(c(15,67)) +  guides(color = FALSE, fill = FALSE) 
Hist_genes <- ggplot(Number_CNA_metadata)+geom_histogram(aes(y=nFeature_RNA,col=Dataset,fill=Dataset), alpha=0.5, position="identity")+ xlab('cell count')+
  theme_manuscript() + theme(axis.title.y=element_blank(),
                           axis.text.y=element_blank(),
                           axis.ticks.y=element_blank(),legend.position = "none") +
  scale_color_manual(values=as.character(c(col_pal["Muraro"],col_pal["Peng"],col_pal["Baron"]))) + scale_fill_manual(values=as.character(c(col_pal["Muraro"],col_pal["Peng"],col_pal["Baron"]))) + guides(color = FALSE, fill = FALSE) #+ylim(c(1900,3500))



# Extract the legend. Returns a gtable
legend <- as_ggplot(cowplot::get_legend(ggplot(Number_CNA_metadata)+geom_point(aes(y=nFeature_RNA,x=Age,col=Dataset)) + scale_color_manual(values=as.character(c(col_pal["Muraro"],col_pal["Peng"],col_pal["Baron"]))))) + theme_manuscript() 
# cvdPlot(ggarrange(Hist_age, legend, main, Hist_genes,  ncol = 2, nrow = 2, font.label = list(size = 30),
#                   widths = c(0.5,0.2), heights = c(0.2,0.5), align = "hv"))

pdf(paste0('plots/RNA2CNA_CBS_Segments_Age_vs_nFeature.pdf'),width = 7,height = 4)
ggarrange(Hist_age, legend, main, Hist_genes,  ncol = 2, nrow = 2, 
          widths = c(0.6,0.4), heights = c(0.2,0.5), align = "hv")
dev.off()


main <- ggplot(Number_CNA_metadata)+geom_point(aes(y=nCount_RNA,x=Age,col=Dataset),alpha=0.5, size=0.5)+geom_smooth(aes(x=Age,y=nCount_RNA),method="lm",color="black")+ylab("Number of counts")+
  theme_manuscript()+  theme(legend.position = "none") + scale_color_manual(values=as.character(c(col_pal["Muraro"],col_pal["Peng"],col_pal["Baron"])))+xlim(c(15,67)) + guides(color = FALSE, fill = FALSE) # ylim(c(1900,3500))+
Hist_age <- ggplot(Number_CNA_metadata)+geom_histogram(aes(x=Age,fill=Dataset,col=Dataset), alpha=0.5, position="identity")  + theme_manuscript()+ theme(axis.title.x=element_blank(),
                                                                                                                                                         axis.text.x=element_blank(),
                                                                                                                                                         axis.ticks.x=element_blank(),
                                                                                                                                                         legend.position = "none") + 
  scale_color_manual(values=as.character(c(col_pal["Muraro"],col_pal["Peng"],col_pal["Baron"]))) + scale_fill_manual(values=as.character(c(col_pal["Muraro"],col_pal["Peng"],col_pal["Baron"])))+xlim(c(15,67)) +  guides(color = FALSE, fill = FALSE) 
Hist_genes <- ggplot(Number_CNA_metadata)+geom_histogram(aes(y=nCount_RNA,col=Dataset,fill=Dataset), alpha=0.5, position="identity")+
  theme_manuscript() + theme(axis.title.y=element_blank(),
                             axis.text.y =element_blank(),
                             axis.ticks.y=element_blank(),legend.position = "none",axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) +
  scale_color_manual(values=as.character(c(col_pal["Muraro"],col_pal["Peng"],col_pal["Baron"]))) + scale_fill_manual(values=as.character(c(col_pal["Muraro"],col_pal["Peng"],col_pal["Baron"]))) + guides(color = FALSE, fill = FALSE) #+ylim(c(1900,3500))


# Extract the legend. Returns a gtable
legend <- as_ggplot(cowplot::get_legend(ggplot(Number_CNA_metadata)+geom_point(aes(y=nCount_RNA,x=Age,col=Dataset)) + theme_manuscript() + scale_color_manual(values=as.character(c(col_pal["Muraro"],col_pal["Peng"],col_pal["Baron"]))))) 
# cvdPlot(ggarrange(Hist_age, legend, main, Hist_genes,  ncol = 2, nrow = 2, font.label = list(size = 30),
#                   widths = c(0.5,0.2), heights = c(0.2,0.5), align = "hv"))

pdf(paste0('plots/RNA2CNA_CBS_Segments_Age_vs_nCount.pdf'),width = 3.5,height = 3.5)
ggarrange(Hist_age, legend, main, Hist_genes,  ncol = 2, nrow = 2, 
          widths = c(0.6,0.4), heights = c(0.2,0.5), align = "hv")
dev.off()

###### nCNA vs lib size
main <- ggplot(Number_CNA_metadata)+geom_point(aes(y=Number_CNAs_per_cell,x=nCount_RNA,col=Dataset),alpha=0.5, size=0.5)+geom_smooth(aes(x=nCount_RNA,y=Number_CNAs_per_cell),method="lm",color="black")+ylab("Number of CNAs per cell")+xlab("Library size")+
  theme_manuscript()+  theme(legend.position = "none") + scale_color_manual(values=as.character(c(col_pal["Muraro"],col_pal["Peng"],col_pal["Baron"])))+ guides(color = FALSE, fill = FALSE) # ylim(c(1900,3500))+
Hist_genes<- ggplot(Number_CNA_metadata)+geom_histogram(aes(x=nCount_RNA,fill=Dataset,col=Dataset), alpha=0.5, position="identity")  + theme_manuscript() + theme(axis.title.x=element_blank(),
                                                                                                                                                                    axis.text.x=element_blank(),
                                                                                                                                                                    axis.ticks.x=element_blank(),
                                                                                                                                                                    legend.position = "none") + ylab("Cell count")+
  scale_color_manual(values=as.character(c(col_pal["Muraro"],col_pal["Peng"],col_pal["Baron"]))) + scale_fill_manual(values=as.character(c(col_pal["Muraro"],col_pal["Peng"],col_pal["Baron"])))+  guides(color = FALSE, fill = FALSE) 
Hist_CNA <- ggplot(Number_CNA_metadata)+geom_histogram(aes(y=Number_CNAs_per_cell,col=Dataset,fill=Dataset), alpha=0.5, position="identity")+
  theme_manuscript() + theme(axis.title.y=element_blank(),
                             axis.text.y=element_blank(),
                             axis.ticks.y=element_blank(),legend.position = "none") +xlab("Cell count")+
  scale_color_manual(values=as.character(c(col_pal["Muraro"],col_pal["Peng"],col_pal["Baron"]))) + scale_fill_manual(values=as.character(c(col_pal["Muraro"],col_pal["Peng"],col_pal["Baron"]))) + guides(color = FALSE, fill = FALSE) #+ylim(c(1900,3500))


# Extract the legend. Returns a gtable
legend <- as_ggplot(cowplot::get_legend(ggplot(Number_CNA_metadata)+geom_point(aes(y=Number_CNAs_per_cell,x=nCount_RNA,col=Dataset)) + scale_color_manual(values=as.character(c(col_pal["Muraro"],col_pal["Peng"],col_pal["Baron"]))))) + theme_manuscript() 

pdf(paste0('plots/RNA2CNA_CBS_Segments_nCNAs_vs_nCount.pdf'),width = 6,height = 5)
ggarrange(Hist_genes, legend, main, Hist_CNA,  ncol = 2, nrow = 2, 
          widths = c(0.6,0.4), heights = c(0.2,0.5), align = "hv")
dev.off()


########## per cell type


table_for_redmine_del_comb <- readRDS("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Aging_Analysis_Segments/Pancreas/Combine_loop/output/Downsampling/Prob/70/table_for_redmine_del_comb.rds")
table_for_redmine_amp_comb <- readRDS("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Aging_Analysis_Segments/Pancreas/Combine_loop/output/Downsampling/Prob/70/table_for_redmine_amp_comb.rds")

table_for_redmine_del_comb$CNA <- 'Deletion'
table_for_redmine_amp_comb$CNA <- 'Amplification'

df <- rbind(table_for_redmine_amp_comb, table_for_redmine_del_comb)
df$Celltype <- as.factor(df$Celltype)
levels(df$Celltype)
df$Celltype <- factor(df$Celltype,levels = c("Endothelial","Delta","Beta","Alpha","Ductal","Acinar","All"))
levels(df$Celltype)

df$FDR <- as.factor(df$FDR)
levels(df$FDR)
#df$FDR <- factor(df$FDR, labels = c("FDR=0.05","FDR=0.1","FDR=0.2","FDR=0.3"))

df$CNA <- as.factor(df$CNA)
levels(df$CNA)

mean_across_celltypes <- c(median(subset(df,df$FDR == "FDR02" & df$CNA == "Amplification")$Coefficient.age),median(subset(df,df$FDR == "FDR02" & df$CNA == "Deletion")$Coefficient.age))
mean_across_celltypes <- c(mean(subset(df,df$FDR == "FDR02" & df$Celltype %in%  celltypes& df$CNA == "Amplification")$Coefficient.age),
                           mean(subset(df,df$FDR == "FDR02" & df$Celltype %in% celltypes & df$CNA == "Deletion")$Coefficient.age))

# FDR 0.2
pdf('plots/Coefficient_Age_per_celltype_FDR02.pdf',width = 7,height = 3.5)
ggplot(subset(df,df$FDR == "FDR02"))+geom_point(aes(y=Celltype,x=Coefficient.age,col=CNA), size=1.5, position=position_dodge(width = 0.5)) +
  geom_errorbar(aes(y=Celltype,xmin=Coefficient.age-error, xmax=Coefficient.age+error,col=CNA), width=.2, position=position_dodge(width = 0.5)) + 
  geom_vline(xintercept = 0,linetype="dashed",col="black")+xlab("Age coefficient (CNA/year)")+labs(colour="") +
  scale_color_manual(values = c("#D16103","#005AB5"))  + theme_manuscript() + theme(axis.text.x = element_text(angle = 45, hjust=1),legend.position = "top") + xlim(c(-0.001,0.0054))+ #xlim(c(-0.005,0.0054))+
  geom_vline(xintercept = mean_across_celltypes,col=c("#D16103","#005AB5"))+
  geom_text(data = subset(df, Celltype == "Endothelial" & CNA == "Deletion"), 
            aes(y = as.numeric(Celltype) - 0.5, x = 0.0054, label = "***"), 
            vjust = -0.5, size = 5)+
  geom_text(data = subset(df, Celltype == "Endothelial" & CNA == "Amplification"), 
            aes(y = as.numeric(Celltype) - 0.75, x = 0.005, label = "**"), 
            vjust = -0.5, size = 5) 
dev.off()

# significance
table_for_redmine_del_comb # ***
table_for_redmine_amp_comb # **

round(subset(table_for_redmine_amp_comb, Celltype == 'Endothelial' & FDR == 'FDR02')$Coefficient.age,4)
round(subset(table_for_redmine_del_comb, Celltype == 'Endothelial' & FDR == 'FDR02')$Coefficient.age,4)

########## Per dataset
# FDR 0.2
table_for_redmine_del_comb <- readRDS("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Aging_Analysis_Segments/Pancreas/Combine_loop/output/Downsampling/Prob/70/table_for_redmine_del_comb_perdataset.rds")
table_for_redmine_amp_comb <- readRDS("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Aging_Analysis_Segments/Pancreas/Combine_loop/output/Downsampling/Prob/70/table_for_redmine_amp_comb_perdataset.rds")


table_for_redmine_del_comb$CNA <- 'Deletion'
table_for_redmine_amp_comb$CNA <- 'Amplification'

df <- rbind(table_for_redmine_amp_comb, table_for_redmine_del_comb)
df$Celltype <- as.factor(df$Celltype)
levels(df$Celltype)
df$Celltype <- factor(df$Celltype,levels = c("Endothelial","Delta","Beta","Alpha","Ductal","Acinar","All"))
levels(df$Celltype)

df$Dataset <- as.factor(df$Dataset)
levels(df$Dataset)

# compute median
median_per_dataset <- c(median(subset(df,df$CNA == "Amplification" & df$Dataset == "Baron")$`Coefficient age`,na.rm = T),median(subset(df,df$CNA == "Deletion" & df$Dataset == "Baron")$`Coefficient age`,na.rm = T),
  median(subset(df,df$CNA == "Amplification" & df$Dataset == "Muraro")$`Coefficient age`,na.rm = T),median(subset(df,df$CNA == "Deletion" & df$Dataset == "Muraro")$`Coefficient age`,na.rm = T),
  median(subset(df,df$CNA == "Amplification" & df$Dataset == "Peng")$`Coefficient age`,na.rm = T),median(subset(df,df$CNA == "Deletion" & df$Dataset == "Peng")$`Coefficient age`,na.rm = T))
median_per_dataset <- data.frame("Median"=median_per_dataset, "Dataset"=c("Baron","Baron","Muraro","Muraro","Peng","Peng"),"CNA"=rep(c("Amplification","Deletion"),3))


pdf(paste0('plots/Coefficient_Age_per_dataset_2.pdf'),width = 3.5,height = 7)
ggplot(df)+geom_point(aes(y=Celltype,x=`Coefficient age`,col=CNA), size=1.5, position=position_dodge(width = 0.5)) +
  geom_errorbar(aes(y=Celltype,xmin=`Coefficient age`-error, xmax=`Coefficient age`+error,col=CNA), width=.2, position=position_dodge(width = 0.5)) + 
  geom_vline(xintercept = 0,linetype="dashed",col="black")+xlab("Age coefficient")+labs(colour="") +
  scale_color_manual(values = c("#D16103","#005AB5"))  + theme_manuscript() + theme(axis.text.x = element_text(angle = 45, hjust=1),legend.position = "top") + xlim(c(-0.0109,0.0109))+ facet_grid(Dataset~.) +
  geom_vline(data=median_per_dataset, aes(xintercept = Median),col=rep(c("#D16103","#005AB5"),3))
dev.off()


