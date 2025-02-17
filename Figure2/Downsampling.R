source.dir_plot <- "../../../Analysis/"
source(paste0(source.dir_plot,"packages.R"))
source(paste0(source.dir_plot,"functions.R"))
library(viridis)

dir_data <- "../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Aging_Analysis_Segments/Pancreas/Combine_loop/output/Downsampling/Prob/FDR02/"

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
col_pal<-c("pink"="#DC267F","blue"="#648FFF","orange"="#FFB000") # 785EF0

Number_CNA_metadata <- subset(Number_CNA_metadata,Celltype=="Endothelial")

main <- ggplot(Number_CNA_metadata)+geom_point(aes(y=nFeature_RNA,x=Age,col=Dataset),alpha=0.5, size=0.5)+geom_smooth(aes(x=Age,y=nFeature_RNA),method="lm",color="black")+ylab("Number of detected genes")+
  theme_manuscript()+  theme(legend.position = "none") + scale_color_manual(values=as.character(c(col_pal["pink"],col_pal["blue"],col_pal["orange"])))+xlim(c(15,67)) + guides(color = FALSE, fill = FALSE) # ylim(c(1900,3500))+
Hist_age <- ggplot(Number_CNA_metadata)+geom_histogram(aes(x=Age,fill=Dataset,col=Dataset), alpha=0.5, position="identity")  + theme_manuscript() + theme(axis.title.x=element_blank(),
                                                                                                                                                          axis.text.x=element_blank(),
                                                                                                                                                          axis.ticks.x=element_blank(),
                                                                                                                                                          legend.position = "none") + 
  scale_color_manual(values=as.character(c(col_pal["pink"],col_pal["blue"],col_pal["orange"]))) + scale_fill_manual(values=as.character(c(col_pal["pink"],col_pal["blue"],col_pal["orange"])))+xlim(c(15,67)) +  guides(color = FALSE, fill = FALSE) 
Hist_genes <- ggplot(Number_CNA_metadata)+geom_histogram(aes(y=nFeature_RNA,col=Dataset,fill=Dataset), alpha=0.5, position="identity")+
  theme_manuscript() + theme(axis.title.y=element_blank(),
                             axis.text.y=element_blank(),
                             axis.ticks.y=element_blank(),legend.position = "none") +
  scale_color_manual(values=as.character(c(col_pal["pink"],col_pal["blue"],col_pal["orange"]))) + scale_fill_manual(values=as.character(c(col_pal["pink"],col_pal["blue"],col_pal["orange"]))) + guides(color = FALSE, fill = FALSE) #+ylim(c(1900,3500))


# Extract the legend. Returns a gtable
legend <- as_ggplot(cowplot::get_legend(ggplot(Number_CNA_metadata)+geom_point(aes(y=nFeature_RNA,x=Age,col=Dataset)) + scale_color_manual(values=as.character(c(col_pal["pink"],col_pal["blue"],col_pal["orange"]))))) + theme_manuscript() 
# cvdPlot(ggarrange(Hist_age, legend, main, Hist_genes,  ncol = 2, nrow = 2, font.label = list(size = 30),
#                   widths = c(0.5,0.2), heights = c(0.2,0.5), align = "hv"))

pdf(paste0('plots/RNA2CNA_CBS_Segments_Age_vs_nFeature_downsampled_Endo.pdf'),width = 3.5,height = 3.5)
ggarrange(Hist_age, legend, main, Hist_genes,  ncol = 2, nrow = 2, 
          widths = c(0.6,0.4), heights = c(0.2,0.5), align = "hv")
dev.off()

summary(glm(nFeature_RNA~Age, data= Number_CNA_metadata,family=quasipoisson))
summary(glm(nFeature_RNA~Age, data= subset(Number_CNA_metadata,Celltype=="Endothelial"),family=quasipoisson))





main <- ggplot(Number_CNA_metadata)+geom_point(aes(y=Number_CNAs_per_cell,x=nFeature_RNA,col=Dataset),alpha=0.5, size=0.5)+geom_smooth(aes(x=nFeature_RNA,y=Number_CNAs_per_cell),method="lm",color="black")+ylab("Number of CNAs per cell")+xlab("Number of detected genes per cell")+
  theme_manuscript()+  theme(legend.position = "none") + scale_color_manual(values=as.character(c(col_pal["pink"],col_pal["blue"],col_pal["orange"])))+ guides(color = FALSE, fill = FALSE) # ylim(c(1900,3500))+
Hist_genes<- ggplot(Number_CNA_metadata)+geom_histogram(aes(x=nFeature_RNA,fill=Dataset,col=Dataset), alpha=0.5, position="identity")  + theme_manuscript() + theme(axis.title.x=element_blank(),
                                                                                                                                                          axis.text.x=element_blank(),
                                                                                                                                                          axis.ticks.x=element_blank(),
                                                                                                                                                          legend.position = "none") + 
  scale_color_manual(values=as.character(c(col_pal["pink"],col_pal["blue"],col_pal["orange"]))) + scale_fill_manual(values=as.character(c(col_pal["pink"],col_pal["blue"],col_pal["orange"])))+  guides(color = FALSE, fill = FALSE) 
Hist_CNA <- ggplot(Number_CNA_metadata)+geom_histogram(aes(y=Number_CNAs_per_cell,col=Dataset,fill=Dataset), alpha=0.5, position="identity")+
  theme_manuscript() + theme(axis.title.y=element_blank(),
                             axis.text.y=element_blank(),
                             axis.ticks.y=element_blank(),legend.position = "none") +
  scale_color_manual(values=as.character(c(col_pal["pink"],col_pal["blue"],col_pal["orange"]))) + scale_fill_manual(values=as.character(c(col_pal["pink"],col_pal["blue"],col_pal["orange"]))) + guides(color = FALSE, fill = FALSE) #+ylim(c(1900,3500))


# Extract the legend. Returns a gtable
legend <- as_ggplot(cowplot::get_legend(ggplot(Number_CNA_metadata)+geom_point(aes(y=Number_CNAs_per_cell,x=nFeature_RNA,col=Dataset)) + scale_color_manual(values=as.character(c(col_pal["pink"],col_pal["blue"],col_pal["orange"]))))) + theme_manuscript() 
# cvdPlot(ggarrange(Hist_age, legend, main, Hist_genes,  ncol = 2, nrow = 2, font.label = list(size = 30),
#                   widths = c(0.5,0.2), heights = c(0.2,0.5), align = "hv"))

pdf(paste0('plots/RNA2CNA_CBS_Segments_nCNAs_vs_nFeature_downsampled.pdf'),width = 3.5,height = 3.5)
ggarrange(Hist_genes, legend, main, Hist_CNA,  ncol = 2, nrow = 2, 
          widths = c(0.6,0.4), heights = c(0.2,0.5), align = "hv")
dev.off()

summary(glm(Number_CNAs_per_cell~nFeature_RNA, data= Number_CNA_metadata,family=quasipoisson))
summary(glm(Number_CNAs_per_cell~nFeature_RNA, data= subset(Number_CNA_metadata,Celltype=="Endothelial"),family=quasipoisson))

###### nCNA vs lib size
main <- ggplot(Number_CNA_metadata)+geom_point(aes(y=Number_CNAs_per_cell,x=nCount_RNA,col=Dataset),alpha=0.5, size=0.5)+geom_smooth(aes(x=nCount_RNA,y=Number_CNAs_per_cell),method="lm",color="black")+ylab("Number of CNAs per cell")+xlab("Library size")+
  theme_manuscript()+  theme(legend.position = "none") + scale_color_manual(values=as.character(c(col_pal["pink"],col_pal["blue"],col_pal["orange"])))+ guides(color = FALSE, fill = FALSE) # ylim(c(1900,3500))+
Hist_genes<- ggplot(Number_CNA_metadata)+geom_histogram(aes(x=nCount_RNA,fill=Dataset,col=Dataset), alpha=0.5, position="identity")  + theme_manuscript() + theme(axis.title.x=element_blank(),
                                                                                                                                                                  axis.text.x=element_blank(),
                                                                                                                                                                  axis.ticks.x=element_blank(),
                                                                                                                                                                  legend.position = "none") + ylab("Cell count")+
  scale_color_manual(values=as.character(c(col_pal["pink"],col_pal["blue"],col_pal["orange"]))) + scale_fill_manual(values=as.character(c(col_pal["pink"],col_pal["blue"],col_pal["orange"])))+  guides(color = FALSE, fill = FALSE) 
Hist_CNA <- ggplot(Number_CNA_metadata)+geom_histogram(aes(y=Number_CNAs_per_cell,col=Dataset,fill=Dataset), alpha=0.5, position="identity")+
  theme_manuscript() + theme(axis.title.y=element_blank(),
                             axis.text.y=element_blank(),
                             axis.ticks.y=element_blank(),legend.position = "none") +xlab("Cell count")+
  scale_color_manual(values=as.character(c(col_pal["pink"],col_pal["blue"],col_pal["orange"]))) + scale_fill_manual(values=as.character(c(col_pal["pink"],col_pal["blue"],col_pal["orange"]))) + guides(color = FALSE, fill = FALSE) #+ylim(c(1900,3500))


# Extract the legend. Returns a gtable
legend <- as_ggplot(cowplot::get_legend(ggplot(Number_CNA_metadata)+geom_point(aes(y=Number_CNAs_per_cell,x=nCount_RNA,col=Dataset)) + scale_color_manual(values=as.character(c(col_pal["pink"],col_pal["blue"],col_pal["orange"]))))) + theme_manuscript() 

pdf(paste0('plots/RNA2CNA_CBS_Segments_nCNAs_vs_nCount_downsampled.pdf'),width = 6,height = 5)
ggarrange(Hist_genes, legend, main, Hist_CNA,  ncol = 2, nrow = 2, 
          widths = c(0.6,0.4), heights = c(0.2,0.5), align = "hv")
dev.off()
summary(glm(Number_CNAs_per_cell~nCount_RNA, data= Number_CNA_metadata,family=quasipoisson))


###### Number of cells and CNA events per cell type
table(Number_CNA_metadata$Celltype)

ggplot(Number_CNA_metadata) + geom_violin(aes(x=Celltype,y=Number_CNAs_per_cell))

# significance (kruskal wallis test oder fisher test?) -> mehrere ungepaarte Gruppen (da Gewebe mit unterschiedlicher Anzahl an Zellen), nicht normalverteilt, Zielgroesse ist verhältnisskaliert, daher Kruskal-Wallis-Test
kruskal.test(Number_CNA_metadata$Number_CNAs_per_cell ~ Number_CNA_metadata$Celltype) #  p-value < 2.2e-16
pairwise.wilcox.test(Number_CNA_metadata$Number_CNAs_per_cell, Number_CNA_metadata$Celltype, p.adjust.method = 'BH')  # Acinar vs Endothelial p-value < 2.2e-16

Number_CNA_metadata<-droplevels(Number_CNA_metadata)
Number_CNAs_comb_split <- split(Number_CNA_metadata,Number_CNA_metadata$Celltype)
lapply(Number_CNAs_comb_split, function(df) mean(df$Number_CNAs_per_cell, na.rm = TRUE))

fraction <- lapply(Number_CNAs_comb_split, function(df) ( length(which(df$Number_CNAs_per_cell != 0)) / length(df$Number_CNAs_per_cell) )*100 )
fraction<- do.call(cbind,fraction)
fraction<-as.vector(fraction)
names(fraction) <- names(Number_CNAs_comb_split)

levels(Number_CNA_metadata$Celltype)
Number_CNA_metadata$Celltype <- factor(Number_CNA_metadata$Celltype,levels=names(sort(fraction,decreasing = T)))
levels(Number_CNA_metadata$Celltype)
levels(Number_CNA_metadata$Celltype)

d=data.frame(x=c(0.99,2,1), y=c(19,19,19), vx=c(1,0,0), vy=c(0,-1,-1))

fraction<-data.frame("Celltype"=names(fraction),"Freq"=fraction)
p1<-ggplot(subset(Number_CNA_metadata,Number_CNA_metadata$Number_CNAs_per_cell !=0)) +
  geom_violin(aes(x=Celltype,y=Number_CNAs_per_cell)) + 
  xlab("Cell type") + ylab('Number of CNAs per cell\n(cells with CNA)')  + theme_manuscript() +  
  theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),legend.position="top",legend.title.position = "top", legend.box = "vertical")+
  geom_point(data=fraction,aes(x=Celltype,y=Freq/1.5,shape="Fraction of cells with CNA"),size=2) + 
  scale_y_continuous(
    # Features of the first axis
    name = "Number of CNAs per cell\n(cells with CNA)",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*1.5, name="Fraction of cells with CNA") 
  ) + 
  guides(col=guide_legend(title='Number of CNAs per cell\n(cells with CNA)'),fill=guide_legend(title='Number of CNAs per cell\n(cells with CNA)'))+
  scale_shape_manual(name = "",
                     labels = c("Fraction of cells with CNA"),
                     values = c(5),
                     guide = guide_legend(theme = theme(legend.title=element_blank(),
                                                        legend.text = element_text(size = 10, face = "bold")))
  )+
  annotate("text",y=19,x=1.5, label = "***",size=4) +
  geom_segment(data=d, mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy), size=0.5, color="black")


legend_Fig1cde <- get_legend(p1)
plot(legend_Fig1cde)

p2 <- ggplot(Number_CNA_metadata) + geom_violin(aes(x=Celltype,y=nCount_RNA)) + theme_manuscript() + xlab('Cell type')+ylab('Library size')+  
  theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),legend.position="top",legend.title.position = "top", legend.box = "vertical")

pdf('plots/Supplement/Number_CNAs_celltype.pdf',width=7,height = 3.5)
ggarrange(p1,p2)
dev.off()
