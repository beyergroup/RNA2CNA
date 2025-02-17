source.dir <- "../../../Analysis/"
source(paste0(source.dir,"packages.R"))
source(paste0(source.dir,"functions.R"))

Tissues <- c("Blood","Intestine","Marrow","Muscle","Node","Pancreas","Prostate","Skin","Spleen","Thymus","Tongue")
Tissue_anot <- c("Blood","Large Intestine","Bone Marrow","Muscle","Lymph Node","Pancreas","Prostate","Skin","Spleen","Thymus","Tongue")
names(Tissue_anot) <- Tissues

cells_keep_50_75<-readRDS('../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Overview_Plots/output/cells_keep_50_75.rds')
Number_CNAs<-readRDS(paste0("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Overview_Plots/output/TS_Number_CNAs.rds"))
Number_CNAs_comb<-do.call(rbind,Number_CNAs)
Number_CNAs_comb <- subset(Number_CNAs_comb,Number_CNAs_comb$Cell %in%  cells_keep_50_75)

Number_CNAs_comb$CNA<-as.factor(Number_CNAs_comb$CNA)
levels(Number_CNAs_comb$CNA)<-c("Amplification",'Deletion')
Number_CNAs_comb$tissue<-as.factor(Number_CNAs_comb$tissue)
levels(Number_CNAs_comb$tissue)<-Tissue_anot[levels(Number_CNAs_comb$tissue)]

sum(Number_CNAs_comb$Number_CNAs_per_cell) # 2596 CNAs
length(unique(Number_CNAs_comb$Cell)) # 26196

table(is.na(Number_CNAs_comb$cell_ontology_class))
####### tissue

### Fraction of Cells 
# before calculating the fraction of cells, I need to combine ampl+del 
TS_Number_CNAs_no_Sep<-readRDS(paste0("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Overview_Plots/output/TS_Number_CNAs_no_Sep.rds"))
Number_CNAs_no_Sep_comb<-do.call(rbind,TS_Number_CNAs_no_Sep)
Number_CNAs_no_Sep_comb <- subset(Number_CNAs_no_Sep_comb,Number_CNAs_no_Sep_comb$Cell %in%  cells_keep_50_75)
Number_CNAs_no_Sep_comb$tissue<-as.factor(Number_CNAs_no_Sep_comb$tissue)
levels(Number_CNAs_no_Sep_comb$tissue)<-Tissue_anot[levels(Number_CNAs_no_Sep_comb$tissue)]

sort(table(Number_CNAs_no_Sep_comb$cell_ontology_class))

# significance (kruskal wallis test oder fisher test?) -> mehrere ungepaarte Gruppen (da Gewebe mit unterschiedlicher Anzahl an Zellen), nicht normalverteilt, Zielgroesse ist verhältnisskaliert, daher Kruskal-Wallis-Test
kruskal.test(Number_CNAs_no_Sep_comb$Number_CNAs_per_cell ~ Number_CNAs_no_Sep_comb$tissue) #  p-value < 2.2e-16
pairwise.wilcox.test(Number_CNAs_no_Sep_comb$Number_CNAs_per_cell, Number_CNAs_no_Sep_comb$tissue, p.adjust.method = 'BH') 

wilcox.test(Number_CNAs_per_cell ~ tissue, data=subset(Number_CNAs_no_Sep_comb, tissue %in% c("Bone Marrow","Blood"))) # p-value < 2.2e-16

# fisher test?

Number_CNAs_comb_split <- split(Number_CNAs_no_Sep_comb,Number_CNAs_no_Sep_comb$tissue)
lapply(Number_CNAs_comb_split, function(df) mean(df$Number_CNAs_per_cell))

fraction <- lapply(Number_CNAs_comb_split, function(df) ( length(which(df$Number_CNAs_per_cell != 0)) / length(df$Number_CNAs_per_cell) )*100 )
fraction<- do.call(cbind,fraction)
fraction<-as.vector(fraction)
names(fraction) <- names(Number_CNAs_comb_split)
  
levels(Number_CNAs_comb$tissue)
Number_CNAs_comb$tissue <- factor(Number_CNAs_comb$tissue,levels=names(sort(fraction,decreasing = T)))
levels(Number_CNAs_comb$tissue)
levels(Number_CNAs_comb$tissue)

d=data.frame(x=c(0.99,11,1), y=c(21,21,21), vx=c(10,0,0), vy=c(0,-1,-1))

fraction<-data.frame("tissue"=names(fraction),"Freq"=fraction)

# fraction_comb <- fraction
# 
# Number_CNAs_comb_all <- Number_CNAs_comb

p1<-ggplot(subset(Number_CNAs_comb,Number_CNAs_comb$Number_CNAs_per_cell !=0)) +
  geom_violin(aes(x=tissue,y=Number_CNAs_per_cell,fill=CNA,col=CNA)) + 
  xlab("Tissue") + ylab('Number of CNAs per cell\n(cells with CNA)')  + theme_manuscript() +  
  scale_color_manual(values = c("#D16103","#005AB5"),"") +
  scale_fill_manual(values = c("#D16103","#005AB5"),"") + 
  theme(axis.text.y.right  = element_blank(),axis.ticks.y.right = element_blank(), axis.title.y.right= element_blank(), axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),legend.position="top",legend.title.position = "top", legend.box = "vertical")+
  geom_point(data=fraction,aes(x=tissue,y=Freq/1.5,shape="Fraction of cells with CNA"),size=2) + 
  scale_y_continuous(
    # Features of the first axis
    name = "Number of CNAs per cell\n(cells with CNA)",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*1.5, name="") 
  ) + 
  guides(col=guide_legend(title='Number of CNAs per cell\n(cells with CNA)'),fill=guide_legend(title='Number of CNAs per cell\n(cells with CNA)'))+
  scale_shape_manual(name = "",
                     labels = c("Fraction of cells with CNA"),
                     values = c(5),
                     guide = guide_legend(theme = theme(legend.title=element_blank(),
                       legend.text = element_text(size = 10, face = "bold")))
                     )+
  annotate("text",y=21,x=6, label = "***",size=4) +
  geom_segment(data=d, mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy), size=0.5, color="black")


legend_Fig1cde <- get_legend(p1)
plot(legend_Fig1cde)

pdf('plots/Figure1cde_Legend.pdf',width=3.5,height =3.5)
plot(legend_Fig1cde)
dev.off()

pdf('plots/Figure1c_Number_CNAs_tissue.pdf',width=3.5,height = 2.8)
p1+theme(legend.position = "none")
dev.off()
library(colorBlindness)
cvdPlot(p1)

#####
# Supplement
library(viridis)
library(cowplot)
Number_CNAs_BM_Blood <- droplevels(subset(Number_CNAs_no_Sep_comb,Number_CNAs_no_Sep_comb$tissue %in% c("Bone Marrow","Blood"))) # Blood:  classical monocyte and  monocyte, BM: granulocyte
pdf('plots/Supplement/Number_celltype_BM_Blood.pdf',width=3.5,height = 3.5)
p1<-ggplot(Number_CNAs_BM_Blood) + geom_bar(aes(fill=cell_ontology_class,x=tissue),position = 'fill') + scale_fill_viridis(discrete = T) +theme_manuscript() + ylab("Cell fraction") + labs(fill="Cell type")+xlab("Tissue")
p1+theme(legend.position = "none")
dev.off()
pdf('plots/Supplement/Number_celltype_BM_Blood_legend.pdf',width=5,height = 5)
plot(get_legend(p1))
dev.off()

table(Number_CNAs_BM_Blood$cell_ontology_class,Number_CNAs_BM_Blood$tissue)

# BM
Number_CNAs_BM <- droplevels(subset(Number_CNAs_comb,Number_CNAs_comb$tissue %in% c("Bone Marrow")))
Number_CNAs_comb_split <- split(Number_CNAs_BM,Number_CNAs_BM$cell_ontology_class)
lapply(Number_CNAs_comb_split, function(df) mean(df$Number_CNAs_per_cell))

fraction <- lapply(Number_CNAs_comb_split, function(df) ( length(which(df$Number_CNAs_per_cell != 0)) / length(df$Number_CNAs_per_cell) )*100 )
fraction<- do.call(cbind,fraction)
fraction<-as.vector(fraction)
names(fraction) <- names(Number_CNAs_comb_split)

Number_CNAs_BM$cell_ontology_class <- factor(Number_CNAs_BM$cell_ontology_class,levels=names(sort(fraction,decreasing = T)))

fraction_BM<-data.frame("cell_ontology_class"=names(fraction),"Freq"=fraction)
p1<-ggplot(subset(Number_CNAs_BM,Number_CNAs_BM$Number_CNAs_per_cell !=0)) +
  geom_violin(aes(x=cell_ontology_class,y=Number_CNAs_per_cell,fill=CNA,col=CNA),drop=F) + facet_wrap(.~tissue, scales = "free_x") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ 
  xlab("Tissue") + ylab('Number of CNAs per cell\n(cells with CNA)')  + theme_manuscript() +
  scale_color_manual(values = c("#D16103","#005AB5"),"") +
  scale_fill_manual(values = c("#D16103","#005AB5"),"") + 
  theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),legend.position="top",legend.title.position = "top", legend.box = "vertical")+
  geom_point(data=fraction_BM,aes(x=cell_ontology_class,y=Freq/1.5,shape="Fraction of cells with CNA"),size=2) + 
  scale_y_continuous(
    # Features of the first axis
    name = "Number of CNAs per cell\n(cells with CNA)",limits=c(0,30),
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*1.5, name="Fraction of cells with CNA") 
  ) + 
  guides(col=guide_legend(title='Number of CNAs per cell\n(cells with CNA)'),fill=guide_legend(title='Number of CNAs per cell\n(cells with CNA)'))+
  scale_shape_manual(name = "",
                     labels = c("Fraction of cells with CNA"),
                     values = c(5),
                     guide = guide_legend(theme = theme(legend.title=element_blank(),
                                                        legend.text = element_text(size = 10, face = "bold")))
  )
pdf('plots/Supplement/Number_CNAs_BM.pdf',width=3.5,height = 3.5)
p1+theme(legend.position = "none")
dev.off()

# Blood
Number_CNAs_Blood <- droplevels(subset(Number_CNAs_comb,Number_CNAs_comb$tissue %in% c("Blood")))
Number_CNAs_comb_split <- split(Number_CNAs_Blood,Number_CNAs_Blood$cell_ontology_class)
lapply(Number_CNAs_comb_split, function(df) mean(df$Number_CNAs_per_cell))

fraction <- lapply(Number_CNAs_comb_split, function(df) ( length(which(df$Number_CNAs_per_cell != 0)) / length(df$Number_CNAs_per_cell) )*100 )
fraction<- do.call(cbind,fraction)
fraction<-as.vector(fraction)
names(fraction) <- names(Number_CNAs_comb_split)
Number_CNAs_Blood$cell_ontology_class <- factor(Number_CNAs_Blood$cell_ontology_class,levels=names(sort(fraction,decreasing = T)))

fraction_Blood<-data.frame("cell_ontology_class"=names(fraction),"Freq"=fraction)
p1<-ggplot(subset(Number_CNAs_Blood,Number_CNAs_Blood$Number_CNAs_per_cell !=0)) +
  geom_violin(aes(x=cell_ontology_class,y=Number_CNAs_per_cell,fill=CNA,col=CNA),drop=F) + facet_wrap(.~tissue, scales = "free_x") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ 
  xlab("Tissue") + ylab('Number of CNAs per cell\n(cells with CNA)')  + theme_manuscript() +  
  scale_color_manual(values = c("#D16103","#005AB5"),"") +
  scale_fill_manual(values = c("#D16103","#005AB5"),"") + 
  theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),legend.position="top",legend.title.position = "top", legend.box = "vertical")+
  geom_point(data=fraction_Blood,aes(x=cell_ontology_class,y=Freq/1.5,shape="Fraction of cells with CNA"),size=2) + 
  scale_y_continuous(
    # Features of the first axis
    name = "Number of CNAs per cell\n(cells with CNA)",limits=c(0,30),
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*1.5, name="Fraction of cells with CNA") 
  ) + 
  guides(col=guide_legend(title='Number of CNAs per cell\n(cells with CNA)'),fill=guide_legend(title='Number of CNAs per cell\n(cells with CNA)'))+
  scale_shape_manual(name = "",
                     labels = c("Fraction of cells with CNA"),
                     values = c(5),
                     guide = guide_legend(theme = theme(legend.title=element_blank(),
                                                        legend.text = element_text(size = 10, face = "bold")))
  )
pdf('plots/Supplement/Number_CNAs_Blood.pdf',width=3.5,height = 4.6)
p1+theme(legend.position = "none")
dev.off()
# 
# Number_CNAs_Blood_BM <- rbind(Number_CNAs_BM,Number_CNAs_Blood)
# fraction_BM$tissue <- "Bone Marrow"
# fraction_Blood$tissue <- "Blood"
# fraction_Blood_BM <- rbind(fraction_BM,fraction_Blood)
# 
# 
# p1<-ggplot(subset(Number_CNAs_Blood_BM,Number_CNAs_Blood_BM$Number_CNAs_per_cell !=0)) +
#   geom_violin(aes(x=cell_ontology_class,y=Number_CNAs_per_cell,fill=CNA,col=CNA),drop=F) + facet_wrap(.~tissue, scales = "free_x") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ 
#   xlab("Tissue") + ylab('Number of CNAs per cell\n(cells with CNA)')  + theme_manuscript() +  
#   scale_color_manual(values = c("#D16103","#005AB5"),"") +
#   scale_fill_manual(values = c("#D16103","#005AB5"),"") + 
#   theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),legend.position="top",legend.title.position = "top", legend.box = "vertical")+
#   geom_point(data=fraction_Blood_BM,aes(x=cell_ontology_class,y=Freq/1.5,shape="Fraction of cells with CNA"),size=2) + 
#   scale_y_continuous(
#     # Features of the first axis
#     name = "Number of CNAs per cell\n(cells with CNA)",
#     # Add a second axis and specify its features
#     sec.axis = sec_axis(~.*1.5, name="Fraction of cells with CNA") 
#   ) + 
#   guides(col=guide_legend(title='Number of CNAs per cell\n(cells with CNA)'),fill=guide_legend(title='Number of CNAs per cell\n(cells with CNA)'))+
#   scale_shape_manual(name = "",
#                      labels = c("Fraction of cells with CNA"),
#                      values = c(5),
#                      guide = guide_legend(theme = theme(legend.title=element_blank(),
#                                                         legend.text = element_text(size = 10, face = "bold")))
#   )
# 
# p1+theme(legend.position = "none")
#####

####### compartment

kruskal.test(Number_CNAs_no_Sep_comb$Number_CNAs_per_cell ~ Number_CNAs_no_Sep_comb$compartment) #  p-value < 2.2e-16
wilcox.test(Number_CNAs_per_cell ~ compartment, data=subset(Number_CNAs_no_Sep_comb, compartment %in% c("immune","stromal"))) # p-value < 2.2e-16
wilcox.test(Number_CNAs_per_cell ~ compartment, data=subset(Number_CNAs_no_Sep_comb, compartment %in% c("endothelial","stromal"))) # pvalue = 2.402e-09
pairwise.wilcox.test(Number_CNAs_no_Sep_comb$Number_CNAs_per_cell, Number_CNAs_no_Sep_comb$compartment, p.adjust.method = 'BH') 

Number_CNAs_comb_split <- split(Number_CNAs_no_Sep_comb,Number_CNAs_no_Sep_comb$compartment)
lapply(Number_CNAs_comb_split, function(df) mean(df$Number_CNAs_per_cell))
fraction <- lapply(Number_CNAs_comb_split, function(df) ( length(which(df$Number_CNAs_per_cell != 0)) / length(df$Number_CNAs_per_cell) )*100 )
fraction<- do.call(cbind,fraction)
fraction<-as.vector(fraction)
names(fraction) <- names(Number_CNAs_comb_split)
  
levels(Number_CNAs_comb$compartment)
Number_CNAs_comb$compartment <- factor(Number_CNAs_comb$compartment,levels=names(sort(fraction,decreasing = T)))
levels(Number_CNAs_comb$compartment)
levels(Number_CNAs_comb$compartment)

d=data.frame(x=c(0.99,3,1,1.99,3,2), y=c(21,21,21,19,19,19), vx=c(2,0,0,1,0,0), vy=c(0,-1,-1,0,-1,-1))

fraction<-data.frame("tissue"=names(fraction),"Freq"=fraction)
p2<-ggplot(subset(Number_CNAs_comb,Number_CNAs_comb$Number_CNAs_per_cell !=0))+geom_violin(aes(x=compartment,y=Number_CNAs_per_cell,fill=CNA,col=CNA))+
  xlab("Compartment")+ylab('') +  scale_color_manual(values = c("#D16103","#005AB5"),"")  + theme_manuscript() +
  scale_fill_manual(values = c("#D16103","#005AB5"),"") + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(),
                                                                axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),legend.position="top")+
  geom_point(data=fraction,aes(x=tissue,y=Freq/1.5,shape="Fraction of cells with CNA"),size=2) + 
   scale_y_continuous(
    
    # Features of the first axis
    name = "",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*1.5, name="")
  )+
  scale_shape_manual(name = "",
                     labels = c("Fraction of cells with CNA"),
                     values = c(5),
                     guide = guide_legend(theme = theme(legend.title=element_blank(),
                                                        legend.text = element_text(size = 10, face = "bold")))
  )+
  annotate("text",y=21,x=2, label = "***",size=4) +
  annotate("text",y=19,x=2.5, label = "***",size=4) +
  geom_segment(data=d, mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy), size=0.5, color="black")



pdf('plots/Figure1c_Number_CNAs_compartment.pdf',width=3.5,height =2.8)
p2+theme(legend.position = "none")
dev.off()
library(colorBlindness)
cvdPlot(p1)


###### stem cells
Number_CNAs_comb$Stem_cell <- "No"
Number_CNAs_comb[grep(pattern = "stem cell", x = Number_CNAs_comb$cell_ontology_class),"Stem_cell"]<-"Yes"

Number_CNAs_no_Sep_comb$Stem_cell <- "No"
Number_CNAs_no_Sep_comb[grep(pattern = "stem cell", x = Number_CNAs_no_Sep_comb$cell_ontology_class),"Stem_cell"]<-"Yes"

wilcox.test(Number_CNAs_no_Sep_comb$Number_CNAs_per_cell ~ Number_CNAs_no_Sep_comb$Stem_cell) # p-value < 2.2e-16

Number_CNAs_comb_split <- split(Number_CNAs_no_Sep_comb,Number_CNAs_no_Sep_comb$Stem_cell)
lapply(Number_CNAs_comb_split, function(df) mean(df$Number_CNAs_per_cell))

fraction <- lapply(Number_CNAs_comb_split, function(df) ( length(which(df$Number_CNAs_per_cell != 0)) / length(df$Number_CNAs_per_cell) )*100 )
fraction<- do.call(cbind,fraction)
fraction<-as.vector(fraction)
names(fraction) <- names(Number_CNAs_comb_split)

levels(Number_CNAs_comb$Stem_cell)
Number_CNAs_comb$Stem_cell <- factor(Number_CNAs_comb$Stem_cell,levels=names(sort(fraction,decreasing = T)))
levels(Number_CNAs_comb$Stem_cell)
levels(Number_CNAs_comb$Stem_cell)

wilcox.test(Number_CNAs_per_cell ~ Stem_cell, data=subset(Number_CNAs_comb,CNA == "Amplification")) # p-value < 2.2e-16
wilcox.test(Number_CNAs_per_cell ~ Stem_cell, data=subset(Number_CNAs_comb,CNA == "Deletion")) # p-value < 2.2e-16


d=data.frame(x=c(0.99,2,1), y=c(21,21,21), vx=c(1.02,0,0), vy=c(0,-1,-1))


fraction<-data.frame("tissue"=names(fraction),"Freq"=fraction)
p3<-ggplot(subset(Number_CNAs_comb,Number_CNAs_comb$Number_CNAs_per_cell !=0))+geom_violin(aes(x=Stem_cell,y=Number_CNAs_per_cell,fill=CNA,col=CNA))+
  xlab("Stem cell") +  scale_color_manual(values = c("#D16103","#005AB5"),"")  + theme_manuscript() +
  scale_fill_manual(values = c("#D16103","#005AB5"),"") + theme(axis.text.y.left = element_blank(),axis.ticks.y.left = element_blank(), axis.text.y.right  = element_text(size = 10),axis.title.y.left = element_blank(), 
                                                                axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position="top")+
  geom_point(data=fraction,aes(x=tissue,y=Freq/1.5,shape="Fraction of cells with CNA"),size=2) + 
   scale_y_continuous(
    
    # Features of the first axis
    name = "",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*1.5, name="Fraction of cells with CNA")
  )+
  scale_shape_manual(name = "",
                     labels = c("Fraction of cells with CNA"),
                     values = c(5),
                     guide = guide_legend(theme = theme(legend.title=element_blank(),
                                                        legend.text = element_text(size = 10, face = "bold")))
  )+
  annotate("text",y=21,x=1.5, label = "***",size=4) +
  geom_segment(data=d, mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy), size=0.5, color="black")
  


pdf('plots/Figure1c_Number_CNAs_Stem_cell.pdf',width=3.5,height =2.8)
p3+theme(legend.position = "none")
dev.off()
library(colorBlindness)
cvdPlot(p1)
rm();gc();mallinfo()

##### combine all plots 
pdf('plots/Figure1c_Number_CNAs_comb.pdf',width=7,height =2.8)
ggarrange((p1+theme(legend.position = "none")),(p2+theme(legend.position = "none")),(p3+theme(legend.position = "none")),ncol=3,nrow=1, align = 'h',widths = c(0.4,0.3,0.2))
dev.off()

#### per cell type

Number_CNAs_comb_cells_with_CNA<- subset(Number_CNAs_comb,Number_CNAs_comb$Number_CNAs_per_cell !=0)
Number_CNAs_comb_cells_with_CNA <- droplevels(Number_CNAs_comb_cells_with_CNA)

Tissues <-c("Blood","Large Intestine","Bone Marrow","Muscle","Lymph Node","Pancreas","Prostate","Skin","Spleen","Thymus","Tongue")
Tissue_anot <- c("Blood","Large Intestine","Bone Marrow","Muscle","Lymph Node","Pancreas","Prostate","Skin","Spleen","Thymus","Tongue")
names(Tissue_anot) <- Tissues
levels(Number_CNAs_comb_cells_with_CNA$tissue)<-Tissue_anot[levels(Number_CNAs_comb_cells_with_CNA$tissue)]

# abkürzungen für celltype

Number_CNAs_comb_cells_with_CNA$cell_ontology_class <- gsub("\\b-positive\\b", "⁺", Number_CNAs_comb_cells_with_CNA$cell_ontology_class)
Number_CNAs_comb_cells_with_CNA$cell_ontology_class <- gsub("\\b-negative\\b", "⁻", Number_CNAs_comb_cells_with_CNA$cell_ontology_class)
Number_CNAs_comb_cells_with_CNA$cell_ontology_class <- gsub("\\bderived\\b", "deri.", Number_CNAs_comb_cells_with_CNA$cell_ontology_class)
Number_CNAs_comb_cells_with_CNA$cell_ontology_class <- gsub("\\bassociated\\b", "asso.", Number_CNAs_comb_cells_with_CNA$cell_ontology_class)
Number_CNAs_comb_cells_with_CNA$cell_ontology_class <- gsub("\\balpha\\b", "\u03B1", Number_CNAs_comb_cells_with_CNA$cell_ontology_class)
Number_CNAs_comb_cells_with_CNA$cell_ontology_class <- gsub("\\bbeta\\b", "\u03B2", Number_CNAs_comb_cells_with_CNA$cell_ontology_class)
Number_CNAs_comb_cells_with_CNA$cell_ontology_class <- gsub("\\bvascular\\b", "vasc.", Number_CNAs_comb_cells_with_CNA$cell_ontology_class)
Number_CNAs_comb_cells_with_CNA$cell_ontology_class <- gsub("\\bsatellite\\b", "sat.", Number_CNAs_comb_cells_with_CNA$cell_ontology_class)
Number_CNAs_comb_cells_with_CNA$cell_ontology_class <- gsub("\\bprostate epithelium\\b", "prostate epithel.", Number_CNAs_comb_cells_with_CNA$cell_ontology_class)
Number_CNAs_comb_cells_with_CNA$cell_ontology_class <- gsub("\\bendothelial cell of", "endothel. cell of", Number_CNAs_comb_cells_with_CNA$cell_ontology_class)
Number_CNAs_comb_cells_with_CNA$cell_ontology_class <- gsub("\\b ,\\b", ",", Number_CNAs_comb_cells_with_CNA$cell_ontology_class)
Number_CNAs_comb_cells_with_CNA$cell_ontology_class <- gsub("\\bα-β\\b", "αβ", Number_CNAs_comb_cells_with_CNA$cell_ontology_class)
Number_CNAs_comb_cells_with_CNA$cell_ontology_class <- gsub("\\bnaive\\b", "naiv.", Number_CNAs_comb_cells_with_CNA$cell_ontology_class)
Number_CNAs_comb_cells_with_CNA$cell_ontology_class <- gsub("\\blymphatic\\b", "lymph.", Number_CNAs_comb_cells_with_CNA$cell_ontology_class)

blood <- subset(Number_CNAs_comb_cells_with_CNA,Number_CNAs_comb_cells_with_CNA$tissue %in% c('Blood'))
df1 <- subset(Number_CNAs_comb_cells_with_CNA,Number_CNAs_comb_cells_with_CNA$tissue %in% c("Large Intestine", "Bone Marrow","Muscle" ))
df2 <- subset(Number_CNAs_comb_cells_with_CNA,Number_CNAs_comb_cells_with_CNA$tissue %in% c("Lymph Node" ,"Pancreas", "Prostate" ,"Skin"))
df3 <- subset(Number_CNAs_comb_cells_with_CNA,Number_CNAs_comb_cells_with_CNA$tissue %in% c("Spleen" ,"Thymus" ,"Tongue"))
df1 <- droplevels(df1)
df2 <- droplevels(df2)
df3 <- droplevels(df3)
blood <- droplevels(blood)

table(df1$cell_ontology_class,df1$tissue) > 2

blood <- blood[!(blood$tissue == 'Blood' & blood$cell_ontology_class == 'plasma cell'),]
blood <- blood[!(blood$tissue == 'Blood' & blood$cell_ontology_class == 'cd8⁺, αβ t cell'),]
df1 <- df1[!(df1$tissue == 'Large Intestine' & df1$cell_ontology_class == 'plasma cell'),]
df1 <- df1[!(df1$tissue == 'Large Intestine' & df1$cell_ontology_class == 'paneth cell of epithelium of large intestine'),]
df1 <- df1[!(df1$tissue == 'Bone Marrow' & df1$cell_ontology_class == 'myeloid progenitor'),]
df1 <- df1[!(df1$tissue == 'Muscle' & df1$cell_ontology_class == 'tendon cell'),]
df1 <- df1[!(df1$tissue == 'Muscle' & df1$cell_ontology_class == 'endothel. cell of artery'),]


table(df2$cell_ontology_class,df2$tissue) > 2
df2 <- df2[!(df2$tissue == 'Prostate' & df2$cell_ontology_class == 'hillock cell of prostate epithel.'),]
df2 <- df2[!(df2$tissue == 'Prostate' & df2$cell_ontology_class == 'hillock-club cell of prostate epithel.'),]
df2 <- df2[!(df2$tissue == 'Lymph Node' & df2$cell_ontology_class == 'regulatory t cell'),]
df2 <- df2[!(df2$tissue == 'Skin' & df2$cell_ontology_class == 'cd8⁺, αβ memory t cell'),]
df3 <- df3[!(df3$tissue == 'Spleen' & df3$cell_ontology_class == 'innate lymphoid cell'),]
df3 <- df3[!(df3$tissue == 'Thymus' & df3$cell_ontology_class == 'naiv. regulatory t cell'),]
df3 <- df3[!(df3$tissue == 'Thymus' & df3$cell_ontology_class == 'medullary thymic epithelial cell'),]
df3 <- df3[!(df3$tissue == 'Thymus' & df3$cell_ontology_class == 'cd8⁺, αβ cytotoxic t cell'),]
df3 <- df3[!(df3$tissue == 'Thymus' & df3$cell_ontology_class == 'nk cell'),]
df3 <- df3[!(df3$tissue == 'Thymus' & df3$cell_ontology_class == 'memory b cell'),]
df3 <- df3[!(df3$tissue == 'Thymus' & df3$cell_ontology_class == 'cd8⁺, αβ t cell'),]


p_blood<-ggplot(blood) +
  geom_violin(aes(x=cell_ontology_class,y=Number_CNAs_per_cell,fill=CNA,col=CNA), show.legend = FALSE) + 
  xlab("") + ylab('Number of CNAs per cell\n(cells with CNA)')  + theme_manuscript() + coord_flip()+  
  scale_color_manual(values = c("#D16103","#005AB5"),"") +
  scale_fill_manual(values = c("#D16103","#005AB5"),"") + 
  theme(axis.text.y = element_text(size = 8),axis.title.x = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1), strip.text = element_text(size = 7), legend.position="top",legend.title.position = "top", legend.box = "vertical")+
  facet_grid(tissue~.,scales = "free", space = "free") + ylim(0,20)
p1<-ggplot(df1) +
  geom_violin(aes(x=cell_ontology_class,y=Number_CNAs_per_cell,fill=CNA,col=CNA), show.legend = FALSE) + 
  xlab("Cell type") + ylab('Number of CNAs per cell\n(cells with CNA)')  + theme_manuscript() + coord_flip()+  
  scale_color_manual(values = c("#D16103","#005AB5"),"") +
  scale_fill_manual(values = c("#D16103","#005AB5"),"") + 
  theme(axis.text.y = element_text(size = 8),axis.title.x = element_text(size = 8), 
        axis.text.x = element_text(size = 9, angle = 45, vjust = 1, hjust = 1), strip.text = element_text(size = 7), legend.position="top",legend.title.position = "top", legend.box = "vertical")+
  facet_grid(tissue~.,scales = "free", space = "free") + ylim(0,8)
p2<-ggplot(df2) +
  geom_violin(aes(x=cell_ontology_class,y=Number_CNAs_per_cell,fill=CNA,col=CNA), show.legend = FALSE) + 
  xlab("") + ylab('Number of CNAs per cell\n(cells with CNA)')  + theme_manuscript() + coord_flip()+  
  scale_color_manual(values = c("#D16103","#005AB5"),"") +
  scale_fill_manual(values = c("#D16103","#005AB5"),"") + 
  theme(axis.text.y = element_text(size = 8),axis.title.x = element_text(size = 8), 
        axis.text.x = element_text(size = 9, angle = 45, vjust = 1, hjust = 1), strip.text = element_text(size = 7), legend.position="top",legend.title.position = "top", legend.box = "vertical")+
  facet_grid(tissue~.,scales = "free", space = "free") + ylim(0,8)
p3<-ggplot(df3) +
  geom_violin(aes(x=cell_ontology_class,y=Number_CNAs_per_cell,fill=CNA,col=CNA), show.legend = FALSE) + 
  xlab("") + ylab('Number of CNAs per cell\n(cells with CNA)')  + theme_manuscript() + coord_flip()+  
  scale_color_manual(values = c("#D16103","#005AB5"),"") +
  scale_fill_manual(values = c("#D16103","#005AB5"),"") + 
  theme(axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 8), 
        axis.text.x = element_text(size = 9, angle = 45, vjust = 1, hjust = 1), strip.text = element_text(size = 7), legend.position="top",legend.title.position = "top", legend.box = "vertical")+
  facet_grid(tissue~.,scales = "free", space = "free") +  ylim(0,8)

library(ggpubr)
library(viridis)

cairo_pdf('plots/Figure1c_Number_CNAs_celltype.pdf',width=7.5, height = 5.7)
ggarrange(ggarrange(p1,p_blood, ncol =1, nrow = 2, align = "hv",heights = c(0.7,0.3)),
          p2,p3, ncol =3, nrow = 1,align = "hv",widths = c(0.325,0.325,0.35))
dev.off()


## fraction

Number_CNAs_no_Sep_comb <- droplevels(Number_CNAs_no_Sep_comb)
Number_CNAs_comb_split <- split(Number_CNAs_no_Sep_comb,Number_CNAs_no_Sep_comb$tissue)
Number_CNAs_comb_split <- lapply(Number_CNAs_comb_split, function(df) split(df, df$cell_ontology_class))

fraction <- lapply(Number_CNAs_comb_split, function(tissue) lapply(tissue, function(celltype) ( length(which(celltype$Number_CNAs_per_cell != 0)) / length(celltype$Number_CNAs_per_cell) )*100 )  )
fraction <- melt(fraction)
fraction[fraction == "NaN"] <- NA

colnames(fraction) <- c('Freq','cell_ontology_class', 'tissue')

# levels(Number_CNAs_comb$cell_ontology_class)
# Number_CNAs_comb$cell_ontology_class <- factor(Number_CNAs_comb$cell_ontology_class,levels=names(sort(fraction,decreasing = T)))
# levels(Number_CNAs_comb$cell_ontology_class)
# levels(Number_CNAs_comb$cell_ontology_class)

d=data.frame(x=c(0.99,11,1), y=c(21,21,21), vx=c(10,0,0), vy=c(0,-1,-1))

p2<-ggplot(subset(Number_CNAs_comb,Number_CNAs_comb$Number_CNAs_per_cell !=0)) +
  geom_violin(aes(x=cell_ontology_class,y=Number_CNAs_per_cell,fill=CNA,col=CNA), show.legend = FALSE) + 
  xlab("Cell type") + ylab('Number of CNAs per cell\n(cells with CNA)')  + theme_manuscript() + coord_flip()+  
  scale_color_manual(values = c("#D16103","#005AB5"),"") +
  scale_fill_manual(values = c("#D16103","#005AB5"),"") + 
  theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.text = element_text(size = 12), legend.position="top",legend.title.position = "top", legend.box = "vertical")+
  facet_grid(tissue~CNA, scales = "free", space = "free")+
  geom_point(data=fraction,aes(x=cell_ontology_class,y=Freq/1.5,shape="Fraction of cells with CNA"),size=2)+ 
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
  )

+
  annotate("text",y=21,x=6, label = "***",size=4) +
  geom_segment(data=d, mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy), size=0.5, color="black")

