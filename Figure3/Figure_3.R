source.dir_plot <- "../../../Analysis/"
source(paste0(source.dir_plot,"packages.R"))
source(paste0(source.dir_plot,"functions.R"))
library(ggpubr)
library(viridis)

Healthy <- readRDS('../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Onkogenic_CNAs/output/FDR01/Baron_Muraro_TS_Peng/tbl_pval_comb.rds')
Cancer <- readRDS('../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Onkogenic_CNAs/output/cancer/FDR01/tbl_pval_comb.rds')

Cancer<-subset(Cancer,Cancer$Celltype %in% unique(Healthy$Celltype))
# Healthy$Celltype<- as.factor(Healthy$Celltype)
# Healthy$Celltype<- factor(Healthy$Celltype, labels=c('Acinar - healthy','Alpha - healthy', 'Beta - healthy', 'Ductal - healthy', 'Endothelial - healthy'))
# Cancer$Celltype<- as.factor(Cancer$Celltype)
# Cancer$Celltype<- factor(Cancer$Celltype, labels=c('Acinar - PDAC','Ductal - PDAC', 'Endothelial - PDAC'))

Cancer$Gene<-as.factor(Cancer$Gene)
Cancer$Gene<-factor(Cancer$Gene,labels = c('Oncogenes','TS'))
Healthy$Gene<-as.factor(Healthy$Gene)
Healthy$Gene<-factor(Healthy$Gene,labels = c('Oncogenes','TS'))

Healthy$Disease <- 'healthy'
Cancer$Disease <- 'cancer'

tbl_pval_comb<-rbind(Healthy,Cancer)
levels(tbl_pval_comb$Celltype)

# Add a column for the original OddsRatio values for labeling
tbl_pval_comb$OddsRatio<-log2(tbl_pval_comb$OddsRatio)

# tbl_pval_comb$Celltype<- factor(tbl_pval_comb$Celltype, levels=c('Acinar - healthy','Acinar - PDAC','Ductal - healthy','Ductal - PDAC', 
# #                                                                  'Endothelial - healthy','Endothelial - PDAC','Alpha - healthy', 'Beta - healthy'))
# text<-as.factor(tbl_pval_comb$FDR < 0.05)
# text<-revalue(text,c('TRUE'='*','FALSE'=''))

# use absolute values as label names
# y_breaks <-c(-0.3,0,0.3,0.6)
# y_labels <-sprintf("%.1f", c(2^-0.3,2^0,2^0.3,2^0.6))
  
df<-subset(tbl_pval_comb,Celltype %in% c('Acinar'))
df$Disease <- as.factor(df$Disease)
levels(df$Disease)
df$Disease <- factor(df$Disease, labels =c('PDAC', 'Healthy'))
text1<-as.factor(df$FDR < 0.05)
text1<-revalue(text1,c('TRUE'='*','FALSE'=''))
plot_acinar<-ggplot(df)+geom_bar(aes(y=OddsRatio,x=Gene,fill=FDR),stat='identity')+geom_hline(yintercept=0,linetype='dashed')+facet_grid(CNA~Disease)+
  xlab('') + geom_text(aes(y=OddsRatio,x=Gene, label = text1),size=7)+theme_manuscript()+ylim(c(-0.5,0.62))+
  scale_fill_viridis(direction = -1, limits=c(0,0.1), na.value="#440154FF")+ theme(axis.text.x = element_text(size=10, angle = 45, hjust=1), axis.title.y = element_text(size=10,face = "bold"),plot.title = element_text(size=15,face ="bold"),
                                                                                   strip.text.y = element_blank(), legend.title =  element_text(size=10))+ggtitle('Acinar')+ylab(c('log2(odds ratio)'))
# scale_y_continuous(breaks=y_breaks, labels = y_labels) # Use the original values for x-axis labels
  
Legende_FDR <- as_ggplot(cowplot::get_legend(plot_acinar))
plot_acinar <- plot_acinar + theme(legend.position = "none")


df<-subset(tbl_pval_comb,Celltype %in% c('Ductal'))
df$Disease <- as.factor(df$Disease)
levels(df$Disease)
df$Disease <- factor(df$Disease, labels =c('PDAC', 'Healthy'))
text2<-df$FDR < 0.05
text2[which(df$FDR < 0.05)] <- rep('*',length(which(df$FDR < 0.05)))
text2[which(df$FDR < 0.01)] <- rep('**',length(which(df$FDR < 0.01)))
text2[which(df$FDR < 0.001)] <- rep('***',length(which(df$FDR < 0.001)))
text2<-revalue(text2,c('FALSE'=''))
plot_ductal<-ggplot(df)+geom_bar(aes(y=OddsRatio,x=Gene,fill=FDR),stat='identity')+geom_hline(yintercept=0,linetype='dashed')+facet_grid(CNA~Disease)+
  xlab('') + geom_text(aes(y=OddsRatio,x=Gene, label = text2),size=7)+theme_manuscript()+ylim(c(-0.5,0.62))+
  scale_fill_viridis(direction = -1, limits=c(0,0.1), na.value="#440154FF")+ theme(axis.text.x = element_text(size=10, angle = 45, hjust=1), axis.title.y = element_blank(),plot.title = element_text(size=15,face ="bold"),
                                                                                   axis.text.y = element_blank(),axis.ticks.y = element_blank(),
                                                                                   strip.text.y = element_blank(), legend.position = "none")+ggtitle('Ductal')

# remove alpha and beta from the plot
#df<-subset(tbl_pval_comb,Celltype %in% c('Alpha','Beta') & Disease %in% 'healthy')
#text<-as.factor(df$FDR < 0.05)
#text<-revalue(text,c('TRUE'='*','FALSE'=''))
#plot_healthy_alpha_beta<-ggplot(df)+geom_bar(aes(y=OddsRatio,x=Gene,fill=FDR),stat='identity')+geom_hline(yintercept=0,linetype='dashed')+facet_grid(CNA~Celltype)+
#  xlab('') + geom_text(aes(y=OddsRatio,x=Gene, label = text),size=10)+theme_manuscript()+ylim(c(-0.5,0.62))+
#  scale_fill_viridis(direction = -1, limits=c(0,0.1), na.value="#440154FF")+ theme(axis.text.x = element_text(size=10, angle = 45, hjust=1),axis.title.y = element_text(size=10,face ="bold"),plot.title = element_text(size=15,face ="bold"),
#                                                                                   strip.text.y = element_blank(),legend.position = "none")+ggtitle('Healthy')+ylab(c('log2(odds ratio)'))


df<-subset(tbl_pval_comb,Celltype %in% c('Endothelial'))
df$Disease <- as.factor(df$Disease)
levels(df$Disease)
df$Disease <- factor(df$Disease, labels =c('PDAC', 'Healthy'))
text3<-df$FDR < 0.05
text3[which(df$FDR < 0.05)] <- rep('*',length(which(df$FDR < 0.05)))
text3[which(df$FDR < 0.01)] <- rep('**',length(which(df$FDR < 0.01)))
text3[which(df$FDR < 0.001)] <- rep('***',length(which(df$FDR < 0.001)))
text3<-revalue(text3,c('FALSE'=''))
plot_endo<-ggplot(df)+geom_bar(aes(y=OddsRatio,x=Gene,fill=FDR),stat='identity')+geom_hline(yintercept=0,linetype='dashed')+facet_grid(CNA~Disease)+
  xlab('') + geom_text(aes(y=OddsRatio,x=Gene, label = text3),size=7)+theme_manuscript()+ylim(c(-0.5,0.62))+
  scale_fill_viridis(direction = -1, limits=c(0,0.1), na.value="#440154FF")+ theme(axis.text.x = element_text(size=10, angle = 45, hjust=1),axis.title.y = element_blank(),axis.text.y = element_blank(),
                                                                                   axis.ticks.y = element_blank(),plot.title = element_text(size=15,face ="bold"), legend.position = "none")+ggtitle('Endothelial')


pdf(paste0('plots/OddsRatio_per_celltype.pdf'),width = 7, height = 3.5)
ggarrange(plot_acinar,plot_ductal,plot_endo,Legende_FDR,
          align = "hv",
          ncol=4, nrow=1, widths = c(0.45,0.35,0.405,0.12),heights = c(0.33,0.33,0.33,0.33))
dev.off()






######## Heatmap Endothelial
tbl_pval_comb <- readRDS(file="../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Onkogenic_CNAs/output/FDR01/tbl_pval_comb_B_M_P_TS_percell_endothelial.rds")

library(ComplexHeatmap)
library(reshape)

# # FDR
tbl_pval_comb_heatmap_FDR <- reshape(tbl_pval_comb[,c("Test","Cell","FDR")],v.names="FDR", timevar = 'Cell',direction = 'wide', idvar=c("Test"))
rownames(tbl_pval_comb_heatmap_FDR) <- tbl_pval_comb_heatmap_FDR$Test
tbl_pval_comb_heatmap_FDR <- tbl_pval_comb_heatmap_FDR[,-1]
colnames(tbl_pval_comb_heatmap_FDR) <- sapply(colnames(tbl_pval_comb_heatmap_FDR), function(cell) strsplit(cell,"\\.")[[1]][2])
tbl_pval_comb_heatmap_FDR<-t(as.matrix(tbl_pval_comb_heatmap_FDR))

# #OR
tbl_pval_comb_heatmap_OR <- reshape(tbl_pval_comb[,c("Test","Cell","OddsRatio")],v.names="OddsRatio", timevar = 'Cell',direction = 'wide', idvar=c("Test"))
rownames(tbl_pval_comb_heatmap_OR) <- tbl_pval_comb_heatmap_OR$Test
tbl_pval_comb_heatmap_OR <- tbl_pval_comb_heatmap_OR[,-1]
colnames(tbl_pval_comb_heatmap_OR) <- sapply(colnames(tbl_pval_comb_heatmap_OR), function(cell) strsplit(cell,"\\.")[[1]][2])
tbl_pval_comb_heatmap_OR<-as.matrix(tbl_pval_comb_heatmap_OR)
tbl_pval_comb_heatmap_OR<- log2(tbl_pval_comb_heatmap_OR+0.000001)
tbl_pval_comb_heatmap_OR[which(tbl_pval_comb_heatmap_OR == log2(0.000001))] <- NA
tbl_pval_comb_heatmap_OR <- t(tbl_pval_comb_heatmap_OR)

annot <- matrix(ifelse( (tbl_pval_comb_heatmap_OR > 0 & tbl_pval_comb_heatmap_FDR < 0.05) , "*", ""), nrow(tbl_pval_comb_heatmap_FDR))
annot[is.na(annot)] <- ""
rownames(annot) <- rownames(tbl_pval_comb_heatmap_FDR)
colnames(annot) <- colnames(tbl_pval_comb_heatmap_FDR)

#nCNA
tbl_pval_comb_heatmap <- reshape(tbl_pval_comb[,c("Test","Cell","nCNA")],v.names="nCNA", timevar = 'Cell',direction = 'wide', idvar=c("Test"))
rownames(tbl_pval_comb_heatmap) <- tbl_pval_comb_heatmap$Test
tbl_pval_comb_heatmap <- tbl_pval_comb_heatmap[,-1]
colnames(tbl_pval_comb_heatmap) <- sapply(colnames(tbl_pval_comb_heatmap), function(cell) strsplit(cell,"\\.")[[1]][2])
tbl_pval_comb_heatmap<-as.matrix(tbl_pval_comb_heatmap)
tbl_pval_comb_heatmap <- t(tbl_pval_comb_heatmap)
tbl_pval_comb_heatmap[,3:4]<-tbl_pval_comb_heatmap[,3:4]*-1

all.equal(colnames(annot),colnames(tbl_pval_comb_heatmap))
all.equal(rownames(annot),rownames(tbl_pval_comb_heatmap))

# dendogram
# row_order <- tbl_pval_comb_heatmap[order(tbl_pval_comb_heatmap[,2])]

dend <- tbl_pval_comb_heatmap
dend[is.na(dend)] <- 0

dend <- hclust(dist(dend))
split <- c(1,1,2,2)
ha<-HeatmapAnnotation(foo = anno_block(gp = gpar(fill =c("white","white")),
                                       labels = c('Amplifications','Deletions'),labels_gp = gpar(fontsize = 10)))
# ha2<-rowAnnotation(foo = anno_block(gp = gpar(fill =c("transparent",'transparent',"#E1BE6A",'transparent','#5D3A9B','transparent','transparent','transparent','transparent','transparent','transparent'),
#                                               col=rep('transparent',10)),
#                                        labels = c('','','protective-\nstate','','cancer-\nstate','','','','',''),
#                                     labels_gp = gpar(fontsize = 8)))
library(colorRamp2)
col_fun<-  colorRamp2(c(min(tbl_pval_comb_heatmap), 0, max(tbl_pval_comb_heatmap)), c("#0071b2", "white", "#FE6100"))

colnames(tbl_pval_comb_heatmap) <- c('Oncogenes','TS','Oncogenes','TS')

pdf("plots/Heatmap_nCNA_endothelial_per_cell.pdf", height = 5.25, width = 3.5)
p<-Heatmap(tbl_pval_comb_heatmap,cluster_rows = dend, cluster_columns =  F, name = "Number of \nCNAs",
        top_annotation = ha , column_title = "Endothelial cells", column_title_gp = gpar(fontsize = 13, fontface = "bold"), col = col_fun,column_names_rot = 45, column_names_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10),
        row_title = "Cells",show_row_names = FALSE, column_split = split, column_gap = unit(2, "mm"), width = unit(5, "cm"), 
        heatmap_legend_param = list(
          title_gp = gpar(fontsize = 10)
        ),
        cell_fun=function(j, i, x, y, width, height, fill) {
          grid.text(annot[i, j], x, y, gp = gpar(fontsize = 10))
        })
print(p)

# decorate_heatmap_body("log\noddsratio", {
#   grid.rect(gp = gpar(fill = "transparent",lwd=4,col='#5D3A9B'))
# }, slice =5,column_slice=1)
# decorate_heatmap_body("log\noddsratio", {
#   grid.rect(gp = gpar(fill = "transparent",lwd=4,col='#5D3A9B'))
# }, slice = 5,column_slice=2)
# decorate_heatmap_body("log\noddsratio", {
#   grid.rect(gp = gpar(fill = "transparent",lwd=4,col='#E1BE6A'))
# }, slice =3,column_slice=1)
# decorate_heatmap_body("log\noddsratio", {
#   grid.rect(gp = gpar(fill = "transparent",lwd=4,col='#E1BE6A'))
# }, slice = 3,column_slice=2)
dev.off()




######## Heatmap Ductal
tbl_pval_comb <- readRDS(file="../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Onkogenic_CNAs/output/FDR01/tbl_pval_comb_B_M_P_TS_percell_ductal.rds")

tbl_pval_comb <- subset(tbl_pval_comb, !(tbl_pval_comb$Donor %in% c("TSP1","TSP9")))
#tbl_pval_comb <- subset(tbl_pval_comb, tbl_pval_comb$Donor %in% c("TSP1","TSP9"))

library(ComplexHeatmap)
library(reshape)

# #FDR
tbl_pval_comb_heatmap_FDR <- reshape(tbl_pval_comb[,c("Test","Cell","FDR")],v.names="FDR", timevar = 'Cell',direction = 'wide', idvar=c("Test"))
rownames(tbl_pval_comb_heatmap_FDR) <- tbl_pval_comb_heatmap_FDR$Test
tbl_pval_comb_heatmap_FDR <- tbl_pval_comb_heatmap_FDR[,-1]
colnames(tbl_pval_comb_heatmap_FDR) <- sapply(colnames(tbl_pval_comb_heatmap_FDR), function(cell) strsplit(cell,"\\.")[[1]][2])
tbl_pval_comb_heatmap_FDR<-t(as.matrix(tbl_pval_comb_heatmap_FDR))

# #OR
tbl_pval_comb_heatmap_OR <- reshape(tbl_pval_comb[,c("Test","Cell","OddsRatio")],v.names="OddsRatio", timevar = 'Cell',direction = 'wide', idvar=c("Test"))
rownames(tbl_pval_comb_heatmap_OR) <- tbl_pval_comb_heatmap_OR$Test
tbl_pval_comb_heatmap_OR <- tbl_pval_comb_heatmap_OR[,-1]
colnames(tbl_pval_comb_heatmap_OR) <- sapply(colnames(tbl_pval_comb_heatmap_OR), function(cell) strsplit(cell,"\\.")[[1]][2])
tbl_pval_comb_heatmap_OR<-as.matrix(tbl_pval_comb_heatmap_OR)
tbl_pval_comb_heatmap_OR<- log2(tbl_pval_comb_heatmap_OR+0.000001)
tbl_pval_comb_heatmap_OR[which(tbl_pval_comb_heatmap_OR == log2(0.000001))] <- NA
tbl_pval_comb_heatmap_OR <- t(tbl_pval_comb_heatmap_OR)

annot <- matrix(ifelse( (tbl_pval_comb_heatmap_OR > 0 & tbl_pval_comb_heatmap_FDR < 0.05) , "*", ""), nrow(tbl_pval_comb_heatmap_FDR))
annot[is.na(annot)] <- ""
rownames(annot) <- rownames(tbl_pval_comb_heatmap_FDR)
colnames(annot) <- colnames(tbl_pval_comb_heatmap_FDR)


#nCNA
tbl_pval_comb_heatmap <- reshape(tbl_pval_comb[,c("Test","Cell","nCNA")],v.names="nCNA", timevar = 'Cell',direction = 'wide', idvar=c("Test"))
rownames(tbl_pval_comb_heatmap) <- tbl_pval_comb_heatmap$Test
tbl_pval_comb_heatmap <- tbl_pval_comb_heatmap[,-1]
colnames(tbl_pval_comb_heatmap) <- sapply(colnames(tbl_pval_comb_heatmap), function(cell) strsplit(cell,"\\.")[[1]][2])
tbl_pval_comb_heatmap<-as.matrix(tbl_pval_comb_heatmap)
tbl_pval_comb_heatmap <- t(tbl_pval_comb_heatmap)
tbl_pval_comb_heatmap[,3:4]<-tbl_pval_comb_heatmap[,3:4]*-1

all.equal(colnames(annot),colnames(tbl_pval_comb_heatmap))
all.equal(rownames(annot),rownames(tbl_pval_comb_heatmap))

dend <- tbl_pval_comb_heatmap
dend[is.na(dend)] <- 0

dend <- hclust(dist(dend))
split <- c(1,1,2,2)
ha<-HeatmapAnnotation(foo = anno_block(gp = gpar(fill =c("white","white")),
                                       labels = c('Amplifications','Deletions'),labels_gp = gpar(fontsize = 10)))
# ha2<-rowAnnotation(foo = anno_block(gp = gpar(fill =c("transparent",'transparent',"#E1BE6A",'transparent','#5D3A9B','transparent','transparent','transparent','transparent','transparent','transparent'),
#                                               col=rep('transparent',10)),
#                                        labels = c('','','protective-\nstate','','cancer-\nstate','','','','',''),
#                                     labels_gp = gpar(fontsize = 8)))
library(colorRamp2)
col_fun<-  colorRamp2(c(min(tbl_pval_comb_heatmap), 0, max(tbl_pval_comb_heatmap)), c("#0071b2", "white", "#FE6100"))

colnames(tbl_pval_comb_heatmap) <- c('Oncogenes','TS','Oncogenes','TS')

pdf("plots/Heatmap_nCNA_ductal_per_cell_healthy.pdf", height = 5.25, width = 3.5)
p<-Heatmap(tbl_pval_comb_heatmap,cluster_rows = dend, cluster_columns =  F, name = "Number of \nCNAs",
           top_annotation = ha , column_title = "Ductal cells", column_title_gp = gpar(fontsize = 13, fontface = "bold"), col = col_fun,column_names_rot = 45, column_names_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10),
           row_title = "Cells",show_row_names = FALSE, column_split = split, column_gap = unit(2, "mm"), width = unit(5, "cm"), row_gap = unit(0, "mm"), row_split = 10, 
           heatmap_legend_param = list(
             title_gp = gpar(fontsize = 10)
           ),
           cell_fun=function(j, i, x, y, width, height, fill) {
             grid.text(annot[i, j], x, y, gp = gpar(fontsize = 10))
           })
print(p)
# 
# decorate_heatmap_body("log\noddsratio", {
#   grid.rect(gp = gpar(fill = "transparent",lwd=4,col='#5D3A9B'))
# }, slice =5,column_slice=1)
# decorate_heatmap_body("log\noddsratio", {
#   grid.rect(gp = gpar(fill = "transparent",lwd=4,col='#5D3A9B'))
# }, slice = 5,column_slice=2)
# decorate_heatmap_body("log\noddsratio", {
#   grid.rect(gp = gpar(fill = "transparent",lwd=4,col='#E1BE6A'))
# }, slice =3,column_slice=1)
# decorate_heatmap_body("log\noddsratio", {
#   grid.rect(gp = gpar(fill = "transparent",lwd=4,col='#E1BE6A'))
# }, slice = 3,column_slice=2)
dev.off()

rm(list=ls());gc();mallinfo()
