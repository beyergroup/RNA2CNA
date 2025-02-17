source.dir_plot <- "../../../Analysis/"
source(paste0(source.dir_plot,"packages.R"))
source(paste0(source.dir_plot,"functions.R"))
library(viridis)
library(ComplexHeatmap)
library(reshape)
library(grid)
library(circlize)
library(ggpubr)


######### GO TERMS

All<- readRDS(file=paste0("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Correlation/Pancreas/output/spearman/Endo_All.rds"))
Neg <- readRDS(file=paste0("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Correlation/Pancreas/output/spearman/Endo_Neg.rds"))
Pos<- readRDS(file=paste0("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Correlation/Pancreas/output/spearman/Endo_Pos.rds")) 

######### Visualization
All <- All[which (All$Significant >= 10 & All$Enrichment >= 1),]
All <- All[order(All$Enrichment),]
All$Term <- as.factor(All$Term)
#levels(Endothelial$Term)
All$Term <- factor(All$Term, levels = All$Term)
Neg <- Neg[which (Neg$Significant >= 10 & Neg$Enrichment >= 1),]
Neg <- Neg[order(Neg$Enrichment),]
Neg$Term <- as.factor(Neg$Term)
#levels(Endothelial$Term)
Neg$Term <- factor(Neg$Term, levels = Neg$Term)

Pos <- Pos[which (Pos$Significant >= 10 & Pos$Enrichment >= 1),]
Pos <- Pos[order(Pos$Enrichment),]
Pos$Term <- as.factor(Pos$Term)
#levels(Endothelial$Term)
Pos$Term <- factor(Pos$Term, levels = Pos$Term)
Pos$Cor <- 'Positive'
Neg$Cor <- 'Negative'


df <- rbind(Pos,Neg)
df$Cor <- as.factor(df$Cor)
df$Cor <- factor(df$Cor,levels = c('Positive','Negative'))

pdf('plots/GO.pdf',height = 2.5, width = 7)
ggplot(df, aes(x=1, y=Term, color = elimFisher, size = Enrichment)) +  theme_manuscript() +
  geom_point()+ scale_color_viridis(direction = -1, limits=c(0,0.05), na.value="#440154FF", breaks=c(0,0.01,0.02,0.03,0.04,0.05)) + labs(size="log2\nenrichment", color="elim Fisher\np-value") +
  facet_wrap(.~Cor,scales='free') + theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(), legend.title = element_text(size=10), 
                                          legend.text = element_text(size=8),legend.position = 'right', axis.text.y = element_text(size=8,face = 'bold'), legend.key.height = unit(3, 'mm')) + 
  scale_size_continuous(breaks = c(1.5,2.0,2.5), range = c(3,8)) # -log10(p.adjust) # Enrichment = log2(Signifcant/Expected)
dev.off()

######## heatmap
rownames(df) <- df$Term

# Define min and max circle sizes
min_val <- round(min(df$Enrichment))  # 1.222762
max_val <- round(max(df$Enrichment))  # 1.562936
min_size <- 3  # Smallest circle (in mm)
max_size <- 6  # Largest circle (in mm)

# Function to scale values
scale_size <- function(value, min_val, max_val, min_size, max_size) {
  ((value - min_val) / (max_val - min_val)) * (max_size - min_size) + min_size
}


# pos 
df_heatmap_pos <- subset(df,df$Cor =='Positive')
df_heatmap_pos <- df_heatmap_pos[,c('elimFisher','Enrichment')]
elimFisher <- df_heatmap_pos[,c('elimFisher'),drop=F]
df_heatmap_pos <- df_heatmap_pos[,c('Enrichment'),drop=F]

col_fun <- colorRamp2(seq(0.05, 0, length = 5), viridis(5), reverse = T)  
nm = rownames(df_heatmap_pos)

# Convert data frame to matrix
mat <- as.matrix(df_heatmap_pos)  

# Define row groups
row_groups <- c("Group 1", "Group 1", "Group 1", "Group 2")

# Define colors for groups
group_colors <- c("Group 1" = "#009E73", "Group 2" = "#CC79A7")

# Define row annotation (fixing the label issue)
row_anno <- rowAnnotation(
  Group = anno_block(
    gp = gpar(fill = c("#009E73","#CC79A7") ) ,  # Assign colors per group
    labels = c("Cell growth", "Proliferation"),  # Labels for groups
    labels_gp = gpar(col = "black", fontsize = 10),
  )
)

# Define heatmap with row annotation
ht_pos <- Heatmap(mat, column_title = "Positive", row_title = "Term", 
                  row_order =c(1,3,4,2),
                  cluster_rows = F, show_column_names = FALSE,
                  cluster_columns = FALSE, col = col_fun,
                  rect_gp = gpar(type = "none"), 
                  left_annotation = row_anno,  # Add row annotation
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.rect(x = x, y = y, width = width, height = height, 
                              gp = gpar(col = "grey", fill = NA))
                    if(i >= j) {
                      circle_size <- scale_size(mat[i, j], min_val, max_val, min_size, max_size)
                      grid.circle(x = x, y = y, 
                                  r = unit(circle_size, "mm") / 2,  
                                  gp = gpar(fill = col_fun(elimFisher[i, j]), col = NA))  
                    }
                  }, 
                  show_heatmap_legend = FALSE,km=2
)

# Draw the heatmap
draw(ht_pos)


#black        orange       skyblue   bluishgreen        yellow 
##     "#000000"     "#E69F00"     "#56B4E9"     "#009E73"     "#F0E442" 
##          blue    vermillion reddishpurple          gray 
##     "#0072B2"     "#D55E00"     "#CC79A7"     "#999999" 

# neg 
df_heatmap_neg <- subset(df,df$Cor =='Negative')
df_heatmap_neg <- df_heatmap_neg[,c('elimFisher','Enrichment')]
elimFisher <- df_heatmap_neg[,c('elimFisher'),drop=F]
df_heatmap_neg <- df_heatmap_neg[,c('Enrichment'),drop=F]

nm = rownames(df_heatmap_neg)

# Convert data frame to matrix
mat <- as.matrix(df_heatmap_neg)  
# Define row groups
row_groups <- c("Group 1", "Group 1", "Group 1", "Group 2")

# Define colors for groups
group_colors <- c("Group 1" = "#E69F00", "Group 2" = "#56B4E9")

# Define row annotation (fixing the label issue)
row_anno <- rowAnnotation(
  Group = anno_block(
    gp = gpar(fill = c("#E69F00","#56B4E9") ) ,  # Assign colors per group
    labels = c("Apoptosis", "Translation"),  # Labels for groups
    labels_gp = gpar(col = "black", fontsize = 10),
  )
)

# Define heatmap with row annotation
ht_neg <- Heatmap(mat, 'elim Fisher \np-value',column_title = "Negative", row_title = "Term", 
                  row_order =c(1,3,4,2),
                  cluster_rows = F, show_column_names = FALSE,
                  cluster_columns = FALSE, col = col_fun,
                  rect_gp = gpar(type = "none"), 
                  left_annotation = row_anno,  # Add row annotation
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.rect(x = x, y = y, width = width, height = height, 
                              gp = gpar(col = "grey", fill = NA))
                    if(i >= j) {
                      circle_size <- scale_size(mat[i, j], min_val, max_val, min_size, max_size)
                      grid.circle(x = x, y = y, 
                                  r = unit(circle_size, "mm") / 2,  
                                  gp = gpar(fill = col_fun(elimFisher[i, j]), col = NA))  
                    }
                  }, 
                  show_heatmap_legend = TRUE,km=2
)
draw(ht_neg)

# Define the size legend with scaled sizes
size_legend = Legend(
  at = c(min_val, (min_val + max_val) / 2, max_val), 
  labels = sprintf("%.0f", c(min_val, (min_val + max_val) / 2, max_val)), 
  title = "log2 \nenrichment", 
  legend_gp = gpar(fill = "black"),
  type = "points",
  grid_height = unit(10, "mm"),
  grid_width = unit(10, "mm"),
  pch = 21,
  size = unit(scale_size(c(min_val, (min_val + max_val) / 2, max_val), 
                         min_val, max_val, min_size, max_size) , "mm")  # Match heatmap scaling
)
draw(ht_neg, heatmap_legend_list = list(size_legend))

grob_pos = grid.grabExpr(draw(ht_pos)) 
grob_neg = grid.grabExpr(draw(ht_neg, heatmap_legend_list = list(size_legend))) 


pdf('plots/GO_comb.pdf',height = 3.5, width = 8)
ggarrange(grob_pos,grob_neg,widths = c(0.45,0.55))
dev.off()


################ Translation stress
CNA_type<-"Deletion"
Expr_stress_melt<- readRDS(paste0("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Correlation/Pancreas/output/Expr_stress_melt_",CNA_type))

Expr_stress_melt$group <- as.factor(Expr_stress_melt$group)
levels(Expr_stress_melt$group)
Expr_stress_melt$group <- factor(Expr_stress_melt$group,labels = c("Cells with \nonly deletions", "Cells without CNAs"))
Expr_stress_melt$group <- factor(Expr_stress_melt$group,levels = c( "Cells without CNAs","Cells with \nonly deletions"))


Expr_stress_melt_del <- subset(Expr_stress_melt, Expr_stress_melt$group == "Cells with \nonly deletions")

CNA_type<-"Amplification"
Expr_stress_melt<- readRDS(paste0("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Correlation/Pancreas/output/Expr_stress_melt_",CNA_type))
head(Expr_stress_melt)

Expr_stress_melt$group <- as.factor(Expr_stress_melt$group)
levels(Expr_stress_melt$group)
Expr_stress_melt$group <- factor(Expr_stress_melt$group,labels = c("Cells with \nonly amplifications", "Cells without CNAs"))
Expr_stress_melt$group <- factor(Expr_stress_melt$group,levels = c( "Cells without CNAs","Cells with \nonly amplifications"))

Expr_stress_melt <- rbind(Expr_stress_melt,Expr_stress_melt_del)

med1<-round(median(subset(Expr_stress_melt,Expr_stress_melt$group == "Cells without CNAs")$value),3)
med2<-round(median(subset(Expr_stress_melt,Expr_stress_melt$group == "Cells with \nonly amplifications")$value),3)
med3<-round(median(subset(Expr_stress_melt,Expr_stress_melt$group == "Cells with \nonly deletions")$value),3)

wilcox.test(value~group, subset(Expr_stress_melt,Expr_stress_melt$group %in% c( "Cells without CNAs","Cells with \nonly amplifications")) ) # p-value = 0.03237
wilcox.test(value~group, subset(Expr_stress_melt,Expr_stress_melt$group %in% c( "Cells without CNAs","Cells with \nonly deletions")) ) # p-value = 1.671e-05

text<-c("*",'***', paste0('median =\n',med1),paste0('median =\n',med2),paste0('median =\n',med3))
d=data.frame(x=c(0.99,0.99,1,2,3,1), y=c(24,26,26,24,26,24), vx=c(1.02,2.02,0,0,0,0), vy=c(0,0,-1,-1,-1,-1))

pdf(paste0('plots/Expresssion_translation_stress_genes_cells_with_CNAs_vs_without_CNAs.pdf'),width = 3.5,height = 3.5)
p1<-ggplot(Expr_stress_melt)+geom_violin(aes(x=group,y=value))+ylab('Translation stress \ngene expression')+xlab("")+ylim(c(-9,26))+
  annotate("text",y=24,x=1.5, label = text[1],size=3)+ 
  annotate("text",y=26,x=2, label = text[2],size=3)+
  annotate("text", y=20, x=1, label=text[3], size=3)+
  annotate("text", y=20, x=2, label=text[4], size=3)+
  annotate("text", y=20, x=3, label=text[5], size=3)+
  geom_segment(data=d, mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy), size=0.5, color="black") +theme_manuscript()+theme(axis.text = element_text(size=10),
                                                                                                                                                           axis.text.x = element_text(size=10,angle = 45,hjust = 1,vjust = 1)
                                                                                                                                                           ,axis.title = element_text(size=10),
                                                                                                                                                           axis.title.y = element_text(size=10),
                                                                                                                                                           title  = element_text(size=10,face='bold'))
print(p1)
dev.off()


####ribosome

CNA_type<-"Deletion"
Expr_ribosomal_melt<- readRDS(paste0("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Correlation/Pancreas/output/Expr_ribosomal_melt_",CNA_type))

Expr_ribosomal_melt$group <- as.factor(Expr_ribosomal_melt$group)
levels(Expr_ribosomal_melt$group)
Expr_ribosomal_melt$group <- factor(Expr_ribosomal_melt$group,labels =  c("Cells with \nonly deletions", "Cells without CNAs"))
Expr_ribosomal_melt$group <- factor(Expr_ribosomal_melt$group,levels = c( "Cells without CNAs","Cells with \nonly deletions"))

Expr_ribosomal_melt_del <- subset(Expr_ribosomal_melt, Expr_ribosomal_melt$group == "Cells with \nonly deletions")

CNA_type<-"Amplification"
Expr_ribosomal_melt<- readRDS(paste0("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Correlation/Pancreas/output/Expr_ribosomal_melt_",CNA_type))

Expr_ribosomal_melt$group <- as.factor(Expr_ribosomal_melt$group)
levels(Expr_ribosomal_melt$group)
Expr_ribosomal_melt$group <- factor(Expr_ribosomal_melt$group,labels =  c("Cells with \nonly amplifications", "Cells without CNAs"))
Expr_ribosomal_melt$group <- factor(Expr_ribosomal_melt$group,levels =  c( "Cells without CNAs","Cells with \nonly amplifications"))

Expr_ribosomal_melt <- rbind(Expr_ribosomal_melt,Expr_ribosomal_melt_del)

med1<-round(median(subset(Expr_ribosomal_melt,Expr_ribosomal_melt$group == "Cells without CNAs")$value),3)
med2<-round(median(subset(Expr_ribosomal_melt,Expr_ribosomal_melt$group == "Cells with \nonly amplifications")$value),3)
med3<-round(median(subset(Expr_ribosomal_melt,Expr_ribosomal_melt$group == "Cells with \nonly deletions")$value),3)

wilcox.test(value~group, subset(Expr_ribosomal_melt,Expr_ribosomal_melt$group %in% c( "Cells without CNAs","Cells with \nonly amplifications")) ) # p-value =  0.001097
wilcox.test(value~group, subset(Expr_ribosomal_melt,Expr_ribosomal_melt$group %in% c( "Cells without CNAs","Cells with \nonly deletions")) ) # p-value  < 2.2e-16

text<-c("**",'***', paste0('median =\n',med1),paste0('median =\n',med2),paste0('median =\n',med3))
d=data.frame(x=c(0.99,0.99,1,2,3,1), y=c(13,14,14,13,14,13), vx=c(1.02,2.02,0,0,0,0), vy=c(0,0,-0.5,-0.5,-0.5,-0.5))


pdf(paste0('plots/Expresssion_translation_ribosomal_genes_cells_with_CNAs_vs_without_CNAs.pdf'),width = 3.5,height = 3.5)
p2<-ggplot(Expr_ribosomal_melt)+geom_violin(aes(x=group,y=value))+ylab('Ribosomal protein \ngene expression')+xlab("")+ylim(c(-3.5,14))+
  annotate("text",y=13,x=1.5, label = text[1],size=3)+ 
  annotate("text",y=14,x=2, label = text[2],size=3)+
  annotate("text", y=11, x=1, label=text[3], size=3)+
  annotate("text", y=11, x=2, label=text[4], size=3)+
  annotate("text", y=11, x=3, label=text[5], size=3)+
  geom_segment(data=d, mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy), size=0.5, color="black") +theme_manuscript()+theme(axis.text = element_text(size=10),
                                                                                                                                                          axis.text.x = element_text(size=10,angle = 45,hjust = 1,vjust = 1)
                                                                                                                                                          ,axis.title = element_text(size=10),
                                                                                                                                                          axis.title.y = element_text(size=10),
                                                                                                                                                          title  = element_text(size=10,face='bold'))
print(p2)
dev.off()

##### CCLE

#### RNA


All<-readRDS( file=paste0("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Correlation/CCLE/output/All_samples_without_del_in_ribosomes.rds"))
Neg<-readRDS( file=paste0("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Correlation/CCLE/output/Neg_samples_without_del_in_ribosomes.rds") )
Pos<-readRDS(file=paste0("../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Correlation/CCLE/output/Pos_samples_without_del_in_ribosomes.rds")) 

Neg <- Neg[which (Neg$Significant >= 10 & Neg$Enrichment >= 1),]
Neg <- Neg[order(Neg$Enrichment,decreasing = T),]
Neg <- Neg[1:21,]
Neg$Term <- as.factor(Neg$Term)
pdf(paste0('plots/Neg_GO_enrichment_samples_without_del_in_ribosomes.pdf'),width = 3.5,height = 3.5)
p4<-ggplot(Neg, aes(x=1, y=Term, color = elimFisher, size = Enrichment)) +  theme_manuscript() +
  geom_point()+ scale_color_viridis(direction = -1, limits=c(0,0.05), na.value="#440154FF", breaks=c(0,0.02,0.04)) + labs(size="log2\nenrichment", color="elim Fisher\np-value") +
   theme(axis.title.y = element_blank(), axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(), legend.title = element_text(size=10), 
                                          legend.text = element_text(size=8),legend.position = 'bottom',    legend.justification = "left", 
         legend.direction = "horizontal", axis.text.y = element_text(size=8), legend.key.height  = unit(3, 'mm'),legend.key.width =  unit(3, 'mm')) +
  scale_size_continuous(breaks = c(1.2,1.6,2.0,2.4), range = c(1,3)) + guides(size=guide_legend(nrow=3,byrow=TRUE))
p5<-ggplot(Neg, aes(x=Enrichment, y=Term, color = elimFisher)) +  theme_manuscript() +
  geom_point()+ scale_color_viridis(direction = -1, limits=c(0,0.05), na.value="#440154FF", breaks=c(0,0.02,0.04)) + labs(size="log2\nenrichment", color="elim Fisher\np-value") +
  theme(axis.title.y = element_text(size=10), axis.text.x = element_text(size=10),axis.title.x = element_text(size=10),legend.title = element_text(size=10), 
        legend.text = element_text(size=8),legend.position = 'bottom',    legend.justification = "left", 
        legend.direction = "horizontal", axis.text.y = element_text(size=8), legend.key.height  = unit(3, 'mm'),legend.key.width =  unit(3, 'mm'))
print(p5)
dev.off()

#### Protein
PT_pval <- readRDS('../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Correlation/CCLE/output/Spearman_cor_Number_Amp_vs_protein_abundance_CCLE_samples_without_del_in_ribosomal_proteins.rds')
PT_pval$Ribo <- as.factor(PT_pval$Ribo)
levels(PT_pval$Ribo)
PT_pval$Ribo <- factor(PT_pval$Ribo,levels = c('Yes','No'))
levels(PT_pval$Ribo)

pdf(paste0('plots/Amp_Dens_coeff.pdf'),width = 3.5,height = 2)
p3<-ggplot(PT_pval)+geom_density(aes(x=coef,fill=Ribo),alpha=0.5, position="identity")+theme_manuscript(legend_position = 'bottom')+xlab('correlation coefficient') + scale_fill_manual(values = c('#0379ff5e','grey10'))+ 
  labs(fill='Ribosomal protein')
print(p3)
dev.off()
#colorBlindness::cvdPlot(p3)

wilcox.test(subset(PT_pval,PT_pval$Ribo=="Yes")$coef,subset(PT_pval,PT_pval$Ribo=="No")$coef) #  p-value < 2.2e-16

## Go enrichment -> no translation stuff because ribosomal genes are not significant, but I cuold perform the GO 
