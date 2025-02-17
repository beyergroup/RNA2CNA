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

df<-subset(tbl_pval_comb,Celltype %in% c('Acinar','Ductal','Endothelial'))
df$Disease <- as.factor(df$Disease)
df$Celltype <- as.factor(df$Celltype)
levels(df$Disease)
df$Disease <- factor(df$Disease, labels =c('PDAC', 'Healthy'))
df$sig<-as.factor(df$FDR < 0.05)
df$sig<-revalue(df$sig,c('FALSE'='n.s.', 'TRUE'='significant'))
# df$sig[which(df$FDR < 0.05)] <- rep('*',length(which(df$FDR < 0.05)))
# df$sig[which(df$FDR < 0.01)] <- rep('**',length(which(df$FDR < 0.01)))
# df$sig[which(df$FDR < 0.001)] <- rep('***',length(which(df$FDR < 0.001)))
# df$sig<-revalue(df$sig,c('FALSE'=''))
df <- na.omit(df)

pdf('plots/OddsRatio_per_celltype2.pdf',width = 7, height = 3.5)
ggplot(df)+geom_point(aes(y=Celltype,x=OddsRatio,col=CNA,shape=sig), size=2.5, position=position_dodge(width = 0.5)) + 
  geom_vline(xintercept = 0,linetype="dashed",col="black")+xlab('Log2(odds ratio)')+ labs(colour="") + labs(shape="") + 
  scale_color_manual(values = c("#D16103","#005AB5"))  + theme_manuscript() + theme(axis.text.x = element_text(angle = 45, hjust=1),legend.position = "top") + xlim(c(-0.75,0.75))+facet_grid(Gene~Disease)
# + 
#   geom_text(data = df, 
#             aes(y = as.numeric(Celltype) - 0.5, x = OddsRatio + 0.2, label = sig), 
#             vjust = -0.5, size = 5)
dev.off()

