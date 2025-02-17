source.dir <- "../../../../"
source(paste0(source.dir,"packages.R"))
source(paste0(source.dir,"functions.R"))
output_dir <- c("output/FDR01/")

### Load all Pancreas datasets
input_path <- c("Baron2016/output/scaledata/genewise/","Muraro2016/output/scaledata/genewise/",
                "Peng2019/output/scaledata/v2/genewise/healthy/") #"Segerstolpe2016/output/scaledata/genewise/","Xin2016/output/scaledata/genewise/","Lawlor2017/output/scaledata/genewise/","Enge2017/output/scaledata/genewise/",


min_cells <- 10


CNA_Matrix_comb <- readRDS(file='output/FDR01/CNA_Matrix_comb.rds')
metadata_comb <- readRDS(file='output/FDR01/metadata_comb.rds')
CNA_Matrix_comb <- CNA_Matrix_comb[,-1]
dim(CNA_Matrix_comb)
dim(metadata_comb)

# add TS Pancreas
TS_Pancreas <- readRDS(file=paste0('output/TS/FDR01/Pancreas_CBS_pvalues_FDR01_Matrix_Imputed.rds'))
data <- readRDS(file = paste0('/cellfile/cellnet/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_CBS/Tabula_Sapiens/Pancreas/output/scaledata/v2/genewise/FixedMovingWindow/data.rds'))
metadata <- data@meta.data
colnames(metadata)[6] <- "Donor"
colnames(metadata)[10] <- "Celltype"
CNA_Matrix_comb <- cbind(CNA_Matrix_comb,TS_Pancreas)
metadata_comb <- rbind(metadata_comb,metadata[,c('orig.ident', 'nCount_RNA', 'nFeature_RNA',  'Donor' , 'Celltype', 'n_counts_UMIs', 'n_genes', 'nCount_SCT', 'nFeature_SCT')])
##
dim(CNA_Matrix_comb)
dim(metadata_comb)

metadata_comb <- droplevels(metadata_comb)
levels(metadata_comb$Celltype)
metadata_comb$Celltype <- revalue(metadata_comb$Celltype, 
                                  c('alpha'='Alpha',"acinar"='Acinar','unsure' = 'Unclassified','delta'='Delta','beta'='Beta',
                                    'ductal'='Ductal','Ductal cell'='Ductal',
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
                                    "None/Other"='Unclassified',"Stellate"="Stellate","Ductal"="Ductal","Gamma/PP"="Gamma/PP",                                              "endothelial cell"="Endothelial", "myeloid cell"="Myeloid","pancreatic acinar cell"="Acinar",
                                    "pancreatic ductal cell"="Ductal", "pancreatic stellate cell"="Stellate", "t cell"="T cell"))

#metadata_comb <- subset(metadata_comb, metadata_comb$Celltype %in% c("Alpha","Beta","Ductal","Acinar","Delta","Endothelial","Macrophage")) 
levels(metadata_comb$Celltype)





metadata_comb$Age <- metadata_comb$Donor
levels(metadata_comb$Age)
metadata_comb$Age <- revalue(metadata_comb$Age, c('Donor1'=17,'Donor2'=51,'Donor3'=38,'Donor4'=59,
                                                  'Donor_1'=23,'Donor_2'=48,'Donor_3'=54,'Donor_4'=59,
                                                  "N1"="30", "N2"="31","N3"="34", "N4"="41","N5"="42", "N6"="50","N7"="52", "N8"="53","N9"="55", "N10"="64","N11"="65", 
                                                  "T1"="36","T2"="44", "T3"="51","T4"="52","T5"="54", "T6"="54","T7"="54","T8"="56", "T9"="58","T10"="58","T11"="59", "T12"="59",
                                                  "T13"="59","T14"="61", "T15"="64","T16"="64","T17"="65", "T18"="66","T19"="67","T20"="67", "T21"="68","T22"="70","T23"="71", 
                                                  "T24"="72",
                                                  'TSP1'=59,'TSP9'=37))



metadata_comb$Age <- as.numeric(as.character(metadata_comb$Age))
summary(metadata_comb$Age)
metadata_comb$AgeGroup <- NA
metadata_comb[which(metadata_comb$Age < median(metadata_comb$Age)),'AgeGroup'] <- 'Young'
metadata_comb[which(metadata_comb$Age >= median(metadata_comb$Age)),'AgeGroup'] <- 'Old'
table(metadata_comb$AgeGroup)


# per donor 
levels(metadata_comb$orig.ident)
levels(metadata_comb$Donor)
table(metadata_comb$orig.ident,metadata_comb$Donor)
metadata_comb$Dataset <- metadata_comb$Donor
metadata_comb$Dataset <- revalue(metadata_comb$Dataset, c('Donor1'='Baron','Donor2'='Baron','Donor3'='Baron','Donor4'='Baron',
                                                          'Donor_1'='Muraro','Donor_2'='Muraro','Donor_3'='Muraro','Donor_4'='Muraro',
                                                          'N1'='Peng','N2'='Peng','N3'='Peng','N4'='Peng','N5'='Peng','N6'='Peng',
                                                          'N7'='Peng','N8'='Peng','N9'='Peng','N10'='Peng','N11'='Peng','TSP1'='TS','TSP9'='TS'))
donors <- levels(metadata_comb$Donor)
table(metadata_comb$Celltype,metadata_comb$Dataset)

celltypes <- table(metadata_comb$Celltype,metadata_comb$Dataset)
celltypes <- rownames(celltypes)[apply(celltypes,1, function(row) sum(row) > 100 & length(which(row !=0)) >= 2)]
celltypes <- celltypes[ !(celltypes %in% c('Macrophage','Stellate','Delta')) ] # rm because after filter for #CNAs only one dataset

table(metadata_comb$Celltype,metadata_comb$Donor)



# Load cancer genes
TumorSuppressorGene_df<- readRDS('output/Pancreas_del_TS.rds')
Oncogenes_df<- readRDS('output/Pancreas_amp_Onkogenes.rds')


# include number of CNAs for filter (not only number cells)

tbl_pval <- list()
for (celltype in celltypes){
  print(celltype)
    # take subset of cell type 
    metadata <- subset(metadata_comb,metadata_comb$Celltype == celltype)
    print(nrow(metadata))
    if(nrow(metadata) >= min_cells){
      cells_keep<-intersect(rownames(metadata),colnames(CNA_Matrix_comb))
      CNA_Matrix_comb_celltype <- CNA_Matrix_comb[,cells_keep, drop=F]
      sum(colSums(CNA_Matrix_comb_celltype,na.rm = T))
      
      print(dim(CNA_Matrix_comb_celltype))
      print(dim(metadata))
      
      #### remove NAs
      N_NAs <- apply(CNA_Matrix_comb_celltype,1, function(gene)  length(which(is.na(gene)))) # in how many cells we have an amplification
      NonNA <- apply(CNA_Matrix_comb_celltype,1, function(gene)  length(na.omit(gene))) # in how many cells we have information (non-NA)  
      CNA_Matrix_comb_celltype <- CNA_Matrix_comb_celltype[NonNA!=0,,drop=F]
      print(dim(CNA_Matrix_comb_celltype))
      
      ### subset the gene sets to the genes we could detect
      Oncogenes <- intersect(Oncogenes_df, rownames(CNA_Matrix_comb_celltype))
      TS <- intersect(TumorSuppressorGene_df, rownames(CNA_Matrix_comb_celltype))
      
      #### aggregate across cells
      # Amplification
      Abs_amp <- apply(CNA_Matrix_comb_celltype,1, function(gene)  length(which(gene > 0))) # in how many cells we have an amplification
      Abs_NOTamp <- apply(CNA_Matrix_comb_celltype,1, function(gene)  length(which(gene <= 0))) # in how many cells we don't have an amplification
      Abs_del <- apply(CNA_Matrix_comb_celltype,1, function(gene)  length(which(gene < 0))) # in how many cells we have a deletion
      Abs_NOTdel <- apply(CNA_Matrix_comb_celltype,1, function(gene)  length(which(gene >= 0))) # in how many cells we don't have an deletion
      NonNA <- apply(CNA_Matrix_comb_celltype,1, function(gene)  length(na.omit(gene))) # in how many cells we have information (non-NA)  
      Abs_normal <- apply(CNA_Matrix_comb_celltype,1, function(gene)  length(which(gene == 0))) # in how many cells we have information (non-NA) and normal CNA state 
      # how many cells do we have in our dataset
      n_cells <- ncol(CNA_Matrix_comb_celltype)
      # how many genes do we have in our dataset (remove genes where he have NAs across all cells)
      n_genes <- nrow(CNA_Matrix_comb_celltype)
      # N = n_genes * n_cells
      N = n_genes * n_cells
      
      
      
      ### Fisher test sum over all cells and genes
      pval <-c()
      odds_ratio<- c()
      
      nonOncogenes <-rownames(CNA_Matrix_comb_celltype)[!(rownames(CNA_Matrix_comb_celltype) %in% Oncogenes)]
      nonTS <-rownames(CNA_Matrix_comb_celltype)[!(rownames(CNA_Matrix_comb_celltype) %in% TS)]
      
      if(length(which(Abs_amp !=0)) > 0){
        # Question: Are amplifications enriched in oncogenes? ductal yes
        Amp_oncogenes <- sum(Abs_amp[Oncogenes]) # Number of cells with Amp in onco genes
        Amp_nononcogenes <- sum(Abs_amp[nonOncogenes]) # Number of cells with Amp in non-onco genes
        NotAmp_oncogenes <- sum(Abs_NOTamp[Oncogenes]) # Number of cells without Amp in onco genes
        NotAmp_NONoncogenes <- sum(Abs_NOTamp[nonOncogenes]) # Number of cells without Amp in onco genes
        
        df_Amp_oncogenes <- matrix(c(Amp_oncogenes,Amp_nononcogenes,NotAmp_oncogenes,NotAmp_NONoncogenes),nrow=2)
        colnames(df_Amp_oncogenes) <- c("Amplification","No amplification")
        rownames(df_Amp_oncogenes) <- c("Oncogenes","Non-oncogenes")
        #fisher.test(df_Amp_oncogenes)
        #table2redmine(df_Amp_oncogenes)
        pval<-append(pval,fisher.test(df_Amp_oncogenes)$p.value)
        odds_ratio<-append(odds_ratio,fisher.test(df_Amp_oncogenes)$estimate)
        
        # Question: Are amplifications enriched in TS? 
        Amp_TS <- sum(Abs_amp[TS]) # Number of cells with Amp in onco genes
        Amp_nonTS <- sum(Abs_amp[nonTS]) # Number of cells with Amp in non-onco genes
        NotAmp_TS <- sum(Abs_NOTamp[TS]) # Number of cells without Amp in onco genes
        NotAmp_NONTS <- sum(Abs_NOTamp[nonTS]) # Number of cells without Amp in onco genes
        
        df_Amp_TS <- matrix(c(Amp_TS,Amp_nonTS,NotAmp_TS,NotAmp_NONTS),nrow=2)
        colnames(df_Amp_TS) <- c("Amplification","No amplification")
        rownames(df_Amp_TS) <- c("TS","Non-TS")
        #fisher.test(df_Amp_TS)
        pval<-append(pval,fisher.test(df_Amp_TS)$p.value)
        odds_ratio<-append(odds_ratio,fisher.test(df_Amp_TS)$estimate)
      }else{
        pval<-append(pval,c(NA,NA))
        odds_ratio<-append(odds_ratio,c(NA,NA))
      }
      if(length(which(Abs_del !=0)) > 0){
        # Question: Are Deletions enriched in Oncogenes? 
        del_oncogenes <- sum(Abs_del[Oncogenes]) # Number of cells with del in onco genes
        del_nononcogenes <- sum(Abs_del[nonOncogenes]) # Number of cells with del in non-onco genes
        Notdel_oncogenes <- sum(Abs_NOTdel[Oncogenes]) # Number of cells without del in onco genes
        Notdel_NONoncogenes <- sum(Abs_NOTdel[nonOncogenes]) # Number of cells without del in onco genes
        
        df_del_oncogenes <- matrix(c(del_oncogenes,del_nononcogenes,Notdel_oncogenes,Notdel_NONoncogenes),nrow=2)
        colnames(df_del_oncogenes) <- c("Deletion","No Deletion")
        rownames(df_del_oncogenes) <- c("Oncogenes","Non-oncogenes")
        #fisher.test(df_del_oncogenes)
        pval<-append(pval,fisher.test(df_del_oncogenes)$p.value)
        odds_ratio<-append(odds_ratio,fisher.test(df_del_oncogenes)$estimate)
        
        # Question: Are Deletions enriched in TS? ductal cells yes
        del_TS <- sum(Abs_del[TS]) # Number of cells with del in onco genes
        del_nonTS <- sum(Abs_del[nonTS]) # Number of cells with del in non-onco genes
        Notdel_TS <- sum(Abs_NOTdel[TS]) # Number of cells without del in onco genes
        Notdel_NONTS <- sum(Abs_NOTdel[nonTS]) # Number of cells without del in onco genes
        
        df_del_TS <- matrix(c(del_TS,del_nonTS,Notdel_TS,Notdel_NONTS),nrow=2)
        colnames(df_del_TS) <- c("Deletion","No Deletion")
        rownames(df_del_TS) <- c("TS","Non-TS")
        #fisher.test(df_del_TS)
        pval<-append(pval,fisher.test(df_del_TS)$p.value)
        odds_ratio<-append(odds_ratio,fisher.test(df_del_TS)$estimate)
      }else{
        pval<-append(pval,c(NA,NA))
        odds_ratio<-append(odds_ratio,c(NA,NA))
      }
      names(pval) <- c("Amplifications + oncogenes","Amplifications + TS","Deletions + Oncogenes","Deletions + TS")
      pvaladjust <- p.adjust(pval,method = 'fdr')
      tbl_pval[[celltype]]<- data.frame("Test"=names(pval),"pvalue" = pval,"FDR" = pvaladjust,'OddsRatio'=odds_ratio,
                                                  'CNA'=c(rep('Amplification',2),rep('Deletion',2)),'Gene'=c(rep(c('Onkogenes','Tumor suppressor genes'),2)),
                                                  'Celltype'=celltype)
      
    
    
  }
}


tbl_pval_comb <- do.call(rbind,tbl_pval)
saveRDS(tbl_pval_comb,'output/FDR01/Baron_Muraro_TS_Peng/tbl_pval_comb.rds')
text<-as.factor(tbl_pval_comb$FDR < 0.05)
text<-revalue(text,c('TRUE'='*','FALSE'=''))

library(viridis)
png(paste0('plots/FDR01/Baron_Muraro_TS_Peng/OddsRatio_per_donor_per_celltype2.png'),width = 600)
ggplot(tbl_pval_comb)+geom_bar(aes(y=OddsRatio,x=Gene,fill=FDR),stat='identity')+geom_hline(yintercept=1,linetype='dashed')+facet_grid(CNA~Celltype)+
  xlab('') + geom_text(aes(y=OddsRatio,x=Gene, label = text),size=10)+theme_readable()+
  scale_fill_viridis(direction = -1, limits=c(0,0.1), na.value="#440154FF")+ theme(axis.text.x = element_text(angle = 45, hjust=1),
                                                                                   strip.text.y = element_text(size = 15))
dev.off()
pdf(paste0('plots/FDR01/Baron_Muraro_TS_Peng/OddsRatio_per_donor_per_celltype2.pdf'),width = 10)
ggplot(tbl_pval_comb)+geom_bar(aes(y=OddsRatio,x=Gene,fill=FDR),stat='identity')+geom_hline(yintercept=1,linetype='dashed')+facet_grid(CNA~Celltype)+
  xlab('') + geom_text(aes(y=OddsRatio,x=Gene, label = text),size=10)+theme_readable()+
  scale_fill_viridis(direction = -1, limits=c(0,0.1), na.value="#440154FF")+ theme(axis.text.x = element_text(angle = 45, hjust=1),
                                                                                   strip.text.y = element_text(size = 15))
dev.off()
