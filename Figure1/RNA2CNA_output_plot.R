source.dir <- "../../../Analysis/"
source(paste0(source.dir,"packages.R"))
source(paste0(source.dir,"functions.R"))
library(BBmisc)

#' Function to sort the CNA profile according to the genomic position of the genes
#'
#' @param CNA_profile CNA_profile
#' @param gencode gene annotation
#'
#' @return sorted CNA profile
#' @export
#'
sort_chromosomes <- function(CNA_profile, gencode){
  CNA_profile <- as.matrix(CNA_profile)
  intersecting_genes <- intersect(row.names(CNA_profile),gencode$gene_name)
  CNA_profile <- CNA_profile[intersecting_genes,]
  gencode <- gencode[match(row.names(CNA_profile),gencode$gene_name),]
  gencode$chr <- factor(gsub("chr","",gencode$seqid),levels = c(seq(1,22),"X","Y"))
  filter <- factor(c(seq(1,22),"X"),levels = c(seq(1,22),"X"))
  gencode <- gencode[which(gencode$chr %in% filter),]
  gencode <- gencode[order(gencode$chr, gencode$start),]
  CNA_profile <- CNA_profile[match(gencode$gene_name,row.names(CNA_profile)),]
  chr_boundaries <- c(0)
  
  for (chr in 1:23){
    boundary <- sum(table(as.numeric(gencode$chr))[1:chr])
    chr_boundaries[chr] <- boundary
  }
  chr_boundaries <- na.omit(chr_boundaries)
  print("Sorted in chromosomal order.")
  gencode<-droplevels(gencode)
  chromosomeorder <- data.frame(chr=gencode$chr,
                                start=gencode$start,
                                end=gencode$end,
                                ensembl=gencode$gene_id,
                                hgnc=gencode$gene_name)
  assign("chromosomeorder", chromosomeorder, .GlobalEnv)
  return(list(CNA_profile,  chr_boundaries))
}

#######
###  Visualisierung

## Gencode
gencode <- readRDS(file='/cellfile/cellnet/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_ZINBWaVE_NetCNA_CBS/Data/gencode_v28.rds')

# output
smoothed <- readRDS(file='/cellfile/cellnet/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_CBS/Tabula_Sapiens/Marrow/output/scaledata/v2/genewise/FixedMovingWindow/CNAs_per_gene.rds')
data <- readRDS(file='/cellfile/cellnet/Ronja/Tumor_Evolution/Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_CBS/Tabula_Sapiens/Marrow/output/scaledata/v2/genewise/FixedMovingWindow/data.rds')
metadata <- data@meta.data

CBS <- readRDS(file="../../../Analysis/Aim1_Baseline_Mutation_Rate/CNA_Calling_Workflow_scTransform_RNA2CNA_fixedWindow/Analysis/Onkogenic_CNAs/output/TS/FDR01/Marrow_CBS_pvalues_FDR01_Matrix_Imputed.rds")

norm_expr <- data@assays$SCT$scale.data

numbers <- c(1:nrow(metadata))
new_names <- sapply(1:nrow(metadata),  function(row) paste( metadata[row,]$donor,  metadata[row,]$organ_tissue, metadata[row,]$cell_ontology_class)) # ,'Nr.',numbers[row]
names(new_names) <- rownames(metadata)
chromosome_names <- data.frame('chr'= c(1:22,'X'), 'posx'= c(0.34,0.37,0.4,0.43,0.46,0.48,0.5,0.52,0.535,0.55, 0.57,0.6,0.63,0.64,0.66,0.68,0.72,0.76,0.8,0.84,0.87,0.9,0.93),
                               'posy'=c(rep(0.04,12),0.02,rep(0.04,7),0.02,rep(0.04,2)))

nrow(norm_expr)
nrow(CBS)
nrow(smoothed)

# do not use imputed values
CBS <- CBS[rownames(smoothed),]
norm_expr <- norm_expr[rownames(smoothed),]

chr_b<-sort_chromosomes(norm_expr,gencode)
all.equal(rownames(norm_expr),rownames(chr_b[[1]]))
CBS <- CBS[rownames(chr_b[[1]]),]
norm_expr <- norm_expr[rownames(chr_b[[1]]),]
smoothed <- smoothed[rownames(chr_b[[1]]),]

color_chr <- c()
for(chr in 1:length(chr_b[[2]])){
  if(chr == 1){  color_chr <- append(color_chr, rep('black',chr_b[[2]][chr]) )}
  else if ((chr %% 2) == 0){
    color_chr <- append(color_chr, rep('darkgrey',chr_b[[2]][chr]-chr_b[[2]][chr-1]) )
  }else(
    color_chr <- append(color_chr, rep('black',chr_b[[2]][chr]-chr_b[[2]][chr-1]) )
  )
}
NA_rep <-nrow(CBS)-length(color_chr)
color_chr <- c(color_chr,rep('NA',NA_rep))

CBS[CBS == 0] <- NA

sort(colSums(CBS,na.rm = T))[1:10]

cells_with_amp_and_del <- apply(CBS,2, function(cell) all(length(which(cell < 0)) !=0, length(which(cell > 0)) !=0) )
which(cells_with_amp_and_del == TRUE)
sort(colSums(CBS,na.rm = T)[which(cells_with_amp_and_del == TRUE)])

CBS_Amp <- CBS
CBS_Amp[CBS_Amp < 0] <- NA
CBS_Del <- CBS
CBS_Del[CBS_Del > 0] <- NA

# CAGGTATAGAAGGGAT_TSP14_BoneMarrow_NA_10X_1_1 TCCTGCATCTCCGCAT_TSP14_BoneMarrow_NA_10X_1_1
example_cell <- "ATCGGCGCAACTGTGT_TSP2_BM_vertebralbody_10X_2_1"

# select single cells
test_cell <- data.frame('CNA_Amp'=CBS_Amp[,"CAGGTATAGAAGGGAT_TSP14_BoneMarrow_NA_10X_1_1"],'CNA_Del'=CBS_Del[,"CAGGTATAGAAGGGAT_TSP14_BoneMarrow_NA_10X_1_1"],
                                     'Gene'=c(1:length(CBS[,"CAGGTATAGAAGGGAT_TSP14_BoneMarrow_NA_10X_1_1"])),
                                     'Raw'=norm_expr[,"CAGGTATAGAAGGGAT_TSP14_BoneMarrow_NA_10X_1_1"], 'Smoothed'=smoothed[,"CAGGTATAGAAGGGAT_TSP14_BoneMarrow_NA_10X_1_1"])
title <- new_names["CAGGTATAGAAGGGAT_TSP14_BoneMarrow_NA_10X_1_1"]

test_cell <- data.frame('CNA_Amp'=CBS_Amp[,"TCCTGCATCTCCGCAT_TSP14_BoneMarrow_NA_10X_1_1"],'CNA_Del'=CBS_Del[,"TCCTGCATCTCCGCAT_TSP14_BoneMarrow_NA_10X_1_1"],
                        'Gene'=c(1:length(CBS[,"TCCTGCATCTCCGCAT_TSP14_BoneMarrow_NA_10X_1_1"])),
                        'Raw'=norm_expr[,"TCCTGCATCTCCGCAT_TSP14_BoneMarrow_NA_10X_1_1"], 'Smoothed'=smoothed[,"TCCTGCATCTCCGCAT_TSP14_BoneMarrow_NA_10X_1_1"])
title <- new_names["TCCTGCATCTCCGCAT_TSP14_BoneMarrow_NA_10X_1_1"]

# select single cells
test_cell <- data.frame('CNA_Amp'=CBS_Amp[,"TGTACAGTCACCCATC_TSP14_BoneMarrow_NA_10X_1_1"],'CNA_Del'=CBS_Del[,"TGTACAGTCACCCATC_TSP14_BoneMarrow_NA_10X_1_1"],
                        'Gene'=c(1:length(CBS[,"TGTACAGTCACCCATC_TSP14_BoneMarrow_NA_10X_1_1"])),
                        'Raw'=norm_expr[,"TGTACAGTCACCCATC_TSP14_BoneMarrow_NA_10X_1_1"], 'Smoothed'=smoothed[,"TGTACAGTCACCCATC_TSP14_BoneMarrow_NA_10X_1_1"])
title <- new_names["TGTACAGTCACCCATC_TSP14_BoneMarrow_NA_10X_1_1"]

test_cell <- data.frame('CNA_Amp'=CBS_Amp[,example_cell],'CNA_Del'=CBS_Del[,example_cell],
                        'Gene'=c(1:length(CBS[,example_cell])),
                        'Raw'=norm_expr[,example_cell], 'Smoothed'=smoothed[,example_cell])
title <- new_names[example_cell]
title <- paste(strsplit(title,"_")[[1]],collapse = " ")

pdf(paste0('plots/','example_cell_smoothed_',title,'.pdf'), width = 4.5,height = 2)
p <- ggplot(test_cell) + geom_point(aes(x=Gene,y=Smoothed,color='Smoothed'),size=0, col=color_chr) + 
  xlab("\n Chromosome")+ ylab('CNA Score') + theme_manuscript() + 
  theme(plot.title = element_text(hjust = 0.5,size=10),legend.text = element_text(size = 10),legend.title = element_blank(),
        legend.position = 'bottom',axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),panel.background = element_rect(fill = "white")) +
  geom_point(aes(x=Gene,y=CNA_Amp,color='CNA_Amp')) + 
  geom_point(aes(x=Gene,y=CNA_Del,color='CNA_Del')) + 
  geom_hline(aes(color='Smoothed',yintercept = 0),alpha=0) +
  geom_hline(yintercept = 0,size=0.5,col='grey',linetype='dashed') +
  scale_color_manual(name = "CNA score",
                     values = c( CNA_Amp = "#D16103",CNA_Del='#005AB5', Smoothed ='gray35' ),
                     labels = c("Amplification", "Deletion","Smoothed"))+
  ggtitle(title)+  
  geom_vline(xintercept = chr_b[[2]][1:23],col="darkgrey",size=0.5) + ylim(min(test_cell$Smoothed),max(test_cell$Smoothed))
print(p)

chromosome_names$posx <- normalize(chr_b[[2]], range=c(0.2,0.94),method = 'range')#range=c(0.335,0.93)
chromosome_names$posy <- c(rep(c(0.38,0.33),11),0.35)
for(chr in 1:length(chromosome_names$chr)){
  grid.text(chromosome_names$chr[chr], x = unit(chromosome_names$posx[chr], "npc"), y = unit(chromosome_names$posy[chr], "npc"),gp=gpar(fontsize=10, col="black"))
  #text(x= chromosome_names$posx[chr], y=chromosome_names$posy[chr],
  #     label=chromosome_names$chr[chr])
}
dev.off()

pdf(paste0('plots/','example_cell_raw_',title,'.pdf'), width = 12,height = 12)
p <- ggplot(test_cell) + geom_point(aes(x=Gene,y=Raw,color='Raw'),size=0.5, col=color_chr) + 
  xlab("\n Chromosome")+ ylab('CNA Score') + theme_readable() + 
  theme(plot.title = element_text(hjust = 0.5,size=20),legend.text = element_text(size = 20),legend.title = element_blank(),
        legend.position = 'bottom',axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_point(aes(x=Gene,y=Smoothed,color='Smoothed'), size=0.5) + 
  geom_point(aes(x=Gene,y=CNA_Amp,color='CNA_Amp')) + 
  geom_point(aes(x=Gene,y=CNA_Del,color='CNA_Del')) + 
  geom_hline(aes(color='Raw',yintercept = 0),alpha=0) +
  geom_hline(yintercept = 0,size=0.5,col='gray35',linetype='dashed') +
  scale_color_manual(name = "CNA score",
                     values = c( CNA_Amp = "#D16103",CNA_Del='#005AB5', Raw ='gray35', Smoothed="red2" ),
                     labels = c("Amplification", "Deletion","Normalized expression","Smoothed"))+
  ggtitle(title)+  
  geom_vline(xintercept = chr_b[[2]][1:23],col="darkgrey") + ylim(min(test_cell$Raw),max(test_cell$Raw))
print(p)

chromosome_names$posx <- (normalize(chr_b[[2]], range=chr_b[[2]][c(1,23)],method = 'range')/10000)+0.07 #range=c(0.335,0.93)
chromosome_names$posy <- c(rep(c(0.1,0.085),11),0.1)
for(chr in 1:length(chromosome_names$chr)){
  grid.text(chromosome_names$chr[chr], x = unit(chromosome_names$posx[chr], "npc"), y = unit(chromosome_names$posy[chr], "npc"),gp=gpar(fontsize=15, col="black"))
  #text(x= chromosome_names$posx[chr], y=chromosome_names$posy[chr],
  #     label=chromosome_names$chr[chr])
}
dev.off()


#NetCNA_Peng2019['TP53',which(metadata$Celltype == 'Ductal cell type 2')]
which(CBS["TP53",] > 0)
CBS["TP53",which(CBS["TP53",] < 0)]
