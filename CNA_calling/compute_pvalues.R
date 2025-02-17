#' @title Compute p-values for each CNA segment
#' @usage compute_pvalues(CBS, metadata, output_dir)
#'
#' @param CBS CBS output
#' @param metadata metadata
#' @param output_dir output_dir
#' @return CNA segments with p-values
#'
#' @description Compute p-values for each CNA segment
#'
#' @export
#'
compute_pvalues <- function(CBS, CBS_permut, counts, output_dir){
  ### detected genes per cell
  detected_genes <- DetectedGenes(counts, Cell_Ref=NULL)

  # create 6 groups
  Q1<-which(detected_genes <= quantile(detected_genes,0.1666667))
  Q2<-which(detected_genes <= quantile(detected_genes,0.3333334) & detected_genes > quantile(detected_genes,0.1666667) )
  Q3<-which(detected_genes <= quantile(detected_genes,0.5000001) & detected_genes > quantile(detected_genes,0.3333334) )
  Q4<-which(detected_genes <= quantile(detected_genes,0.6666668) & detected_genes > quantile(detected_genes,0.5000001) )
  Q5<-which(detected_genes <= quantile(detected_genes,0.8333335) & detected_genes > quantile(detected_genes,0.6666668) )
  Q6<-which(detected_genes > quantile(detected_genes,0.8333335))

  ## use the CBS output
  # change names
  cell_names<-colnames(counts)
  samples_names<-paste0("Sample.",1: length(cell_names))
  sample_cell_names<-data.frame(cell_names,samples_names)
  colnames(sample_cell_names)<-c('Cell','ID')

  if(all.equal(sample_cell_names$Cell, names(detected_genes))){sample_cell_names$detected_genes <- detected_genes}else{print("sample names not same")}

  sample_cell_names$group <- NA
  sample_cell_names[which(sample_cell_names$Cell %in% names(Q1)),"group"]<-"Q1"
  sample_cell_names[which(sample_cell_names$Cell %in% names(Q2)),"group"]<-"Q2"
  sample_cell_names[which(sample_cell_names$Cell %in% names(Q3)),"group"]<-"Q3"
  sample_cell_names[which(sample_cell_names$Cell %in% names(Q4)),"group"]<-"Q4"
  sample_cell_names[which(sample_cell_names$Cell %in% names(Q5)),"group"]<-"Q5"
  sample_cell_names[which(sample_cell_names$Cell %in% names(Q6)),"group"]<-"Q6"

  sample_cell_names$group<-as.factor(sample_cell_names$group)

  CNA_scores_list<-list()
  for (group in 1:6){
    #print(group)
    cell_names <- sample_cell_names[which(sample_cell_names$group == levels(sample_cell_names$group)[group]),"ID"]
    CNA_scores<-list()
    for (permut in 1:10){
      CNA_scores[[permut]] <- CBS_permut[[permut]]$output[which(CBS_permut[[permut]]$output$ID %in% cell_names),]
    }
    if(nrow(do.call(rbind,CNA_scores)) != 0){
      CNA_scores_list[[group]]<- do.call(rbind,CNA_scores)
      CNA_scores_list[[group]]<-cbind(CNA_scores_list[[group]],"group"=rep(levels(sample_cell_names$group)[group],nrow(CNA_scores_list[[group]])))
    }
  }


  # again split by num.mark
  # compute quantiles based on size within each group of detected genes
  for (group in 1:length(CNA_scores_list)){
    n_genes<-CNA_scores_list[[group]]$num.mark
    n_genes[n_genes==0] <- NA
    CNA_scores_list[[group]]$n_genes <- NA
    CNA_scores_list[[group]]$Qn_genes <- NA
    CNA_scores_list[[group]]$min_genes <- NA
    CNA_scores_list[[group]]$max_genes <- NA

    Q1_n_genes <- which(n_genes < quantile(n_genes,0.25,na.rm=T))
    Q2_n_genes <- which(n_genes < quantile(n_genes,0.5,na.rm=T) & n_genes >= quantile(n_genes,0.25,na.rm=T))
    Q3_n_genes <- which(n_genes < quantile(n_genes,0.75,na.rm=T) & n_genes >= quantile(n_genes,0.5,na.rm=T))
    Q4_n_genes <- which(n_genes >= quantile(n_genes,0.75,na.rm=T))

    CNA_scores_list[[group]][Q1_n_genes,"n_genes"]<-paste(as.character(summary(CNA_scores_list[[group]][Q1_n_genes,"num.mark"])[c(1,6)]),collapse = "-")
    CNA_scores_list[[group]][Q2_n_genes,"n_genes"]<-paste(as.character(summary(CNA_scores_list[[group]][Q2_n_genes,"num.mark"])[c(1,6)]),collapse = "-")
    CNA_scores_list[[group]][Q3_n_genes,"n_genes"]<-paste(as.character(summary(CNA_scores_list[[group]][Q3_n_genes,"num.mark"])[c(1,6)]),collapse = "-")
    CNA_scores_list[[group]][Q4_n_genes,"n_genes"]<-paste(as.character(summary(CNA_scores_list[[group]][Q4_n_genes,"num.mark"])[c(1,6)]),collapse = "-")

    CNA_scores_list[[group]][Q1_n_genes,"Qn_genes"]<-"Q1_n_genes"
    CNA_scores_list[[group]][Q2_n_genes,"Qn_genes"]<-"Q2_n_genes"
    CNA_scores_list[[group]][Q3_n_genes,"Qn_genes"]<-"Q3_n_genes"
    CNA_scores_list[[group]][Q4_n_genes,"Qn_genes"]<-"Q4_n_genes"

    CNA_scores_list[[group]][Q1_n_genes,"min_genes"]<-summary(CNA_scores_list[[group]][Q1_n_genes,"num.mark"])[1]
    CNA_scores_list[[group]][Q2_n_genes,"min_genes"]<-summary(CNA_scores_list[[group]][Q2_n_genes,"num.mark"])[1]
    CNA_scores_list[[group]][Q3_n_genes,"min_genes"]<-summary(CNA_scores_list[[group]][Q3_n_genes,"num.mark"])[1]
    CNA_scores_list[[group]][Q4_n_genes,"min_genes"]<-summary(CNA_scores_list[[group]][Q4_n_genes,"num.mark"])[1]
    CNA_scores_list[[group]][Q1_n_genes,"max_genes"]<-summary(CNA_scores_list[[group]][Q1_n_genes,"num.mark"])[6]
    CNA_scores_list[[group]][Q2_n_genes,"max_genes"]<-summary(CNA_scores_list[[group]][Q2_n_genes,"num.mark"])[6]
    CNA_scores_list[[group]][Q3_n_genes,"max_genes"]<-summary(CNA_scores_list[[group]][Q3_n_genes,"num.mark"])[6]
    CNA_scores_list[[group]][Q4_n_genes,"max_genes"]<-summary(CNA_scores_list[[group]][Q4_n_genes,"num.mark"])[6]

  }
  names(CNA_scores_list) <- unique(sample_cell_names$group)
  CNA_scores_list<-do.call(rbind,CNA_scores_list)
  CNA_scores_list$num.mark[CNA_scores_list$num.mark==0] <- NA


  #plots
  png(paste0("plots/",output_dir,'Permuted_CNA_score_dis_Segment_size.png'))
  p<-ggplot(CNA_scores_list)+geom_density(aes(x=seg.mean,col=Qn_genes))+xlab("CNA scores")+theme_readable()
  print(p)
  dev.off()

  png(paste0("plots/",output_dir,'Permuted_CNA_score_dis_detected_genes.png'))
  p<-ggplot(CNA_scores_list)+geom_density(aes(x=seg.mean,col=group))+xlab("CNA scores")+theme_readable()
  print(p)
  dev.off()

  png(paste0("plots/",output_dir,'Permuted_CNA_score_dis_detected_genes_Segment_size2.png'),width = 800)
  p<-ggplot(na.omit(CNA_scores_list))+geom_density(aes(x=seg.mean))+facet_grid(Qn_genes~group)+xlab("CNA scores")+theme_readable()
  print(p)
  dev.off()

  png(paste0("plots/",output_dir,'Permuted_CNA_score_hist_Segment_size.png'))
  p<-ggplot(CNA_scores_list)+geom_histogram(aes(x=num.mark))+facet_wrap(.~group)+xlab("number genes per segment")+theme_readable()
  print(p)
  dev.off()

  # TRUE CNA
  CBS$output$SegmentID<-1:nrow(CBS$output)
  # change names
  CBS_pvalues<-merge(CBS$output, sample_cell_names,by='ID')
  CBS_pvalues<-CBS_pvalues[order(CBS_pvalues$SegmentID),] # merge changes the order, so I change it back to the inputs order because subsequent function CreateMatrix_fromLargeDNAcopy needs the order from the CBS output

  pval_amp <- c()
  pval_del <- c()
  for (segment in 1:nrow(CBS_pvalues)){
    #print(segment)
    seg <- CBS_pvalues[segment,]
    index<- which(seg$ID == CNA_scores_list$ID)
    # select the groups
    Group <- unique(CNA_scores_list[index,"group"])
    N_genes <- unique(CNA_scores_list[index,][which(CNA_scores_list[index,"min_genes"] <= seg$num.mark & CNA_scores_list[index,"max_genes"] >= seg$num.mark),"n_genes"])
    # choose permut distribution
    Permut_dis <- subset(CNA_scores_list,CNA_scores_list$group == Group & CNA_scores_list$n_genes == N_genes)
    pval_amp <- c(pval_amp, length(which(Permut_dis$seg.mean > seg$seg.mean))/length(na.omit(Permut_dis$seg.mean)) )
    pval_del <-  c(pval_del, length(which(Permut_dis$seg.mean < seg$seg.mean))/length(na.omit(Permut_dis$seg.mean)) )
  }
  CBS_pvalues$pval_amp <- pval_amp
  CBS_pvalues$pval_del <- pval_del

  CBS_pvalues[is.na(CBS_pvalues$seg.mean),"num.mark"]<-NA
  CBS_pvalues[is.na(CBS_pvalues$seg.mean),"pval_amp"]<-NA
  CBS_pvalues[is.na(CBS_pvalues$seg.mean),"pval_del"]<-NA

  # select for every segment hte smalles pvalue
  CBS_pvalues$pval <- 2*as.numeric(sapply(1:nrow(CBS_pvalues), function(x) min(CBS_pvalues[x,c("pval_amp",  "pval_del")])))
  FDR <- p.adjust(CBS_pvalues[!is.na(CBS_pvalues$pval),"pval"],method = "fdr")
  CBS_pvalues$FDR <- NA
  CBS_pvalues[!is.na(CBS_pvalues$pval),"FDR"] <- FDR

  cell_names <- unique(CBS_pvalues$Cell)
  CBS_pvalues$FDRpercell <- NA
  #fdr correction per cell
  for(cell in cell_names){
    pvalues <- CBS_pvalues[which(CBS_pvalues$Cell %in% cell),"pval"]
    pvalue_index <- which(!is.na(pvalues))
    CBS_pvalues[which(CBS_pvalues$Cell %in% cell)[pvalue_index],"FDRpercell"] <- p.adjust(pvalues[pvalue_index],method = "fdr")
  }

  png(paste0("plots/",output_dir,'Hist_pvalues.png'))
  p <- ggplot(CBS_pvalues) + geom_histogram(aes(x=pval),position = "identity",color = "black") + xlab("p-values") + theme_readable()
  print(p)
  dev.off()
  png(paste0("plots/",output_dir,'Hist_pvalues_adjusted.png'))
  p <-ggplot(CBS_pvalues) + geom_histogram(aes(x=FDR),position = "identity",color = "black") + xlab("FDR") + theme_readable()
  print(p)
  dev.off()
  png(paste0("plots/",output_dir,'Hist_pvalues_adjusted_per_cell.png'))
  p <-ggplot(CBS_pvalues) + geom_histogram(aes(x=FDRpercell),position = "identity",color = "black") + xlab("FDR per cell") + theme_readable()
  print(p)
  dev.off()

  CBS_pvalues <- CBS_pvalues[order(CBS_pvalues$SegmentID),]
  return(CBS_pvalues)
}
