#' @title Create Heatmap
#' @usage create_Heatmaps(CNA_profile, metadata, gencode, path, file_name)
#'
#' @param CNA_profile CNA_profile
#' @param metadata data.frame; col 1 = Donors, col 2 = celltypes, rownames = cell_names
#' @param gencode gene annotation
#' @param path path
#' @param file_name file_name
#'
#' @return heatmap
#'
#' @description Create Heatmap
#'
#' @export
#'
create_Heatmaps <- function(CNA_profile, metadata, gencode, path, file_name, limit_values=TRUE){

  ###### Annotation bars
  # order cells

  if(ncol(metadata) > 2){
    order_index <- order(metadata[[3]],metadata[[1]],metadata[[2]])
  }else{
    order_index <- order(metadata[[1]],metadata[[2]])
  }

  metadata <- metadata[order_index,]
  CNA_profile <- CNA_profile[,order_index]
  Cell_names_in_order <- rownames(metadata)
  # change color
  bar1 <- viridis(length(unique(metadata[[1]])))
  names(bar1) <- unique(metadata[[1]])
  bar2 <- magma(length(unique(metadata[[2]])))
  names(bar2) <- unique(metadata[[2]])

  bar1_name <- paste(colnames(metadata)[1])
  bar2_name <- paste(colnames(metadata)[2])

    if(ncol(metadata) > 2){
    bar3 <- c("black","#c11509")
    names(bar3) <- unique(metadata[[3]])
    bar3_name <- paste(colnames(metadata)[3])

    # create Annotation bar
    row_ha = rowAnnotation("bar3" = metadata[[3]],"bar1" = metadata[[1]], "bar2" = metadata[[2]],
                           col = list("bar1" = bar1, "bar2" = bar2, "bar3" = bar3),
                           annotation_name_gp= gpar(fontsize = 18),
                           show_annotation_name=FALSE,
                           simple_anno_size = unit(5, "mm"),
                           annotation_legend_param = list("bar3"  = list(title = bar3_name, title_gp = gpar(fontsize = 15,fontface = "bold"),
                                                                         labels_gp = gpar(fontsize = 15),annotation_legend_side = "top"),
                                                          "bar1"  = list(title = bar1_name, title_gp = gpar(fontsize = 15,fontface = "bold"),
                                                                         labels_gp = gpar(fontsize = 15),annotation_legend_side = "top"),
                                                          "bar2"  = list(title = bar2_name, title_gp = gpar(fontsize = 15,fontface = "bold"),
                                                                         labels_gp = gpar(fontsize = 15),annotation_legend_side = "top")))
  }else{
  # create Annotation bar
  row_ha = rowAnnotation("bar1" = metadata[[1]], "bar2" = metadata[[2]],
                         col = list("bar1" = bar1, "bar2" = bar2),
                         annotation_name_gp= gpar(fontsize = 18),
                         show_annotation_name=FALSE,
                         simple_anno_size = unit(5, "mm"),
                         annotation_legend_param = list("bar1"  = list(title = bar1_name, title_gp = gpar(fontsize = 15,fontface = "bold"),
                                                                       labels_gp = gpar(fontsize = 15),annotation_legend_side = "top"),
                                                        "bar2"  = list(title = bar2_name, title_gp = gpar(fontsize = 15,fontface = "bold"),
                                                                        labels_gp = gpar(fontsize = 15),annotation_legend_side = "top")))
  }
  if(limit_values == TRUE){
    CNA_profile[CNA_profile > 5] <- 5
    CNA_profile[CNA_profile < -5] <- -5
  }
  # make sure that the CNAs are sorted by chromosome and gene location
  if(str_detect(gencode$gene_id[1], "ENSG")){n_autosome=22}else if(str_detect(gencode$gene_id[1], "ENSMUSG")){n_autosome=19}
  CNA_profile <- sort_chromosomes(CNA_profile, gencode)
  gencode <- CNA_profile[[2]]
  CNA_profile <- CNA_profile[[1]]
  # border chromosomes
  gencode$seqid <- sapply(as.character(gencode$seqid), function(x) strsplit(x,"chr")[[1]][2])
  gencode$seqid <- as.factor(gencode$seqid)
  gencode$seqid <- factor(gencode$seqid,levels=c(1:n_autosome,"X"))

  # turn off warning message
  ht_opt$message <- FALSE
  # change the position of the column and row titles
  ht_opt$TITLE_PADDING = unit(c(20, 10), "points")

  if(min(CNA_profile,na.rm = T) == 0){
    col_fun = colorRamp2(c(min(CNA_profile,na.rm = T), max(CNA_profile,na.rm = T)), c("white", "red"))
  }else if(max(CNA_profile,na.rm = T) == 0){
    col_fun = colorRamp2(c(min(CNA_profile,na.rm = T), max(CNA_profile,na.rm = T)), c("blue", "white"))
  }else{
    col_fun = colorRamp2(c(min(CNA_profile,na.rm = T), 0, max(CNA_profile,na.rm = T)), c("blue", "white", "red"))
  }
  png(paste0(path, file_name, ".png"),width = 800)
  h<-Heatmap(t(CNA_profile), cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE,
             border = TRUE, column_split=gencode$seqid , row_title = "Cells",
             use_raster=TRUE,row_order=Cell_names_in_order,
             column_order = gencode$gene_name, left_annotation = row_ha, column_title_side = "bottom", row_title_side = "right",
             column_title_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 20),col=col_fun,
             heatmap_legend_param = list(labels_gp = gpar(fontsize = 15), title_gp = gpar(fontsize = 15,fontface = "bold"),
                                         title = "CNA score", grid_width = unit(1, "cm"), legend_height = unit(4, "cm"))) #%v% NULL
  draw(h,annotation_legend_side = "right")
  dev.off()
}

#' Function to sort the CNA profile according to the genomic position of the genes
#'
#' @param CNA_profile CNA_profile
#' @param gencode gene annotation
#'
#' @return sorted CNA profile
#' @export
#'
sort_chromosomes <- function(CNA_profile, gencode){
  if(str_detect(gencode$gene_id[1], "ENSG")){n_autosome=22}else if(str_detect(gencode$gene_id[1], "ENSMUSG")){n_autosome=19}
  CNA_profile <- as.matrix(CNA_profile)
  intersecting_genes <- intersect(row.names(CNA_profile),gencode$gene_name)
  CNA_profile <- CNA_profile[intersecting_genes,]
  gencode <- gencode[match(row.names(CNA_profile),gencode$gene_name),]
  gencode$chr <- factor(gsub("chr","",gencode$seqid),levels = c(seq(1,n_autosome),"X"))
  filter <- factor(c(seq(1,n_autosome),"X"),levels = c(seq(1,n_autosome),"X"))
  gencode <- gencode[which(gencode$chr %in% filter),]
  gencode <- gencode[order(gencode$chr, gencode$start),]
  CNA_profile <- CNA_profile[match(gencode$gene_name,row.names(CNA_profile)),]
  print("Sorted in chromosomal order.")
  return(list(CNA_profile,gencode))
}
