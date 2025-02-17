library("dplyr")

Merge_sequential_CNAs <- function(CNA_profile,gencode,Segments = FALSE,give_genes_on_boundaries_same_value=give_genes_on_boundaries_same_value){

  # Set NAs to zero and create a TRUE/FAlSE Matrix for Amp and Del separately
  CNA_profile_Logical <- CNA_profile
  #CNA_profile_Logical_NA <- is.na(CNA_profile_Logical)
  CNA_profile_Logical_zeros <- CNA_profile_Logical == 0
  CNA_profile_Logical_zeros[is.na(CNA_profile_Logical_zeros)] <- FALSE
  # check for amp and del
  CNA_profile_Logical[is.na(CNA_profile_Logical)] <- 0
  CNA_profile_Logical_Amp <- CNA_profile_Logical > 0
  CNA_profile_Logical_Del <- CNA_profile_Logical < 0
  #message('Compute Changepoints for NAs')

  ## set FALSE on chromosome boundaries
  CNAPro <- genome_sort(CNA_profile, gencode)
  chr_split<-split(chromosomeorder,chromosomeorder$chr)
  for(chr in 1:(length(chr_split)-1)){
    last_gene_chr<-chr_split[[chr]][nrow(chr_split[[chr]]),'hgnc']
    if(chr == 23){
      CNA_profile_Logical_Del[which(rownames(CNA_profile_Logical_Del) == last_gene_chr),] <- FALSE # set FALSE on chromosome boundaries
      CNA_profile_Logical_Amp[which(rownames(CNA_profile_Logical_Amp) == last_gene_chr),] <- FALSE
      CNA_profile_Logical_zeros[which(rownames(CNA_profile_Logical_zeros) == last_gene_chr),] <- FALSE
    }else{
      CNA_profile_Logical_Del[which(rownames(CNA_profile_Logical_Del) == last_gene_chr)+1,] <- FALSE # set FALSE on chromosome boundaries
      CNA_profile_Logical_Amp[which(rownames(CNA_profile_Logical_Amp) == last_gene_chr)+1,] <- FALSE
      CNA_profile_Logical_zeros[which(rownames(CNA_profile_Logical_zeros) == last_gene_chr)+1,] <- FALSE
    }
  }

  #changepoints_NA <- apply(CNA_profile_Logical_NA,2, function(col) Compute_Changepoints(col))
  message('Compute Changepoints for zeros')
  changepoints_zeros <- apply(CNA_profile_Logical_zeros,2, function(col) Compute_Changepoints(col))
  message('Compute Changepoints for Amplifications')
  changepoints_amp <- apply(CNA_profile_Logical_Amp,2, function(col) Compute_Changepoints(col))
  message('Compute Changepoints for Deletions')
  changepoints_del <- apply(CNA_profile_Logical_Del,2, function(col) Compute_Changepoints(col))

  message('Compute Changepoints done!')

  message('Merge Sequential CNAs')
  CNA_profile_joined <- CNA_profile

  #CNA_profile_joined <- Merge_sequential_CNAs_per_Cell(changepoints_NA,CNA_profile_joined)
  if(length(changepoints_zeros) != 0){  CNA_profile_joined <- Merge_sequential_CNAs_per_Cell(changepoints_zeros,CNA_profile_joined)}
  if(length(changepoints_amp) != 0){  CNA_profile_joined <- Merge_sequential_CNAs_per_Cell(changepoints_amp,CNA_profile_joined)}
  if(length(changepoints_del) != 0){  CNA_profile_joined <- Merge_sequential_CNAs_per_Cell(changepoints_del,CNA_profile_joined)}
  message('Done!')

  if(Segments == TRUE){
    message('Compute Segments')
    if(length(changepoints_amp) != 0){  CNA_profile_Segments_Amp <- Create_Segment(changepoints_amp,CNA_profile_joined,gencode,'Amp')}else{CNA_profile_Segments_Amp <- NULL}
    if(length(changepoints_del) != 0){  CNA_profile_Segments_Del <- Create_Segment(changepoints_del,CNA_profile_joined,gencode,'Del')}else{CNA_profile_Segments_Del <- NULL}
    #CNA_profile_Segments_NA <- Create_Segment(changepoints_NA,CNA_profile_joined,gencode,'NA')
    if(length(changepoints_zeros) != 0){  CNA_profile_Segments_zeros <- Create_Segment(changepoints_zeros,CNA_profile_joined,gencode,'Zero')}else{CNA_profile_Segments_zeros <- NULL}

    CNA_profile_Segments <- rbind(CNA_profile_Segments_Amp,CNA_profile_Segments_Del, CNA_profile_Segments_zeros)
    message('Done!')
    return(list(CNA_profile_joined,CNA_profile_Segments))

  }else{
    return(CNA_profile_joined)}
}


Compute_Changepoints <- function(CNA_profile_Logical) {
  changepoints = integer()
  if(CNA_profile_Logical[1]){
    changepoints = c(changepoints, 1)
  }
  for(gene in 2:length(CNA_profile_Logical)){
    if(CNA_profile_Logical[gene-1] != CNA_profile_Logical[gene]){
      #browser()

      changepoints = c(changepoints, gene)
    }
  }
  return(changepoints)
}

Merge_sequential_CNAs_per_Cell <- function(changepoints, CNA_profile_joined) {

  for(cell in 1:ncol(CNA_profile_joined)){
    #print(cell)
    changepoints_cell <- changepoints[cell][[1]]
    if(length(changepoints_cell) > 1 ){
      for (change_index in 1:(length(changepoints_cell)) ){
        if (!(change_index %% 2) == 0){
          CNA_index <- changepoints_cell[change_index]:(changepoints_cell[change_index+1]-1)
          CNA_value <- mean(CNA_profile_joined[CNA_index,cell],na.rm=T)
          CNA_profile_joined[CNA_index,cell] <- CNA_value
        }
      }
    }
  }
  return(CNA_profile_joined)
}


Create_Segment <- function(changepoints, CNA_profile_joined,gencode,CNA) {
  CNAPro <- genome_sort(CNA_profile_joined, gencode)
  Segment_List <- list()
  #browser()
  for(cell in 1:ncol(CNA_profile_joined)){

    changepoints_cell <- changepoints[cell][[1]]
    if(length(changepoints_cell) != 0 ){
      CNA_start <- c()
      CNA_end <- c()
      CNA_value <- c()
      Gene_start<-c()
      Gene_end<-c()
      Chr <- c()
      start <- c()
      end <- c()
      Cell <- c()
      for (change_index in 1:length(changepoints_cell)){
        if (!(change_index %% 2) == 0){
          CNA_start <- append(CNA_start, changepoints_cell[change_index])
          CNA_end <- append(CNA_end,(changepoints_cell[change_index+1]-1))
          CNA_index <- changepoints_cell[change_index]:(changepoints_cell[change_index+1]-1)
          CNA_value <- append(CNA_value, mean(CNA_profile_joined[CNA_index,cell]))
          Gene_start <- append(Gene_start, rownames(CNA_profile_joined)[changepoints_cell[change_index]])
          Gene_end <- append(Gene_end, rownames(CNA_profile_joined)[(changepoints_cell[change_index+1]-1)])
          Chr <- append(Chr, subset(chromosomeorder,hgnc==rownames(CNA_profile_joined)[changepoints_cell[change_index]])$chr)
          start <- append(start, subset(chromosomeorder,hgnc==rownames(CNA_profile_joined)[changepoints_cell[change_index]])$start)
          end <- append(end, subset(chromosomeorder,hgnc==rownames(CNA_profile_joined)[(changepoints_cell[change_index+1]-1)])$end)
          Cell <- append(Cell, colnames(CNA_profile_joined)[cell])

        }
      }

      Segment_List[[cell]] <- data.frame('Cell'=Cell, 'Chr'=Chr, 'Gene_start'=Gene_start,
                                         'Gene_end'=Gene_end,
                                         'Index_start'=CNA_start,'Index_end'=CNA_end,
                                         'loc.start'=start,
                                         'loc.end'=end,
                                         'CNA_value'=CNA_value,'CNA'=CNA)
    }
  }
  #names(Segment_List) <- colnames(CNA_profile_joined)
  Segment_List <- Segment_List[!sapply(Segment_List,is.null)]
  Segment_df <- bind_rows(Segment_List)
  return(Segment_df)
  #return(Segment_List)
}

# Function to sort data frames according to the genomic position of the genes
## Input:       data = data matrix; rows -> genes_hgnc
##              gencode = gencode object
## Output:      matrix of genes sorted by chromosomal order
genome_sort <- function(data, gencode){
  data <- as.matrix(data)
  isect <- intersect(row.names(data),gencode$gene_name)
  data <- data[isect,]
  gencode <- gencode[match(row.names(data),gencode$gene_name),]
  gencode$chr <- factor(gsub("chr","",gencode$seqid),levels=c(seq(1,22),"X","Y"))
  filter <- factor(c(seq(1,22),"X","Y"),levels=c(seq(1,22),"X","Y"))
  gencode <- gencode[which(gencode$chr %in% filter),]
  gencode <- gencode[order(gencode$chr,gencode$start),]
  gencode$ensembl <- sapply(gencode$gene_id,function(x) strsplit(x,"\\.")[[1]][1])
  data <- data[match(gencode$gene_name,row.names(data)),]
  chromosomeorder <- data.frame(chr=gencode$chr,
                                start=gencode$start,
                                end=gencode$end,
                                ensembl=gencode$ensembl,
                                hgnc=gencode$gene_name)
  assign("chromosomeorder", chromosomeorder, .GlobalEnv)
  return(data)
}
