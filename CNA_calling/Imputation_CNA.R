Imputation_CNA <- function(CNA_Matrix, CNA_Segments,gencode_gene_types){

  CNA_Segments$Chr<-  as.factor(CNA_Segments$Chr)
  CNA_Segments$Chr<-  factor(CNA_Segments$Chr,levels = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','19','20','22','X'))

  levels(CNA_Segments$Chr)
  #Split for each Cell
  CNA_Segments_split <- split(CNA_Segments, CNA_Segments$Cell)

  # create empty dataframe
  Imputed_CNA_Matrix <- data.frame(matrix(NA, ncol = ncol(CNA_Matrix), nrow = nrow(gencode_gene_types)))
  colnames(Imputed_CNA_Matrix) <- colnames(CNA_Matrix)
  rownames(Imputed_CNA_Matrix) <- gencode_gene_types$gene_name

  # iterate
  for (cell in names(CNA_Segments_split)){
    CNA_cell <- CNA_Segments_split[[cell]]
    CNA_cell<-droplevels(CNA_cell)
    #print(cell)
    #Imputed_CNA_Matrix[,cell]

    # iterate chr
    #gencode_gene_types_split <- split(gencode_gene_types, gencode_gene_types$chr)
    CNA_cell_split <- split(CNA_cell, CNA_cell$Chr)

    for (CNA_event_chr in 1:length(CNA_cell_split) ){
      #print(CNA_event_chr)
      #print(cell)
      CNA_cell_split_chr  <- CNA_cell_split[[CNA_event_chr]]
      chromosome <- unique(CNA_cell_split_chr$Chr)

      if(nrow(CNA_cell_split_chr)!=0){
        gencode_gene_types_split_chr <- subset(gencode_gene_types, gencode_gene_types$chr == as.character(chromosome))

        for (CNA_event in 1:nrow(CNA_cell_split_chr)){

          CNA_event_row <- CNA_cell_split_chr[CNA_event,]
          if(!is.na(CNA_event_row$CNA_value)){
            index_start <- which(gencode_gene_types_split_chr$start >= CNA_event_row$loc.start)[1]
            index_end <- which(gencode_gene_types_split_chr$end <= CNA_event_row$loc.end)[length(which(gencode_gene_types_split_chr$end <= CNA_event_row$loc.end))]

            gene_start <- gencode_gene_types_split_chr$gene_name[index_start]
            gene_end <- gencode_gene_types_split_chr$gene_name[index_end]


            # if( CNA_event_row$CNA_value !=0){
            #   #browser()
            # }
            #browser()
            #print(cell)
            #print(CNA_event_chr)
            #print(CNA_event)
            Matrix_index_start <- which(rownames(Imputed_CNA_Matrix) ==  gene_start)
            Matrix_index_end <- which(rownames(Imputed_CNA_Matrix) ==  gene_end)

            # index is based on chr
            Imputed_CNA_Matrix[Matrix_index_start:Matrix_index_end,cell] <- CNA_event_row$CNA_value
          }
        }
      }
    }
  }
  return(Imputed_CNA_Matrix)
}
