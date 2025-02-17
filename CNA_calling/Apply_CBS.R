#' @title Apply CBS algorithm
#'
#' @usage DNAcopy_NetCNA(NetCNA, gencode, alpha = 0.01, undo.SD = 25,
#' smooth.region = 10, outlier.SD.scale=4,smooth.SD.scale=2, trim=0.025)
#'
#' @param NetCNA; data frame with relative smoothed CNA scores from all input cells
#' @param gencode; gene annotation file
#' @param alpha; significance levels for the test to accept change-points.
#' @param undo.SD; the number of SDs between means to keep a split if undo.splits="sdundo".
#' @param smooth.region; number of points to consider on the left and the right of a point to detect it as an outlier. (default=10)
#' @param outlier.SD.scale; the number of SDs away from the nearest point in the smoothing region to call a point an outlier.
#' @param smooth.SD.scale; the number of SDs from the median in the smoothing region where a smoothed point is positioned.
#' @param trim; proportion of data to be trimmed for variance calculation for smoothing outliers and undoing splits based on SD.
#'
#' @return data frame with relative smoothed CNA scores from all input cells
#'
#' @description Apply CBS algorithm
#'
#' @export
#'

DNAcopy_NetCNA <- function(NetCNA, gencode, alpha = 0.01, undo.SD = 25, smooth.region = 10, outlier.SD.scale=4,smooth.SD.scale=2, trim=0.025){

  NetCNA <- gffsort(data = NetCNA,
                      symbol = "hgnc_symbol",
                      target = "hgnc_symbol",
                      gencode = gencode)

        data <- CNA(NetCNA,
                    chromosomeorder$chr,
                    chromosomeorder$start,
                    data.type = "logratio",
                    presorted = TRUE) # warning: chr11 two genes with same start positions
        data <- smooth.CNA(data,smooth.region = smooth.region, outlier.SD.scale = outlier.SD.scale,
                           smooth.SD.scale=smooth.SD.scale ,trim=trim)

        data <- segment(data,
                        verbose = F,
                        alpha = alpha,
                        undo.splits = "sdundo",
                        undo.SD = undo.SD,
                        min.width = 5)
    return(data)
}

## Function to sort data frames according to the genomic position of the genes
#
# ARGS:
# data = data matrix
# names = gene names (default the row names)
# symbol = either hgnc_symbol or ensemble_gene_id
# target = target identifier
# gencode = gencode object
gffsort <- function(data,names=row.names(data),symbol,target,gencode){
  if(str_detect(gencode$gene_id[1], "ENSG")){n_autosome=22}else if(str_detect(gencode$gene_id[1], "ENSMUSG")){n_autosome=19}
  if(symbol=="hgnc_symbol"){
    data <- as.matrix(data)
    print(dim(gencode))
    genes.have <- intersect(names,gencode$gene_name)
    data <- data[genes.have,]
    gencode <- gencode[match(row.names(data),gencode$gene_name),]
    gencode$chr <- factor(gsub("chr","",gencode$seqid),levels=c(seq(1,n_autosome),"X","Y"))
    filter <- factor(c(seq(1,n_autosome),"X","Y"),levels=c(seq(1,n_autosome),"X","Y"))
    gencode <- gencode[which(gencode$chr %in% filter),]
    gencode <- gencode[order(gencode$chr,gencode$start),]
    gencode$ensembl <- sapply(gencode$gene_id,function(x) strsplit(x,"\\.")[[1]][1])
    data <- data[match(gencode$gene_name,row.names(data)),]
    if(target=="ensembl_gene_id"){
      row.names(data) <- gencode$ensembl
    } else {
      row.names(data) <- gencode$gene_name
    }
    chromosomeorder <- data.frame(chr=gencode$chr,
                                  start=gencode$start,
                                  end=gencode$end,
                                  ensembl=gencode$ensembl,
                                  hgnc=gencode$gene_name)
    assign("chromosomeorder",chromosomeorder,.GlobalEnv)
  } else if(symbol=="ensembl_gene_id"){
    data <- as.matrix(data)
    gencode$ensembl <- sapply(gencode$gene_id,function(x) strsplit(x,"\\.")[[1]][1])
    genes.have <- intersect(names,gencode$ensembl)
    data <- data[genes.have,]
    gencode <- gencode[match(row.names(data),gencode$ensembl),]
    gencode$chr <- factor(gsub("chr","",gencode$seqid),levels=c(seq(1,n_autosome),"X","Y"))
    filter <- factor(c(seq(1,n_autosome),"X","Y"),levels=c(seq(1,n_autosome),"X","Y"))
    gencode <- gencode[which(gencode$chr %in% filter),]
    gencode <- gencode[order(gencode$chr,gencode$start),]
    data <- data[match(gencode$ensembl,row.names(data)),]
    # row.names(data) <- gencode$ensembl
    if(target=="ensembl_gene_id"){
      row.names(data) <- gencode$ensembl
    } else {
      row.names(data) <- gencode$gene_name
    }
    chromosomeorder <- data.frame(chr=gencode$chr,
                                  start=gencode$start,
                                  end=gencode$end,
                                  ensembl=gencode$ensembl,
                                  hgnc=gencode$gene_name)
    assign("chromosomeorder",chromosomeorder,.GlobalEnv)
  }
  return(data)
}

