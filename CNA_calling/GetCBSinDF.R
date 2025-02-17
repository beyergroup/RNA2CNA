#' @title Function to create Matrix from a large DNAcopy object
#'
#' @usage CreateMatrix_fromLargeDNAcopy(DNAcopy, CNAProf, gencode)
#'
#' @param DNAcopy; large DNAcopy object
#' @param CNAProf; matrix/dataframe with relative CNA profiles from which CNAs where called
#' @param gencode; gene annotation file
#'
#' @return CNA profile matrix for every gene
#'
#' @description Create a CNA profile matrix for every gene based on the CBS output
#'
#' @export
#'
CreateMatrix_fromLargeDNAcopy <- function(DNAcopy, CNAProf, gencode){
    CNAProf <- genome_sort(CNAProf, gencode)
    CBSList <-  split(DNAcopy$output, factor(DNAcopy$output$ID, levels = unique(DNAcopy$output$ID)))
    CBSList <- lapply(CBSList, function(x){
        unname(unlist(apply(x, 1, function(y){
            as.numeric(rep(y[[6]], y[[5]]))
        })))
    })
    CBSProf <- t(do.call(rbind, lapply(as.list(1:length(CBSList)), function(i){
        c(CBSList[[i]], rep(NA, sum(is.na(CNAProf[,i]))))[order(as.numeric(c(which(!is.na(CNAProf[,i])), which(is.na(CNAProf[,i])))))]
    })))
    row.names(CBSProf) <- row.names(CNAProf)
    colnames(CBSProf) <- colnames(CNAProf)
    return(CBSProf)
}


#' Function to sort data frames according to the genomic position of the genes
#'
#' @param data; data matrix; rows -> genes_hgnc
#' @param gencode; gencode object
#'
#' @return matrix of genes sorted by chromosomal order
#' @export
#'
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
