#' @title Preprocessing of the data to apply scTransform
#'
#' @description \code{Preprocessing_scTransform} filters out low quality cells, by
#' removing cells with less than 2.000 detected genes and we remove lowly
#' expressed genes that are present in less than 10\% of the cells.
#'
#' @usage Preprocessing_scTransform(data, metadata)
#'
#' @param data matrix; Raw read counts from scRNA-seqdata with cells in
#' columns and genes in rows.
#' @param metadata data.frame; The names of the cells should be in the rownames of \code{metadata}.
#'
#' @return SeuratObject object containing the filtered \code{data} and \code{metadata}
#'
#' @export

Preprocessing_scTransform <- function(data, metadata, filter_gene){

  if (missing(data) | is.null(data)) {
    stop("Please provide an input data matrix\n")}
  if (missing(metadata) | is.null(data)) {
    stop("Please provide an input metadata table\n")}
  #if (!is.matrix(data)) {
  #  stop("data should be a matrix!")}
  #if (!is.data.frame(metadata)) {
  #  stop("metadata should be a data.frame!")}
 
  data <- as.matrix(data)
  
  	if(!is.null(metadata)){
  	  metadata <- as.data.frame(metadata)


	  if( !identical(colnames(data), rownames(metadata))){
	    stop("colnames of data should be identical to rownames of metadata!")}
	}

  # Filter cells
  data <- data[, apply(data, 2, function(x) length(which(x != 0)))> 2000]
  
  # Filter genes 
  if(filter_gene){
    genes_to_keep <- row.names(data[apply(data, 1, function(x) length(which(x!=0))/length(x) >=.1),])
  }else{
    genes_to_keep <- row.names(data[apply(data, 1, function(x) sum(x) != 0),]) # rm genes that are 0 in all cells
  }

  data <- data[genes_to_keep,]
  # Filter cells again
  data <- data[, apply(data, 2, function(x) length(which(x != 0)))> 2000]
    
  # Subset metadata to the same cells that are in data
  metadata <- metadata[colnames(data),]
  all.equal(colnames(data),rownames(metadata))

  # Create a Seurat object
  counts_se <- Seurat::CreateSeuratObject(counts = data, meta.data = metadata)

  print("keep number genes:")
  print(length(genes_to_keep))


  print("keep number cells:")
  print(ncol(data))
  return(counts_se)
}
