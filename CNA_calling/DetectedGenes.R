#' @title Calculation of the number of detected genes per cell
#'
#' @param test_data matrix; Raw read counts from scRNA-seqdata with cells in
#' columns and genes in rows.
#' @param Cell_Ref; Cell_Ref
#'
#' @return Number of detected genes per cell
#'
#' @description Compute number of detected genes per cell
#'
#' @export

DetectedGenes <- function(test_data, Cell_Ref){
  if(is.null(Cell_Ref)){
    number_detected_genes <- apply(test_data, 2, function(x) length(which(x != 0)))
  }else{
    number_detected_genes <- list()
    for(donor in 1:length(test_data)){
      number_detected_genes[[donor]] <- apply(test_data[[donor]], 2, function(x) length(which(x != 0)))
    }
    number_detected_genes <- unlist(number_detected_genes)
  }
  cat("Summary of number of detected genes:\n")
  print(summary(number_detected_genes))
  return(number_detected_genes)
}
