#' @title Inference of copy number alterations
#'
#' @usage CNA_caller(data, metadata, gencode, window_size, covariates = NULL,
#' plot_metadata="covariates",normalizedValues = NULL, output_dir="output/",  Merge=TRUE, n_permut = 10)
#'
#' @param data matrix; Raw read counts from scRNA-seqdata with cells in
#' columns and genes in rows.
#' @param metadata data.frame; Metadata information containing the Donor and
#' celltype information about the cells. The cell names should be in the
#' rownames of \code{metadata}.
#' @param gencode gencode
#' @param window_size window_size
#' @param covariates design matrix; Sample-level covariates. If missing, it will
#' contain only an intercept. See \code{X} parameter of \code{zinbwave}.
#' @param plot_metadata plot_metadata
#' @param normalizedValues if normalized values are already computed, insert it here
#' @param output_dir string; define the output directory
#' @param Merge Merge
#' @param n_permut Number of permutations for CNA background profile
#'
#' @return CNA Profile per cell
#'
#' @description
#' A more detailed description of each input is provided below:
#'
#' The \code{data}:
#'
#'            SRR5959721 SRR5959722 SRR5959726 SRR5959727 SRR5959732
#' DVL1                  0          0          0          0          0
#' MXRA8                 1         16          1          0          1
#' AURKAIP1              0          0          0          0          1
#' CCNL2                 0         19         12          1          1
#' RP4-758J18.2          0          0          0          0          0
#' MRPL20                0          0          0         43          0
#' ...
#'
#' The \code{metadata}:
#'
#'            Donor Celltype Treatment Age
#' SRR5959721    P1       X       post    44
#' SRR5959722    P1       X       post    44
#' SRR5959726    P1       X       post    44
#' SRR5959727    P1       X       post    44
#' SRR5959732    P1       X       post    44
#' SRR5959733    P1       X       post    44
#'
#' @export
#'

CNA_caller <- function(data, metadata, gencode, window_size, covariates = NULL, plot_metadata="covariates",
                       normalizedValues = NULL, output_dir="output/", Merge=TRUE, UMI=TRUE, gencode_gene_types, filter_gene=TRUE, n_permut = 10){
message("output directory:")
message(output_dir)


if (!file.exists(output_dir)){
    dir.create(file.path(output_dir))
}

if(plot_metadata == "covariates"){plot_metadata <- covariates}

  if(is.null(normalizedValues)){

    if (UMI){
        message("Preprocess data")
	    data <- Preprocessing_scTransform(data, metadata, filter_gene)

	    message("Step1: Normalizing data using scTransform")
	    # scTransform
	    data <- Seurat::SCTransform(data, verbose = T, vars.to.regress = covariates, return.only.var.genes = FALSE, vst.flavor = "v2")
	    saveRDS(data, file = paste0(output_dir,"normalized_data.rds"))
	}else{
	    message("Preprocess data")
	    counts_se <- Preprocessing_ZINBWaVE(data, metadata)

	    message("Step1: Normalizing data using ZINBWaVE")
	    # Zinbwave
	    if(is.null(covariates)){
	        counts_zinb_K0 <- zinbFit(Y = counts_se, K = 0, verbose=T,BPPARAM = BiocParallel::SerialParam())
	    }else{
	        counts_zinb_K0 <- zinbFit(Y = counts_se, X = covariates, K = 0, verbose=T,BPPARAM = BiocParallel::SerialParam())
        }
	    data <- zinbwave(counts_se, fitted_model = counts_zinb_K0, normalizedValues=TRUE, residuals=TRUE, verbose=T)
	    saveRDS(data, file = paste0(output_dir,"normalized_data.rds"))
	}

    }else if(is.object(normalizedValues)){
        data <- normalizedValues
    }else{
        stop("'Seurat' object including scTransform-normalized data is missing!")
    }

  # Fixed moving window
  message("Step2: Fixed moving window")

  rm(metadata)
  gc()
  mallinfo::malloc.trim()
  browser()

  # check input and compute number of detected genes
  if (UMI){
    if( is.null(data@assays$SCT@scale.data) ){stop("assay 'normalizedValues' is missing in 'SingleCellExperiment' object!")}
    if( is.null(data@meta.data) ){stop("assay 'colData' is missing in 'SingleCellExperiment' object!")}
    input_data <- data@assays$SCT@scale.data
  }else{
    if( is.null(assay(data,'normalizedValues')) ){stop("assay 'normalizedValues' is missing in 'SingleCellExperiment' object!")}
    if( is.null(data@colData) ){stop("assay 'colData' is missing in 'SingleCellExperiment' object!")}
    input_data <- assay(data,'normalizedValues')
  }

    if (missing(gencode) | is.null(gencode)) {
        stop("Please provide gencode")
    }

  CNAs_per_gene <- FixedMovingWindow(data = input_data,
                                     gencode = gencode,
                                     detected_genes = detected_genes,
                                     window_size=151)
  saveRDS(CNAs_per_gene, file = paste0(output_dir,"CNAs_per_gene.rds"))


  #remove cells with no smoothed CNA scores or just one CNA score for one gene (otherwise there is an error in CBS algorithm)
  n_cna_per_cell <- apply(CNAs_per_gene,2, function(cell) length(which(cell !=0)))
  cells_keep<-which(n_cna_per_cell>1)
  CNAs_per_gene<-CNAs_per_gene[,cells_keep]



  # CBS
  message("Step3: CNA boundary detection")
  CBS_CNAs_per_gene <- DNAcopy_NetCNA(CNAs_per_gene, gencode, undo.SD = 20)
  saveRDS(CBS_CNAs_per_gene, file=paste0(output_dir,'CBS_CNAs_per_gene.rds'))

  rm(CBS_CNAs_per_gene); gc(); mallinfo::malloc.trim() # empty memory

    # Permutation
    for (i in 1:n_permut){

      gCNA_per <- t(apply(input_data, 1, sample)) #permut for each gene between cells
      colnames(gCNA_per) <- colnames(input_data)

      gCNA_per<-gCNA_per[,cells_keep]#remove cells with no CNAs in the true CNA calls
      saveRDS(gCNA_per, file = paste0(output_dir,"gCNA_per_tmp",i,".rds"))
      CNAs_per_gene_permut_tmp <- FixedMovingWindow(data = gCNA_per,
                                                                          gencode = gencode,
                                                                          detected_genes = detected_genes,
                                                                          genes_per_window = window_size)
      saveRDS(CNAs_per_gene_permut_tmp, file = paste0(output_dir,"CNAs_per_gene_permut_tmp",i,".rds"))
      rm(CNAs_per_gene_permut_tmp)
      gc()
      mallinfo::malloc.trim()
    }
    # save everything in one list
    CNAs_per_gene_permut <-list()
    gCNA_per <-list()
    for (i in 1:10){
    CNAs_per_gene_permut[[i]] <- readRDS(file = paste0(output_dir,"CNAs_per_gene_permut_tmp",i,".rds"))
    gCNA_per[[i]] <- readRDS(file = paste0(output_dir,"gCNA_per_tmp",i,".rds"))
    }
    saveRDS(CNAs_per_gene_permut, file = paste0(output_dir,"CNAs_per_gene_permut.rds"))
    saveRDS(gCNA_per, file = paste0(output_dir,"Residuals_permuted.rds"))

    #remove tmp files
    for (i in 1:10){
      if (file.exists(paste0(output_dir,"CNAs_per_gene_permut_tmp",i,".rds"))) {
        unlink(paste0(output_dir,"CNAs_per_gene_permut_tmp",i,".rds"))
        print("File is deleted..")
      } else{
       print("File not exists..")
      }
      if (file.exists(paste0(output_dir,"gCNA_per_tmp",i,".rds"))) {
        unlink(paste0(output_dir,"gCNA_per_tmp",i,".rds"))
        print("File is deleted..")
      } else{
       print("File not exists..")
      }
    }

    rm(CNAs_per_gene);rm(gCNA_per); gc(); mallinfo::malloc.trim() # empty memory

    # CBS
    CBS_CNAs_per_gene_permut<-list()
    for (i in 1:10){
      CBS_CNAs_per_gene_permut_tmp <- DNAcopy_NetCNA(CNAs_per_gene_permut[[i]], gencode, undo.SD = 20)
      saveRDS(CBS_CNAs_per_gene_permut_tmp, file = paste0(output_dir,"CBS_CNAs_per_gene_permut_tmp",i,".rds"))
      rm(CBS_CNAs_per_gene_permut_tmp); gc(); mallinfo::malloc.trim() # empty memory
    }
    # save everything in one list
    for (i in 1:10){
    CBS_CNAs_per_gene_permut[[i]] <- readRDS(file = paste0(output_dir,"CBS_CNAs_per_gene_permut_tmp",i,".rds"))
    }
    saveRDS(CBS_CNAs_per_gene_permut, file=paste0(output_dir,'CBS_CNAs_per_gene_permut.rds'))

    #remove tmp files
    for (i in 1:10){
      if (file.exists(paste0(output_dir,"CBS_CNAs_per_gene_permut_tmp",i,".rds"))) {
        unlink(paste0(output_dir,"CBS_CNAs_per_gene_permut_tmp",i,".rds"))
        print("File is deleted..")
      } else{
       print("File not exists..")
      }
    }

    rm(CNAs_per_gene_permut);rm(CBS_CNAs_per_gene_permut); rm(data); gc(); mallinfo::malloc.trim() # empty memory


  CNAs_per_gene <- readRDS(file=paste0(output_dir,'CNAs_per_gene.rds'))
  CBS_CNAs_per_gene <- readRDS(file=paste0(output_dir,'CBS_CNAs_per_gene.rds'))

  ### Create Matrix from CBS output
  CNAs_Matrix <- CreateMatrix_fromLargeDNAcopy(CBS_CNAs_per_gene, CNAs_per_gene, gencode)
  saveRDS(CNAs_Matrix,file=paste0(output_dir,'CNAs_Matrix_not_merged.rds'))

  rm(CNAs_per_gene); gc(); mallinfo::malloc.trim() # empty memory



  # merge
  if (Merge){
    CNAs_Matrix_wo_threshold <- Merge_sequential_CNAs(CNAs_Matrix,gencode,Segments=TRUE)
    saveRDS(CNAs_Matrix_wo_threshold,file=paste0(output_dir,'CNAs_Matrix_wo_threshold.rds'))
  }else{
    saveRDS(CNAs_Matrix,file=paste0(output_dir,'CNAs_Matrix_wo_threshold_not_merged.rds'))
  }

      ########### Compute p-values
      CBS_CNAs_per_gene_permut <- readRDS(file=paste0(output_dir,'CBS_CNAs_per_gene_permut.rds'))
      CBS_CNAs_per_gene <- readRDS(file=paste0(output_dir,'CBS_CNAs_per_gene.rds'))

      message("Compute p-values")
      data <- readRDS(file = paste0(output_dir,"normalized_data.rds"))

      if(UMI){
        metadata<-data@meta.data
        counts <- data@assays$RNA@counts
      }else{
        counts <- assay(data,'counts')
        metadata<-data@colData
      }

      CBS_pvalues<-compute_pvalues(CBS=CBS_CNAs_per_gene,CBS_permut=CBS_CNAs_per_gene_permut,counts=counts, output_dir=output_dir)
      saveRDS(CBS_pvalues,file=paste0(output_dir,'CBS_pvalues.rds'))

      #CBS_pvalues_filtered <- rm_chromosomes_CBS_output(CBS_pvalues, output_dir = output_dir, name="CBS")
      #saveRDS(CBS_pvalues_filtered,file=paste0(output_dir,'CBS_pvalues_filtered.rds'))

      CNAs_per_gene <- readRDS( file = paste0(output_dir,"CNAs_per_gene.rds"))

      CBS_pvalues_filtered <- rm_chromosomes_smoothed_matrix(CNAs_per_gene , CBS_pvalues, metadata = metadata, gencode_gene_types = gencode_gene_types, output_dir = output_dir, name="CBS")
      saveRDS(CBS_pvalues_filtered,file=paste0(output_dir,'CBS_pvalues_filtered.rds'))



  ############ Imputation

  CNA_Segments <- CNAs_Matrix[[2]]
  CNA_Matrix <- CNAs_Matrix[[1]]

  CNA_Matrix_Imputed<-list()
  CNA_Matrix_Imputed[[1]] <- Imputation_CNA(CNA_Matrix,CNA_Segments,gencode_gene_types)
  CNA_Matrix_Imputed[[2]] <- CNA_Segments

  saveRDS(CNA_Matrix_Imputed,file=paste0(output_dir,'CNA_Matrix_Imputed.rds'))

  ######## Remove Chromosomes where I don't have enough CNAs across the cells
  data <- readRDS(file = paste0(output_dir,"normalized_data.rds"))
  if(UMI){metadata<-data@meta.data}else{metadata<-data@colData}

  CNA_Matrix_Imputed_filtered <- rm_chromosomes(CNA_Imputed = CNA_Matrix_Imputed, metadata = metadata, gencode_gene_types = gencode_gene_types,output_dir=output_dir,name="Imputed")
  saveRDS(CNA_Matrix_Imputed_filtered,file=paste0(output_dir,'CNA_Matrix_Imputed_filtered.rds'))
  rm(CNA_Matrix_Imputed);rm(CNA_Matrix_Imputed_filtered)
  gc()
  mallinfo::malloc.trim()

  CNA_Matrix_filtered <- rm_chromosomes(CNA_Imputed = CNAs_Matrix, metadata = metadata, gencode_gene_types = gencode_gene_types,output_dir=output_dir,name="")
  saveRDS(CNA_Matrix_filtered,file=paste0(output_dir,'CNA_Matrix_filtered.rds'))
  # Segments are already subsetted in this function

  rm(CNAs_Matrix);rm(CNA_Matrix_filtered)
  gc()
  mallinfo::malloc.trim()

  ########## VISUALIZATION
  message("Start plotting")
  print("Start plotting")
  data <- readRDS(file = paste0(output_dir,"normalized_data.rds"))

  if(UMI){metadata<-data@meta.data}else{metadata<-data@colData}

  gc()
  mallinfo::malloc.trim()

  # read
  CNAs_Matrix<- readRDS(file=paste0(output_dir,'CNAs_Matrix.rds'))
  CNAs_Matrix <- CNAs_Matrix[[1]]
  all.equal(rownames(metadata),colnames(CNAs_Matrix))
  metadata <- droplevels(metadata)

  ### Heatmap
  create_Heatmaps(CNA_profile = CNAs_Matrix,
                  metadata = metadata[colnames(CNAs_Matrix),plot_metadata],
                  gencode = gencode,
                  path = paste0("plots/",output_dir),
                  file_name = "Heatmap_RNA2CNA")

  gc()
  mallinfo::malloc.trim()

    ##### Permutation Plots
    CNAs_per_gene_permut <-  readRDS(file = paste0(output_dir,"CNAs_per_gene_permut.rds"))

    Mutation_Frequency_per_Cell_Amp <- apply(CNAs_per_gene_permut[[1]][,rownames(metadata)], 2, function(row) length( which(row !=0 & row > 0))/(length(na.omit(row) ))*100 )
    Mutation_Frequency_per_Cell_Del <- apply(CNAs_per_gene_permut[[1]][,rownames(metadata)], 2, function(row) length( which(row !=0 & row < 0))/(length(na.omit(row) ))*100 )

    Mutation_Frequency_per_Cell <- data.frame('Amp'=Mutation_Frequency_per_Cell_Amp, 'Del'= Mutation_Frequency_per_Cell_Del)

    if(UMI){
    detected_genes <- DetectedGenes(data@assays$RNA@counts, Cell_Ref=NULL)
    }else{
    detected_genes <- DetectedGenes(assay(data,'counts'),Cell_Ref=NULL)
    }
    all.equal(names(detected_genes),rownames(Mutation_Frequency_per_Cell))
    Mutation_Frequency_per_Cell$detected_genes <- detected_genes

    png(paste0("plots/",output_dir,'AfterRNA2CNA_Percentage_CNAs_Del_Amp_per_Cell_vs_detected_genes.png'),width = 800,height = 500)
    p <- ggplot2::ggplot(reshape::melt(Mutation_Frequency_per_Cell,id.vars = c("detected_genes"))) + geom_point(aes(x=.data$detected_genes,y=.data$value))+ theme_readable(base_size = 14)+
    facet_wrap(.~variable)+xlab("Number of detected genes")+ylab("Percentage of CNAs")
    print(p)
    dev.off()

    rm(data)
    gc()
    mallinfo::malloc.trim()

}



