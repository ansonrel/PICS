#Aug 28 2020


clr_function <- function(x) {
  return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
}



#' CLR seurat
#' 
#' Use seurat normalization method for ADT counts
#'
#' @param xlist An object of class `SingleCellExperiment`
#' @param margin Normalize across features (1) or cells (2)
#' @param counts assay used to normalize
#' 
#' @return A SCE object.
#' @importFrom SingleCellExperiment assay 
#' @export 
clr.seurat <- function(x, margin = NA, counts = "counts"){
  
  if (!margin %in% c(1,2)){
    stop(paste("margin was misspecified in clr.seurat, expected \"1\" (ADTs) or \"2\" (cells) and received \"",
               margin, "\""))
  }

  d.seurat <- Seurat::as.Seurat(x,
                                counts = "counts")
  
  d.seurat <- Seurat::NormalizeData(d.seurat, verbose = FALSE, 
                                    normalization.method = "CLR", 
                                    margin = margin)
  
  adt.data <- Seurat::GetAssayData(d.seurat, assay = "originalexp")
  idx <- match(colnames(adt.data),colnames(x))
  rownames(adt.data) <- rownames(x)
  
  if (margin == 1) {
    assay(x, "clr_adts") <- adt.data[,idx]
  } else if (margin == 2) {
    assay(x, "clr_cells") <- adt.data[,idx]
  } 
  
  return(x)
}


#' Use seurat normalization method for ADT counts
#' across adts (margin = 1)
clr.seurat_adts <- function(xlist, counts = "counts"){
  return(clr.seurat(xlist, 1, counts))
}
#' Use seurat normalization method for ADT counts
#' across cells (margin = 2)
clr.seurat_cells <- function(xlist, counts = "counts"){
  return(clr.seurat(xlist, 2, counts))
}

#' ScaleData Seurat
#' 
#' CLR normalize and Scale data as in Seurat pipeline
#'
#' @param xlist An object of class `SingleCellExperiment`
#' @param margin Normalize across features (1) or cells (2)
#' @param center Whether to center the data
#' @param scale Whether to scale the data
#' @param batch column in colData specifying the batch
#' 
#' @return A SCE object.
#' @importFrom SingleCellExperiment assay 
#' @export 
scaledata.seurat <- function(xlist, margin = 2, center = TRUE, scale = TRUE,
                             batch = "batch"){
  
  x <- xlist[[1]]
  
  d.seurat <- Seurat::as.Seurat(x,
                                counts = "counts",
                                data = NULL,
                                assay = "ADT")
  
  d.seurat <- Seurat::NormalizeData(d.seurat, verbose = FALSE, 
                                    normalization.method = "CLR", 
                                    assay = "ADT",
                                    margin = margin)
  
  d.seurat <- Seurat::ScaleData(d.seurat, 
                                assay = "ADT", verbose = FALSE,
                                do.scale = scale,
                                do.center = center,
                                split.by = batch
  )
  
  adt.data <- Seurat::GetAssayData(d.seurat, slot = "scale.data", assay = "ADT")
  
  idx <- match(colnames(adt.data),colnames(x))
  rownames(adt.data) <- rownames(x)
  assay(x, "scaleData") <- adt.data[,idx]
  
  return(x)
}

#' rescaleBatches batchelor
#' 
#' Modified rescaleBatches() function to use CLR normalizaton instead of
#' logcounts
#'
#' @param xlist An object of class `SingleCellExperiment`
#' @param margin Normalize across features (1) or cells (2)
#' @param batch column in colData specifying the batch
#' @param counts assay used to normalize
#' 
#' @return A SCE object.
#' @importFrom SingleCellExperiment assay 
#' @export 
scaledata.batchelor <- function(xlist, margin = 2, counts = "counts",
                                batch = "batch"){
  
  x <- xlist[[1]]
  
  cs <- split(seq_len(ncol(x)), x[[batch]])
  batches <- lapply(cs, function(u){
    assay(x, counts)[,u]
  })
  
  nbatches <- length(batches)
  
  averages <- vector("list", nbatches)
  for (b in seq_along(batches)) {
    curbatch <- batches[[b]]
    averages[[b]] <- rowMeans(curbatch)
  }
  
  # Defining the reference.
  reference <- do.call(pmin, averages)
  for (b in seq_along(batches)) {
    rescale <- reference / averages[[b]] 
    rescale[!is.finite(rescale)] <- 0
    batches[[b]] <- batches[[b]] * rescale
  }
  
  batch.labels <- names(batches)
  ncells.per.batch <- vapply(batches, ncol, FUN.VALUE=0L)
  batch.names <- rep(batch.labels, ncells.per.batch)
  
  #join all batches and add to SCE
  corrected <- do.call(cbind, batches)
  assay(x, "norm") <- corrected
  
  return(x)
}

#' dsb normalization
#' 
#' background normalization with dsb
#'
#' @param xlist An object of class `SingleCellExperiment`
#' @param batch column in colData specifying the batch
#' @param iso character to grep isotype markers
#' 
#' @return A SCE object.
#' @importFrom SingleCellExperiment assay 
#' @export 
normalize.dsb <- function(xlist, iso = "type",
                          batch = "batch"){
  
  pos <- xlist[[1]]
  neg <- xlist[[2]]
  
  #isotypes <- rownames(pos)[grepl(iso,rownames(pos))]
  isotypes <- rownames(pos)[rowData(kot)$is.isotype]
  
  cs <- split(seq_len(ncol(pos)), pos[[batch]])
  cs_neg <- split(seq_len(ncol(neg)), neg[[batch]])
  
  pos_batches <- lapply(cs, function(u){
    assay(pos, "counts")[,u]
  })
  
  neg_batches <- lapply(cs_neg, function(u){
    assay(neg, "counts")[,u]
  })
  
  nbatches <- length(pos_batches)
  norm_mats <- vector("list", nbatches)
  
  for (b in seq_along(pos_batches)) {
    pos_b <- pos_batches[[b]]
    neg_b <- neg_batches[[b]]
    
    #apply dsb function
    norm_mats[[b]] <- dsb::DSBNormalizeProtein(cell_protein_matrix = pos_b,
                                               empty_drop_matrix = neg_b, 
                                               use.isotype.control = TRUE, 
                                               isotype.control.name.vec = isotypes,
                                               define.pseudocount = TRUE, 
                                               pseudocount.use = 1)
  }
  
  
  #join all batches and add to SCE
  corrected <- do.call(cbind, norm_mats)
  idx <- match(colnames(pos), colnames(corrected))
  #corrected <- corrected
  
  assay(pos, "dsb") <- corrected[,idx]
  
  return(pos)
}

#' log norm count normalization
#' 
#' naive normalization by taking the log over normalized counts
#'
#' @param xlist An object of class `SingleCellExperiment`
#' 
#' @return A SCE object.
logNormCounts <- function(x) {

  save_logcounts <- NULL
  # to make sure we don't interfere with methods relying on values stored
  # in "logcounts" assay
  if (any(names(assays(x)) == "logcounts")) {
    save_logcounts <- logcounts(x)
  }
  
  x <- scater::logNormCounts(x)
  assay(x, "logNormCounts") <- as.matrix(logcounts(x))
  
  logcounts(x) <- save_logcounts
  
  return(x)
}