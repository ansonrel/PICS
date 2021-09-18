
# function used to remove isotype ADTs from the analysis
get_markers <- function(x) {
  if (any(names(rowData(x)) == "is.isotype")) {
    # Select out the isotypes
    markers <- rownames(x)[!rowData(x)$is.isotype] 
  }
  
  return(markers)
}

# function used to get the name of the normalization method used
get_assay <- function(x) {
  
  known_assays <- c("clr_adts", "clr_cells", "scaleData", "dsb", "logNormCounts")
  found_assay <- names(assays(x)) %in% known_assays
  
  if( any(found_assay) ) {
    asa <- names(assays(x))[found_assay]
  } else {
    stop("No assay found")
  }
  
  return(asa)
}

# function used to save sce after normalization for low dimensionality check
save_sce <- function(x, output_sce, clustmethod) {
  # Check which normalization technique was used in order to store the sce
  if (any(names(assays(x)) != "counts")) {
    dir.create(output_sce, showWarnings = FALSE)
    for (a in names(assays(x))[-1]) {
      dir.create(paste0(output_sce, "/", x@metadata$name, "_sce"), showWarnings = FALSE)
      
      dir.name <- paste0(output_sce, "/", x@metadata$name, "_sce/", clustmethod)
      dir.create(dir.name, showWarnings = FALSE)
      saveRDS(x, paste0(dir.name,"/", a, ".rds"))
    }
  } else {
    stop(paste("Could not find the normalization technique used in", x@metadata$name))
  }
  
}


# functions used to calculate the normalization score based on negative and 
# positive control signals

get_negScore <- function(d) {
  # Shift values toward zero
  d <- d - min(d)
  # Scale so that distribution is between 0 and 1
  d <- d/max(d)
  
  # Returns 10 if median or sd is 0
  return( 1/(abs(median(d)) * sd(d) + .1) )
}

get_posScore <- function(d_neg, d_pos) {
  
  # Shift values toward zero
  d_neg <- d_neg - min(d_neg)
  d_pos <- d_pos - min(d_pos)
  
  # Scale so that both distributions are between 0 and 1
  d_neg <- d_neg/max(d_neg)
  d_pos <- d_pos/max(d_pos)
  
  # Compare
  sd_neg <- sd(d_neg)
  sd_pos <- sd(d_pos)
  
  t.res <- t.test(d_neg, d_pos, sigma.x=sd_neg, sigma.y=sd_pos, alternative="two.sided") 
  #would make more sens to use "greater" but the results don't change so not sure what's up
  
  return( abs(t.res$statistic) )
}

get_normScore <- function(d_neg, d_pos) {
  assertthat::not_empty(d_neg)
  assertthat::not_empty(d_pos)
  
  
  if (is.null(dim(d_neg)) & is.null(dim(d_pos))) {
    # both object only have a single element
    negScore <- get_negScore(d_neg)
    posScore <- get_posScore(d_neg, d_pos)
    
    
  } else if (is.null(dim(d_neg))) {
    # there is only one isotype
    negScore <- get_negScore(d_neg)
    
    posScore.list <- lapply(d_pos, function(x) {
      get_posScore(d_neg, x)
    })
    
    posScore <- mean(unlist(posScore.list))
    
  } else if (is.null(dim(d_pos))) {
    # there is only one positive control
    negScore.list <- lapply(rownames(d_neg), function(x) {
      get_negScore(d_neg[x,])
    })
    
    posScore.list <- lapply(d_neg, function(x) {
      get_posScore(x, d_pos)
    })
    
    negScore <- mean(unlist(negScore.list))
    posScore <- mean(unlist(posScore.list))
    
  } else {
    # there are mutliple isotypes and positive controls
    
    comb <- expand.grid(1:dim(d_neg)[1], 1:dim(d_pos)[1])
    
    negScore.list <- lapply(rownames(d_neg), function(x) {
      get_negScore(d_neg[x,])
    })
    
    posScore.list <- apply(comb, MARGIN = 1, FUN=function(x) {
      get_posScore( d_neg[x[1],], d_pos[x[2],])
    })
    
    negScore <- mean(unlist(negScore.list))
    posScore <- mean(unlist(posScore.list))
  }
  
  #TODO: return the score of every combination of isotypes and positive controls
  # instead of aggregating them using mean()
  return(negScore * posScore)
}

# function used to save negative and positive controls to later assess the
# performance of normalization
saveControlScore <- function(x) {
  
  score <- 0 
  #if the dataset has both isotypes and positive controls, score returned 
  # should be larger than 0
  
  if (!any(names(rowData(x)) == "is.isotype") | all(!rowData(x)$is.isotype)) {
    print(paste("Warning: No isotypes are available in", x@metadata$name, "\n the evaluation of negative and positive controls separation is skipped"))
  } else if (!any(names(rowData(x)) == "is.posControl") | all(!rowData(x)$is.posControl)) {
    print(paste("Warning: No ipositive controls are available in", x@metadata$name, "\n the evaluation of negative and positive controls separation is skipped"))
  } else {
    
    asa <- get_assay(x)
    
    names_iso <- rownames(x)[rowData(x)$is.isotype]
    names_pos <- rownames(x)[rowData(x)$is.posControl]
    
    iso <- assays(x)[[asa]][names_iso,]
    pos <- assays(x)[[asa]][names_pos,]
    
    score <- get_normScore(iso, pos) 
    
    control.list <- list("name"=x@metadata$name, "method"=asa, "score"=score)
    
    dir.create(paste0("pipelineOutput/controlScore/",x@metadata$name), showWarnings = FALSE)
    saveRDS(control.list, fil=paste0("pipelineOutput/controlScore/",x@metadata$name,"/controlScore.Rds"))
    
  }
  
  return(score)
}

picsScore <- function(x, assaynam = NULL) {
  
  if(is.null(assaynam)){
    asa <- get_assay(x)
  } else {
    asa <- assaynam
  }
  
  #cat(paste0("PICS score on \t", x@metadata$name, "\t", asa, "\t ..."))
  
  score <- 0
  #if the dataset has both isotypes and positive controls, score returned 
  # should be larger than 0
  
  if (!any(names(rowData(x)) == "is.isotype") | all(!rowData(x)$is.isotype)) {
    cat(paste("\n Warning: No isotypes available in \t", x@metadata$name, "\n PICS score computation skipped \n"))
  } else if (!any(names(rowData(x)) == "is.posControl") | all(!rowData(x)$is.posControl)) {
    cat(paste("\n Warning: No positive controls available in \t", x@metadata$name, "\n PICS score computation skipped \n"))
  } else {
    
    names_iso <- rownames(x)[rowData(x)$is.isotype]
    names_pos <- rownames(x)[rowData(x)$is.posControl]
    
    iso <- assays(x)[[asa]][names_iso,]
    pos <- assays(x)[[asa]][names_pos,]
    
    score <- get_normScore(iso, pos) 
    
    #cat(paste("\t", round(score, 3),"\t Done! \n"))
  }
  
  return(score)
}