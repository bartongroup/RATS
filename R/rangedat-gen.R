#===============================================================================
#' Generate artificial read-count data for simulated 2-transcript genes, 
#' systematically covering a broad range of proportions, magnitudes and fold-changes.
#' 
#' @return sleuth-like object. Contains \code{"kal"} ans \code{"bootstrap"} with the simulated data and the non-sleuth 
#'         \code{"rangesim"} with the simulation details.
#'
#' @param propfrom, propto, propby: Range and step of proportions for 1st transcript. Sibling transcript complemented from these.
#' @param magfrom, magto, magby: Range of count magnitudes, in powers of 10.
#' @param foldfrom, foldto, foldby: Range of transcript ratio fold changes, in powers of 2.
#' @return list containing: \code{data} a minimal sleuth-like object.
#'                          \code{anno} a dataframe matching \code{target_id} to \code{parent_id}.
#'                          \code{sim} a list of dataframes with the simulation parameters per transcript per condition.
#'
#' @export
rangedat_gen <- function(propfrom=0, propto=1, propby=0.05,
                         magfrom=-5, magto=5, magby=1,
                         foldfrom=-7, foldto=7, foldby=1) {
  # Range of proportions. (21 levels)
  # Although the relationship between transcript is symmetric and the ranges 0-0.5 and 0.5-1 are equivalent, 
  # they are not redundant: The two sibling transcripts receive different treatment when fold-changes are applied.
  proportions_A_T1 <- seq(from=propfrom, to=propto, by=propby)
  proportions_A_T2 <- 1-proportions_A_T1
  
  # Range of magnitudes. (11 levels)
  magnitudes <- seq(from=-magfrom, to=magto, by=magby)
  
  # Range of ratio fold-changes. (15 levels)
  folds <- seq(from=foldfrom, to=foldto, by=foldby)
  
  
  # Combine proportions and magnitudes. (3465 rows)
  # Start by griding the values for easier visualization and naming, instead of compbing the outer products directly.
  a1 <- expand.grid(proportions_A_T1, magnitudes, folds)  # condition A, transcript 1, to remain the same
  a2 <- expand.grid(proportions_A_T2, magnitudes, folds)  # conditions A, transcript 2, to remain the same
  b1 <- expand.grid(proportions_A_T1, magnitudes, folds)  # condition B, transcript 1, to be decreased
  b2 <- expand.grid(proportions_A_T2, magnitudes, folds)  # condition B, transcript 2, to be increased
  
  names(a1) <- c("prop", "mag", "fold")
  names(a2) <- c("prop", "mag", "fold")
  names(b1) <- c("prop", "mag", "fold")
  names(b2) <- c("prop", "mag", "fold")
  
  
  # Calculate "counts".
  a1["est_counts"] <- a1$prop * (10 ^ a1$mag)
  a2["est_counts"] <- a2$prop * (10 ^ a2$mag)
  # Apply foldchange symmetrically to the transcripts. The total ratio change is the nominal foldchange.
  b1["est_counts"] <- b1$prop * (10 ^ b1$mag) * (1 / sqrt(2 ^ b1$fold))
  b2["est_counts"] <- b2$prop * (10 ^ b2$mag) * sqrt(2 ^ b2$fold)
  
  
  # Name the "transcripts".
  opts <- options()
  options(scipen=0, digits=3)
  
  a1["parent_id"] <- paste("p", a1$prop, "_m", a1$mag, "_x", a1$fold, sep="")
  a2["parent_id"] <- paste("p", a1$prop, "_m", a2$mag, "_x", a2$fold, sep="")
  b1["parent_id"] <- paste("p", b1$prop, "_m", b1$mag, "_x", b1$fold, sep="")
  b2["parent_id"] <- paste("p", b1$prop, "_m", b2$mag, "_x", b2$fold, sep="")
  
  a1["target_id"] <- paste(a1$parent_id, "_t1", sep="")
  a2["target_id"] <- paste(a2$parent_id, "_t2", sep="")
  b1["target_id"] <- paste(b1$parent_id, "_t1", sep="")
  b2["target_id"] <- paste(b1$parent_id, "_t2", sep="")
  
  options(scipen=opts$scipen, digits=opts$digits)
  
  
  # Make a sleuth-like object.
  sl <- list()
  sl["kal"] <- list()
  sl[["sample_to_covariates"]] <- data.frame("sample"=c("A-1", "B-1"), "condition"=c("A", "B"))
  
  sl$kal[[1]] <- list()
  sl$kal[[1]]["bootstrap"] <- list()
  sl$kal[[1]]$bootstrap[[1]] <- data.frame("target_id"=c(a1$target_id, a2$target_id), "est_counts"=c(a1$est_counts, a2$est_counts))
  
  sl$kal[[2]] <- list()
  sl$kal[[2]]["bootstrap"] <- list()
  sl$kal[[2]]$bootstrap[[1]] <- data.frame("target_id"=c(b1$target_id, b2$target_id), "est_counts"=c(b1$est_counts, b2$est_counts))
  
  # Create annotation:
  anno <- data.frame("target_id"=c(a1$target_id, a2$target_id), "parent_id"=c(a1$parent_id, a2$parent_id))
  
  # Store simulation for ease.
  results <- list("data"=sl, "anno"=anno, "sim" = list("A_t1"=a1 , "A_t2"=a2, "B_t1"=b1, "B_t2"=b2))
  
  return(results)
}