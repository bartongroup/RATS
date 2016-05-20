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
rangedat_gen <- function(propfrom=0, propto=1, propby=0.01,
                         magfrom=0, magto=3, magby=0.25,
                         foldfrom=0, foldto=7, foldby=0.1) {
  # Range of proportions. (21 levels)
  proportions_A_T1 <- seq(from=propfrom, to=propto, by=propby)
  proportions_A_T2 <- 1-proportions_A_T1
  
  # Range of magnitudes. (11 levels)
  magnitudes <- seq(from=magfrom, to=magto, by=magby)
  
  # Range of ratio fold-changes. (15 levels)
  folds <- seq(from=foldfrom, to=foldto, by=foldby)
  
  
  # Combine proportions and magnitudes. (3465 rows)
  # Start by griding the values for easier visualization and naming, instead of compbing the outer products directly.
  a1 <- expand.grid(proportions_A_T1, magnitudes, folds)  # condition A, transcript 1, to remain the same
  a2 <- expand.grid(proportions_A_T2, magnitudes, folds)  # conditions A, transcript 2, to remain the same
  names(a1) <- c("prop", "mag", "fold")
  names(a2) <- c("prop", "mag", "fold")
  # Calculate "counts".
  a1["est_counts"] <- a1$prop * (10 ^ a1$mag)
  a2["est_counts"] <- a2$prop * (10 ^ a2$mag)
  
  # Apply foldchange to proportion of A_t1 and re-scale into range 0-1: b = xa / (1 + a(x-1))
  b1 <- data.frame("prop" = (2 ^ a1$fold) * a1$prop / (1 + a1$prop * ((2 ^ a1$fold) - 1)),
                   "mag" = a1$mag,
                   "fold" = a1$fold)
  b2 <- data.frame("prop" = 1 - b1$prop, "mag" = b1$mag, "fold" = b1$fold)
  b1["est_counts"] <- b1$prop * (10 ^ b1$mag)
  b2["est_counts"] <- b2$prop * (10 ^ b2$mag)
  
  
  # Name the "transcripts".
  opts <- options()
  options(scipen=0, digits=3)
  
  a1["parent_id"] <- paste("p", a1$prop, "_m", a1$mag, "_x", a1$fold, sep="")
  a2["parent_id"] <- a1$parent_id
  b1["parent_id"] <- a1$parent_id
  b2["parent_id"] <- a1$parent_id
  
  a1["target_id"] <- paste(a1$parent_id, "_t1", sep="")
  a2["target_id"] <- paste(a2$parent_id, "_t2", sep="")
  b1["target_id"] <- a1$target_id
  b2["target_id"] <- a2$target_id
  
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
  results <- list("data"=sl, "anno"=anno, "sim" = list("A_t1"=a1, "A_t2"=a2, "B_t1"=b1, "B_t2"=b2))
  
  return(results)
}


#===============================================================================
#' Combine simulatation details with DTU calls, in preparation for plotting.
#' 
#' @param sim An object generated with rangedata_gen(), containing the simulated data details with which the \code{dtu} object was calculated.
#' @param dtu The object with the DTU results corresponding to the \code{sim} object's data.
#' @return dataframe
#'
#' @export
combine_sim_dtu <- function(sim, dtu) {
  results <- data.table(dtu$Genes[, list(parent_id, pval_AB,  pval_BA, dtu_AB, dtu_BA, dtu)])
  setkey(results, parent_id)
  tmp <- data.table(sim$sim$A_t1[, c("parent_id", "prop", "mag", "fold")])
  setkey(tmp, parent_id)
  results <- merge(results, tmp, by="parent_id")
  #results[, mag := 10^mag]
  #results[, fold := 2^fold]
  return(results)
}


#===============================================================================
#' Plot relationship of parameters.
#' 
#' @param data the output from combine_sim_dtu()
#' @param type 1: dtu, 2: pval_AB, 3: pval_BA, 4: diff dtu_AB - dtu_BA
#' 
#' @export
plot_sim <- function(data, type = 1) {
  if(type ==1){
    ggplot2::ggplot(data, ggplot2::aes(x=fold, y=prop, color=dtu)) +
      ggplot2::labs(x="Ratio A1/A2 fold-change (2^x)", y = "Proportion of A1") +
      ggplot2::scale_color_manual(values=c("blue","red")) +
      ggplot2::geom_point() + 
      ggplot2::facet_grid(. ~ mag)
  } else if(type == 2) {
    ggplot2::ggplot(data, ggplot2::aes(x=fold, y=prop, color=pval_AB)) +
      ggplot2::labs(x="Ratio A1/A2 fold-change (2^x)", y = "Proportion of A1") +
      ggplot2::scale_color_gradientn(colors=c("red", "white", "blue"), values=c(0,0.04999,0.05001,1)) +
      ggplot2::geom_point() + 
      ggplot2::facet_grid(. ~ mag)
  } else if(type == 3) {
    ggplot2::ggplot(data, ggplot2::aes(x=fold, y=prop, color=pval_BA)) +
      ggplot2::labs(x="Ratio A1/A2 fold-change (2^x)", y = "Proportion of A1") +
      ggplot2::scale_color_gradientn(colors=c("red", "white", "blue"), values=c(0,0.04999,0.05001,1)) +
      ggplot2::geom_point() + 
      ggplot2::facet_grid(. ~ mag)
  } else if(type == 4) {
    ggplot2::ggplot(data, ggplot2::aes(x=fold, y=prop, color=as.factor(ifelse(dtu_AB,1,0)+ifelse(dtu_BA,10,0)))) +
      ggplot2::labs(x="Ratio A1/A2 fold-change (2^x)", y = "Proportion of A1") +
      ggplot2::scale_color_manual(values=c("blue", "cyan", "lightblue", "red")) +
      ggplot2::geom_point() + 
      ggplot2::facet_grid(. ~ mag)
  }
}

