#===============================================================================
#' Generate artificial read-count data for simulated 2-transcript genes, 
#' systematically covering a broad range of proportions and magnitudes.
#' 
#' For each level of transcript proportions in condition A, all proportion levels are
#' generated for condition B as possible new states. Therefore the number of transcripts 
#' generated is 2 * N^2 * M, where 2 is the number of gene isoforms per paramater combination, 
#' N is the number of proportion levels and M the number of magnitude levels.
#' 
#' @param propfrom, propto, propby: Range and step of proportions for 1st transcript. Sibling transcript complemented from these.
#' @param magnitudes: vector of read-count sizes to use.
#' @return list containing: \code{data} a minimal sleuth-like object.
#'                          \code{anno} a dataframe matching \code{target_id} to \code{parent_id}.
#'                          \code{sim} a dataframe with the simulation parameters per transcript.
#'
#' @export
countrange_sim <- function(proportions=seq(0, 1, 0.01),
                      magnitudes=c(1,5,10,30,100,500,1000)) {
  # Range of proportions.
  proportions <- 
  
  # Combine proportions and magnitudes.
  # Start by grid-ing the values for easier visualization and naming, instead of computing the outer products directly.
  grd <- expand.grid(c("t1", "t2"), proportions, proportions, magnitudes)
  names(grd) <- c("sibl", "propA", "propB", "mag")
  # Name the "transcripts".
  opts <- options()
  options(scipen=0, digits=3)
  grd["parent_id"] <- paste(grd$mag, "_a", grd$propA,"_b", grd$propB, sep="")
  grd["target_id"] <- paste(grd$parent_id, grd$sibl, sep="_")
  options(scipen=opts$scipen, digits=opts$digits)
  # Complementary proportions
  grd[grd$sibl == "t2", c("propA", "propB")] <- c(1-grd$propA[grd$sibl == "t2"], 1-grd$propB[grd$sibl == "t2"])
  
  # Make a sleuth-like object.
  sl <- list()
  sl["kal"] <- list()
  sl[["sample_to_covariates"]] <- data.frame("sample"=c("A-1", "B-1"), "condition"=c("A", "B"))
  # condition A
  sl$kal[[1]] <- list()
  sl$kal[[1]]["bootstrap"] <- list()
  sl$kal[[1]]$bootstrap[[1]] <- data.frame("target_id"=grd$target_id, "est_counts"=grd$propA * grd$mag)
  # condition B
  sl$kal[[2]] <- list()
  sl$kal[[2]]["bootstrap"] <- list()
  sl$kal[[2]]$bootstrap[[1]] <- data.frame("target_id"=grd$target_id, "est_counts"=grd$propB * grd$mag)
  
  # Create annotation:
  anno <- data.frame("target_id"=grd$target_id, "parent_id"=grd$parent_id)
  
  # Store simulation for ease.
  results <- list("data"=sl, "anno"=anno, "sim"=grd)
  
  return(results)
}


#===============================================================================
#' Combine simulatation details with DTU calls, in preparation for plotting.
#' 
#' @param sim An object generated with countrange_sim(), containing the simulated data details with which the \code{dtu} object was calculated.
#' @param dtu The object with the DTU results corresponding to the \code{sim} object's data.
#' @return dataframe
#'
#' @export
combine_sim_dtu <- function(sim, dtu) {
  results <- data.table::data.table(dtu$Genes[, list(parent_id, pval_AB,  pval_BA, dtu_AB, dtu_BA, dtu, pval_prop_min, dtu_prop)])
  setkey(results, parent_id)
  tmp <- data.table::data.table(sim$sim[sim$sim["sibl"]=="t2", c("parent_id", "propA", "mag", "propB")])
  setkey(tmp, parent_id)
  results <- merge(results, tmp, by="parent_id")
  return(results)
}


#===============================================================================
#' Plot relationship of parameters.
#' 
#' @param data the output from combine_sim_dtu()
#' @param type (str) "AvBvM", "A/BvM", "AvBvMprop", "A/BvMprop"
#' 
#' @export
plot_sim <- function(data, type = "AvBvM") {
  
  if(type =="AvBvM"){
    ggplot2::ggplot(data[order(dtu, propA, propB), ], ggplot2::aes(x=propA, y=propB, color=dtu)) +
      ggplot2::labs(x="Prop t1 in A", y = "Prop t1 in B") +
      ggplot2::scale_color_manual(values=c("lightblue","red")) +
      ggplot2::geom_point() + 
      ggplot2::facet_grid(. ~ mag)
  } else if(type == "A/BvM") {
    ggplot2::ggplot(data[order(dtu, mag), ], ggplot2::aes(x=propB/propA, y=mag, color=dtu)) +
      ggplot2::labs(x="Prop t1 in B / Prop t1 in A", y = "Gene magnitude") +
      ggplot2::scale_color_manual(values=c("lightblue","red")) +
      ggplot2::geom_point()
  } else if(type =="AvBvMprop"){
    ggplot2::ggplot(data[order(dtu, propA, propB), ], ggplot2::aes(x=propA, y=propB, color=dtu_prop)) +
      ggplot2::labs(x="Prop t1 in A", y = "Prop t1 in B") +
      ggplot2::scale_color_manual(values=c("lightblue","red")) +
      ggplot2::geom_point() + 
      ggplot2::facet_grid(. ~ mag)
  } else if(type == "A/BvMprop") {
    ggplot2::ggplot(data[order(dtu, mag), ], ggplot2::aes(x=propB/propA, y=mag, color=dtu_prop)) +
      ggplot2::labs(x="Prop t1 in B / Prop t1 in A", y = "Gene magnitude") +
      ggplot2::scale_color_manual(values=c("lightblue","red")) +
      ggplot2::geom_point()
  }
}

