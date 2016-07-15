#==============================================================================
#' Generate an artificial minimal sleuth-like structure, for code-testing and for examples.
#' 
#' The default values here should match the default values expected by calculate_DTU().
#' 
#' @param varname Name for the variables by which to compare.
#' @param COUNTS_COL Name for the bootdtraps column containing counts.
#' @param TARGET_COL Name for annotation column containing transcript identification.
#' @param PARENT_COL Name for annotation column containing respective gene identification.
#' @param BS_TARGET_COL Name for bootstraps column containing transcript identification.
#' @param cnames A vector of (two) name values for the comparison variable.
#' @param errannot_inconsistent Logical. Introduces an inconsistency in the transcript IDs, for testing of sanity checks.
#' @return a list with \code{slo} a minimal sleuth-like object, \code{annot} a corresponding annotation data.frame, 
#'         \code{isx} a vector of trancripts common between generated annotation and bootstraps (useful for code testing).
#' 
#' The simulated data will have non-uniform number of bootstraps per sample and non-uniform order of
#' transcript id's among samples and bootstraps. These conditions are unlikely to occur in real data,
#' but this allows testing that the code can handle it.
#' 
#' @export
#' 
sim_sleuth_data <- function(varname="condition", COUNTS_COL="est_counts", TARGET_COL="target_id" , 
                            PARENT_COL="parent_id", BS_TARGET_COL="target_id", cnames=c("A","B"), 
                            errannot_inconsistent=FALSE)
{
  tx <- data.frame("target_id"= c("NIB.1", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.1", "1B1C.2", "CC_a", "CC_b", "1NN", "2NN", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d", "LC1", "LC2"), 
                   "parent_id"= c("NIB", "1A1N", "1D1C", "1D1C", "1B1C", "1B1C", "CC", "CC", "NN", "NN", "MIX6", "MIX6", "MIX6", "MIX6", "MIX6", "MIX6", "LC", "LC"))
  names(tx) <- c(TARGET_COL, PARENT_COL)
  
  sl <- list()
  sl[["sample_to_covariates"]] <- data.frame("foo"=c(cnames[1], cnames[2], cnames[1], cnames[2]), 
                                             "bar"=c("ba", "ba", "bb", "bb"))
  names(sl[["sample_to_covariates"]]) <- c(varname, "batch")
  
  sl[["kal"]] <- list()
  sl$kal[[1]] <- list()
  sl$kal[[1]]["bootstrap"] <- list()
  sl$kal[[1]]$bootstrap[[1]] <- data.frame("target"= c("LC1", "NIA1", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.2", "CC_a", "CC_b", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d", "1NN", "2NN", "LC2"), 
                                           "my_counts"= c(3,   333,    666,    10,       20,       0,          76,         52,       20,     50,     123,       321,       0,         100,       33,        0,        10,    30,    4), stringsAsFactors=FALSE)
  sl$kal[[1]]$bootstrap[[2]] <- data.frame("target"= c("LC1", "NIA1", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.2", "CC_a", "CC_b", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d", "1NN", "2NN", "LC2"), 
                                           "my_counts"= c(2,   330,    660,    11,       21,       0,          80,         55,       22,     52,     125,       320,       0,         104,       30,        0,        11,    29,    5), stringsAsFactors=FALSE)
  sl$kal[[1]]$bootstrap[[3]] <- data.frame("target"= c("LC1", "NIA1", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.2", "CC_a", "CC_b", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d", "1NN", "2NN", "LC2"), 
                                           "my_counts"= c(0,  340,    670,    7,       18,       0,          72,         50,        21,     49,     120,       325,       0,         99,        32,        0,         9,    28,     4), stringsAsFactors=FALSE)
  
  sl$kal[[2]] <- list()
  sl$kal[[2]]["bootstrap"] <- list()
  sl$kal[[2]]$bootstrap[[1]] <- data.frame("target"= c("NIA1", "LC1", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.2", "CC_a", "CC_b", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d", "1NN", "2NN", "LC2"), 
                                           "my_counts"= c(333,  6,     666,    12,       23,       0,          25,         150,      15,     80,     325,       125,       40,         0,        33,        0,        15,    45,    7), stringsAsFactors=FALSE)
  sl$kal[[2]]$bootstrap[[2]] <- data.frame("target"= c("NIA1", "LC1", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.2", "CC_a", "CC_b", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d", "1NN", "2NN", "LC2"), 
                                           "my_counts"= c(323, 1,      656,    15,       22,       0,          23,         160,      20,     90,     320,       135,       50,         0,        30,        0,        20,    60,    9), stringsAsFactors=FALSE)

  sl$kal[[3]] <- list()
  sl$kal[[3]]$bootstrap[[1]] <- data.frame("target"= c("NIA1", "1A1N-2", "LC1", "NIA2", "1D1C:one", "1D1C:two", "1B1C.2", "MIX6.c1", "1A1N-1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.d", "MIX6.nc", "CC_a", "CC_b", "1NN", "2NN", "LC2"), 
                                           "my_counts"= c(333,  20,       3,     666,    0,          76,         52,       123,       10,       321,       0,         100,       0,        33,        25,     60,     7,     27,    7), stringsAsFactors=FALSE)
  sl$kal[[3]]$bootstrap[[2]] <- data.frame("target"= c("MIX6.c4", "NIA1", "LC1", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "1B1C.2", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.nc", "1D1C:two", "MIX6.d", "CC_b", "CC_a", "1NN", "2NN", "LC2"), 
                                           "my_counts"= c(104,     330,    4,     660,    11,       21,       0,          55,       125,       320,       0,         30,        80,         0,        55,     23,     10,    31,    2), stringsAsFactors=FALSE)

  sl$kal[[4]] <- list()
  sl$kal[[4]]["bootstrap"] <- list()
  sl$kal[[4]]$bootstrap[[1]] <- data.frame("target"= c("NIA2", "1A1N-1", "NIA1", "1A1N-2", "1D1C:one", "1D1C:two", "LC1", "1B1C.2", "MIX6.c2", "MIX6.c3", "MIX6.c1", "MIX6.c4", "MIX6.nc", "MIX6.d", "CC_b", "CC_a", "1NN", "2NN", "LC2"), 
                                           "my_counts"= c(666,  12,       333,    23,       0,          25,         0,     150,      125,       40,        325,       0,         33,        0,        93,     22,     22,    61,    9), stringsAsFactors=FALSE)
  sl$kal[[4]]$bootstrap[[2]] <- data.frame("target"= c("NIA1", "MIX6.d", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "LC1", "1D1C:two", "MIX6.c1", "1B1C.2", "MIX6.c3", "MIX6.c2", "MIX6.c4", "MIX6.nc", "CC_a", "CC_b", "1NN", "2NN", "LC2"), 
                                           "my_counts"= c(323,  0,        656,    15,       22,       0,          2,     23,         320,       160,      45,        120,       0,         32,        18,     89,     19,    58,    7), stringsAsFactors=FALSE)
  sl$kal[[4]]$bootstrap[[3]] <- data.frame("target"= c("NIA1", "NIA2", "1A1N-1", "1A1N-2", "MIX6.nc", "1B1C.2", "LC1", "MIX6.c1", "MIX6.c2", "MIX6.c3", "1D1C:one", "1D1C:two", "MIX6.c4", "MIX6.d", "CC_a", "CC_b", "1NN", "2NN", "LC2"), 
                                           "my_counts"= c(343,  676,    13,       21,       34,        145,      3,     323,       123,       41,         0,          20,         0,        0,        23,     80,     17,    50,    10), stringsAsFactors=FALSE)
  
  names(sl$kal[[1]]$bootstrap[[1]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[1]]$bootstrap[[2]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[1]]$bootstrap[[3]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[2]]$bootstrap[[1]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[2]]$bootstrap[[2]]) <- c(BS_TARGET_COL, COUNTS_COL)
#   names(sl$kal[[2]]$bootstrap[[3]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[3]]$bootstrap[[1]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[3]]$bootstrap[[2]]) <- c(BS_TARGET_COL, COUNTS_COL)
#   names(sl$kal[[3]]$bootstrap[[3]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[4]]$bootstrap[[1]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[4]]$bootstrap[[2]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[4]]$bootstrap[[3]]) <- c(BS_TARGET_COL, COUNTS_COL)
  
  if (errannot_inconsistent) {
    sl$kal[[3]]$bootstrap[[1]][10, BS_TARGET_COL] <- "Unexpected-name"
  }
  
  return(list("annot" = tx, "slo" = sl, "isx" = intersect(tx[[TARGET_COL]], sl$kal[[1]]$bootstrap[[1]][[BS_TARGET_COL]]) ))
}


#===============================================================================
#' Generate artificial read-count data for simulated 2-transcript genes, 
#' systematically covering a broad range of proportions and magnitudes.
#' 
#' For each level of transcript proportions in condition A, all proportion levels are
#' generated for condition B as possible new states. Therefore the number of transcripts 
#' generated is 2 * N^2 * M, where 2 is the number of gene isoforms per paramater combination, 
#' N is the number of proportion levels and M the number of magnitude levels.
#' 
#' @param proportions Range and step of proportions for 1st transcript. Sibling transcript complemented from these.
#' @param magnitudes vector of read-count sizes to use.
#' @return list containing: \code{data} a minimal sleuth-like object.
#'                          \code{anno} a dataframe matching \code{target_id} to \code{parent_id}.
#'                          \code{sim} a dataframe with the simulation parameters per transcript.
#'
#' @export
countrange_sim <- function(proportions= seq(0, 1, 0.01),
                           magnitudes= c(1,5,10,30,100,500,1000)) {
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
#' @param type (str) "AvBvM", "B/AvM", "AvBvMprop", "B/AvMprop", "B-AvM"
#' 
#' @export
plot_sim <- function(data, type = "AvBvM") {
  
  if(type =="AvBvM"){
    ggplot2::ggplot(data[order(dtu, propA, propB), ], ggplot2::aes(x=propA, y=propB, color=dtu)) +
      ggplot2::labs(x="Prop t1 in A", y = "Prop t1 in B") +
      ggplot2::scale_color_manual(values=c("lightblue","red")) +
      ggplot2::geom_point() + 
      ggplot2::facet_grid(. ~ mag)
  } else if(type == "B/AvM") {
    ggplot2::ggplot(data[order(dtu, mag), ], ggplot2::aes(x=propB/propA, y=mag, color=dtu)) +
      ggplot2::labs(x="Prop t1 in B / Prop t1 in A", y = "Gene magnitude") +
      ggplot2::scale_color_manual(values=c("lightblue","red")) +
      ggplot2::geom_point()
  } else if(type =="AvBvMprop"){
    ggplot2::ggplot(data[order(dtu_prop, propA, propB), ], ggplot2::aes(x=propA, y=propB, color=dtu_prop)) +
      ggplot2::labs(x="Prop t1 in A", y = "Prop t1 in B") +
      ggplot2::scale_color_manual(values=c("lightblue","red")) +
      ggplot2::geom_point() + 
      ggplot2::facet_grid(. ~ mag)
  } else if(type == "B/AvMprop") {
    ggplot2::ggplot(data[order(dtu_prop, mag), ], ggplot2::aes(x=propB/propA, y=mag, color=dtu_prop)) +
      ggplot2::labs(x="Prop t1 in B / Prop t1 in A", y = "Gene magnitude") +
      ggplot2::scale_color_manual(values=c("lightblue","red")) +
      ggplot2::geom_point()
  }
}
