#==============================================================================
#' Generate an artificial minimal sleuth-like structure, for code-testing or examples.
#' 
#' The default values here should match the default values expected by calculate_DTU().
#' 
#' @param varname Name for the variables by which to compare.
#' @param COUNTS_COL Name for the bootdtraps column containing counts.
#' @param TARGET_COL Name for annotation column containing transcript identification.
#' @param PARENT_COL Name for annotation column containing respective gene identification.
#' @param BS_TARGET_COL Name for bootstraps column containing transcript identification.
#' @param cnames A vector of (two) name values for the comparison variable.
#' @param errannot_inconsistent Logical. Introduces an inconsistency in the transcript IDs, for testing of sanity checks. (FALSE)
#' @param cv_dt Logical. Whether covariates table should be a data.table (FALSE).
#' @return A list with \code{slo} a minimal sleuth-like object, \code{annot} a corresponding annotation data.frame, 
#'         \code{isx} a vector of trancripts common between generated annotation and bootstraps (useful for code testing).
#' 
#' The simulated data will have non-uniform number of bootstraps per sample and non-uniform order of
#' transcript id's among samples and bootstraps. These conditions are unlikely to occur in real data,
#' but this allows testing that the code can handle it.
#' 
#' @import data.table
#' @export
#' 
sim_sleuth_data <- function(varname="condition", COUNTS_COL="est_counts", TARGET_COL="target_id" , 
                            PARENT_COL="parent_id", BS_TARGET_COL="target_id", cnames=c("A","B"), 
                            errannot_inconsistent=FALSE, cv_dt=FALSE)
{
  # !!! Some of the tests of the package are tightly connected to the specifics of the object returned by this function.
  # !!! Entry additions might be tolerated. Changes or removals of entries or structure will certainly cause failures.
  
  tx <- data.frame("target_id"= c("NIB.1", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.1", "1B1C.2", "CC_a", "CC_b", "1NN", "2NN", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d", "LC1", "LC2", "ALLA1", "ALLB1", "ALLB2"), 
                   "parent_id"= c("NIB", "1A1N", "1D1C", "1D1C", "1B1C", "1B1C", "CC", "CC", "NN", "NN", "MIX6", "MIX6", "MIX6", "MIX6", "MIX6", "MIX6", "LC", "LC", "ALLA", "ALLB", "ALLB"))
  names(tx) <- c(TARGET_COL, PARENT_COL)
  
  sl <- list()
  sl[["sample_to_covariates"]] <- NaN 
  if (cv_dt) {
    sl[["sample_to_covariates"]] <- data.table("foo"=c(cnames[1], cnames[2], cnames[1], cnames[2]), 
                                               "bar"=c("ba", "ba", "bb", "bb"))
  } else {
    sl[["sample_to_covariates"]] <- data.frame("foo"=c(cnames[1], cnames[2], cnames[1], cnames[2]), 
                                             "bar"=c("ba", "ba", "bb", "bb"))
  }
  names(sl[["sample_to_covariates"]]) <- c(varname, "batch")
  
  sl[["kal"]] <- list()
  sl$kal[[1]] <- list()
  sl$kal[[1]]["bootstrap"] <- list()
  sl$kal[[1]]$bootstrap[[1]] <- data.frame("target"= c("LC1", "NIA1", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.2", "CC_a", "CC_b", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d", "1NN", "2NN", "LC2", "ALLA1", "ALLB1", "ALLB2"), 
                                           "my_counts"= c(3,   333,    666,    10,       20,       0,          76,         52,       20,     50,     103,       321,       0,         100,       90,        0,        10,    30,    4,     50,      0,       0), stringsAsFactors=FALSE)
  sl$kal[[1]]$bootstrap[[2]] <- data.frame("target"= c("LC1", "NIA1", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.2", "CC_a", "CC_b", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d", "1NN", "2NN", "LC2", "ALLA1", "ALLB1", "ALLB2"), 
                                           "my_counts"= c(2,   310,    680,    11,       21,       0,          80,         55,       22,     52,     165,       320,       0,         130,       80,        0,        11,    29,    5,     40,      0,       0), stringsAsFactors=FALSE)
  sl$kal[[1]]$bootstrap[[3]] <- data.frame("target"= c("LC1", "NIA1", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.2", "CC_a", "CC_b", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d", "1NN", "2NN", "LC2", "ALLA1", "ALLB1", "ALLB2"), 
                                           "my_counts"= c(0,  340,    610,    7,       18,       0,          72,         50,        21,     49,     150,       325,       0,         120,        70,        0,         9,    28,     4,     60,      0,       0), stringsAsFactors=FALSE)
  
  sl$kal[[2]] <- list()
  sl$kal[[2]]["bootstrap"] <- list()
  sl$kal[[2]]$bootstrap[[1]] <- data.frame("target"= c("NIA1", "LC1", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.2", "CC_a", "CC_b", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d", "1NN", "2NN", "LC2", "ALLA1", "ALLB1", "ALLB2"), 
                                           "my_counts"= c(333,  6,     666,    12,       23,       0,          25,         150,      15,     80,     325,       105,       40,         0,        200,        0,        15,    45,    12,     0,      80,       200), stringsAsFactors=FALSE)
  sl$kal[[2]]$bootstrap[[2]] <- data.frame("target"= c("NIA1", "LC1", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.2", "CC_a", "CC_b", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d", "1NN", "2NN", "LC2", "ALLA1", "ALLB1", "ALLB2"), 
                                           "my_counts"= c(323, 1,      606,    15,       22,       0,          23,         190,      20,     90,     270,       115,       30,         0,        150,        0,        20,    60,    15,     0,      120,       250), stringsAsFactors=FALSE)

  sl$kal[[3]] <- list()
  sl$kal[[3]]$bootstrap[[1]] <- data.frame("target"= c("NIA1", "1A1N-2", "LC1", "NIA2", "1D1C:one", "1D1C:two", "1B1C.2", "MIX6.c1", "1A1N-1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.d", "MIX6.nc", "CC_a", "CC_b", "1NN", "2NN", "LC2", "ALLA1", "ALLB1", "ALLB2"), 
                                           "my_counts"= c(333,  20,       3,     666,    0,          76,         52,       100,       10,       360,       0,         100,       0,        180,       25,     60,     7,     27,    13,     35,      0,       0), stringsAsFactors=FALSE)
  sl$kal[[3]]$bootstrap[[2]] <- data.frame("target"= c("MIX6.c4", "NIA1", "LC1", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "1B1C.2", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.nc", "1D1C:two", "MIX6.d", "CC_b", "CC_a", "1NN", "2NN", "LC2", "ALLA1", "ALLB1", "ALLB2"), 
                                           "my_counts"= c(80,      330,    4,     560,    11,       21,       0,          55,       90,        380,       0,         240,        80,         0,        55,     23,     10,    31,    2,     55,      0,       0), stringsAsFactors=FALSE)

  sl$kal[[4]] <- list()
  sl$kal[[4]]["bootstrap"] <- list()
  sl$kal[[4]]$bootstrap[[1]] <- data.frame("target"= c("NIA2", "1A1N-1", "NIA1", "1A1N-2", "1D1C:one", "1D1C:two", "LC1", "1B1C.2", "MIX6.c2", "MIX6.c3", "MIX6.c1", "MIX6.c4", "MIX6.nc", "MIX6.d", "CC_b", "CC_a", "1NN", "2NN", "LC2", "ALLA1", "ALLB1", "ALLB2"), 
                                           "my_counts"= c(666,  12,       333,    23,       0,          25,         0,     150,      155,       40,        300,       0,         33,        0,        93,     22,     22,    61,    11,     0,      110,      210), stringsAsFactors=FALSE)
  sl$kal[[4]]$bootstrap[[2]] <- data.frame("target"= c("NIA1", "MIX6.d", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "LC1", "1D1C:two", "MIX6.c1", "1B1C.2", "MIX6.c3", "MIX6.c2", "MIX6.c4", "MIX6.nc", "CC_a", "CC_b", "1NN", "2NN", "LC2", "ALLA1", "ALLB1", "ALLB2"), 
                                           "my_counts"= c(323,  0,        656,    15,       22,       0,          2,     23,         280,       160,      35,        120,       0,         95,        18,     119,     19,    58,    7,     0,      150,      220), stringsAsFactors=FALSE)
  sl$kal[[4]]$bootstrap[[3]] <- data.frame("target"= c("NIA1", "NIA2", "1A1N-1", "1A1N-2", "MIX6.nc", "1B1C.2", "LC1", "MIX6.c1", "MIX6.c2", "MIX6.c3", "1D1C:one", "1D1C:two", "MIX6.c4", "MIX6.d", "CC_a", "CC_b", "1NN", "2NN", "LC2", "ALLA1", "ALLB1", "ALLB2"), 
                                           "my_counts"= c(343,  676,    13,       21,       55,        145,      3,     270,       133,       30,         0,          20,        0,         0,        23,     80,     17,    50,    14,    0,       130,     200), stringsAsFactors=FALSE)
  
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

#==============================================================================
#' Generate an artificial dataset of bootstrapped abundance estimates, for code-testing or examples.
#' 
#' Based on sim_sleuth_data().
#' 
#' @param cnames A vector of (two) name values for the comparison variable.
#' @return  A list with 3 elements. First a data.frame with the corresponding annotation. Then, 2 lists of data.tables. One list per condition. Each data table represents a sample and contains the estimates from the bootstrap iterations.
#' 
#' @export
#'
sim_boot_data <- function(cnames=c("A", "B")) {
  sim <- sim_sleuth_data(cnames=cnames)
  # Emulate non-sleuth bootstrap data.
  data_A <- denest_sleuth_boots(sim$slo, sim$annot, c(1,3), "est_counts", "target_id")
  data_B <- denest_sleuth_boots(sim$slo, sim$annot, c(2,4), "est_counts", "target_id")
  return(list('annot'= sim$annot, 'boots_A'= data_A, 'boots_B'= data_B))
}

#==============================================================================
#' Generate an artificial dataset of bootstrapped abundance estimates, for code-testing or examples.
#' 
#' Based on sim_sleuth_data().
#' 
#' @param cnames A vector of (two) name values for the comparison variable.
#' @return A list with 3 elements. First, a data.frame with the corresponding annotation table. 
#' Then, 2 data.tables, one per condition. Each data table contains the abundance estimates for all the samples for the respective condition.
#' 
#' @export
#'
sim_count_data <- function(cnames=c("A","B")) {
  sim <- sim_sleuth_data(cnames=cnames)
  # Emulate non-sleuth bootstrap data.
  data_A <- denest_sleuth_boots(sim$slo, sim$annot, c(1,3), "est_counts", "target_id")
  data_B <- denest_sleuth_boots(sim$slo, sim$annot, c(2,4), "est_counts", "target_id")
  # Emulate non-bootstrap data.
  counts_A <- data_A[[1]]
  counts_B <- data_B[[1]]
  return(list(sim$annot, 'counts_A'= counts_A, 'counts_B'= counts_B))
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
#' @param dtuo The object with the DTU results corresponding to the \code{sim} object's data.
#' @return dataframe
#'
#' @import data.table
#' @export
combine_sim_dtu <- function(sim, dtuo) {
  with(dtuo, {
    results <- data.table(Genes[, list(parent_id, pval_AB,  pval_BA, dtu_AB, dtu_BA, dtu, pval_prop_min, dtu_prop)])
    setkey(results, parent_id)
    tmp <- data.table(sim$sim[sim$sim["sibl"]=="t2", c("parent_id", "propA", "mag", "propB")])
    setkey(tmp, parent_id)
    results <- merge(results, tmp, by="parent_id")
    return(results)
  })
}


#===============================================================================
#' Plot relationship of parameters.
#' 
#' @param data the output from combine_sim_dtu()
#' @param type (str) "AvBvM", "B/AvM", "AvBvMprop", "B/AvMprop", "B-AvM"
#' 
#' @import ggplot2
#' @import data.table
#' @export
plot_sim <- function(data, type = "AvBvM") {
  with(data, {
    if(type =="AvBvM"){
      return(
        ggplot(data[order(dtu, propA, propB), ], aes(x=propA, y=propB, color=dtu)) +
          labs(x="Prop t1 in A", y = "Prop t1 in B") +
          scale_color_manual(values=c("lightblue","red")) +
          geom_point() + 
          facet_grid(. ~ mag)
      )
    } else if(type == "B/AvM") {
      return(
        ggplot(data[order(dtu, mag), ], aes(x=propB/propA, y=mag, color=dtu)) +
          labs(x="Prop t1 in B / Prop t1 in A", y = "Gene magnitude") +
          scale_color_manual(values=c("lightblue","red")) +
          geom_point()
        )
    } else if(type =="AvBvMprop"){
      return(
        ggplot(data[order(dtu_prop, propA, propB), ], aes(x=propA, y=propB, color=dtu_prop)) +
          labs(x="Prop t1 in A", y = "Prop t1 in B") +
          scale_color_manual(values=c("lightblue","red")) +
          geom_point() + 
          facet_grid(. ~ mag)
      )
    } else if(type == "B/AvMprop") {
      return(
        ggplot(data[order(dtu_prop, mag), ], aes(x=propB/propA, y=mag, color=dtu_prop)) +
          labs(x="Prop t1 in B / Prop t1 in A", y = "Gene magnitude") +
          scale_color_manual(values=c("lightblue","red")) +
          geom_point()
      )
    }
  })
}
