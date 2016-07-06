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
#' @return a list with: a minimal sleuth-like object, a corresponding annotation data.frame.
#' 
#' @export
#' 
sim_sleuth_data <- function(varname="condition", COUNTS_COL="est_counts", TARGET_COL="target_id" , 
                            PARENT_COL="parent_id", BS_TARGET_COL="target_id", cnames=c("one","two")) {
  
  tx <- data.frame("target_id"= c("NIB.1", "NIB.2", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.1", "1B1C.2", "CC_a", "CC_b", "1NN", "2NN", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d"), 
                   "parent_id"= c("NIB", "NIB", "1A1N", "1D1C", "1D1C", "1B1C", "1B1C", "CC", "CC", "NN", "NN", "MIX6", "MIX6", "MIX6", "MIX6", "MIX6", "MIX6"))
  names(tx) <- c(TARGET_COL, PARENT_COL)
  
  sl <- list()
  sl[["sample_to_covariates"]] <- data.frame("foo"=c(cnames[1], cnames[2], cnames[1], cnames[2]), 
                                             "bar"=c("ba", "ba", "bb", "bb"))
  names(sl[["sample_to_covariates"]]) <- c(varname, "foodbar")
  
  sl[["kal"]] <- list()
  sl$kal[[1]] <- list()
  sl$kal[[1]]["bootstrap"] <- list()
  sl$kal[[1]]$bootstrap[[1]] <- data.frame("target"= c("NIA1", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.2", "CC_a", "CC_b", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d", "1NN", "2NN"), 
                                           "my_counts"= c(333,  666,    10,       20,       0,          76,         52,       20,     50,     123,       321,       0,         100,       33,        0,        10,    30))
  sl$kal[[1]]$bootstrap[[2]] <- data.frame("target"= c("NIA1", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.2", "CC_a", "CC_b", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d", "1NN", "2NN"), 
                                           "my_counts"= c(330,  660,    11,       21,       0,          80,         55,       22,     52,     125,       320,       0,         104,       30,        0,        11,    29))
  sl$kal[[1]]$bootstrap[[3]] <- data.frame("target"= c("NIA1", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.2", "CC_a", "CC_b", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d", "1NN", "2NN"), 
                                           "my_counts"= c(340,  670,    7,       18,       0,          72,         50,        21,     49,     120,       325,       0,         99,        32,        0,         9,    28))
  
  sl$kal[[2]] <- list()
  sl$kal[[2]]["bootstrap"] <- list()
  sl$kal[[2]]$bootstrap[[1]] <- data.frame("target"= c("NIA1", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.2", "CC_a", "CC_b", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d", "1NN", "2NN"), 
                                           "my_counts"= c(333,  666,    12,       23,       0,          25,         150,      15,     80,     325,       125,       40,         0,        33,        0,        15,    45))
  sl$kal[[2]]$bootstrap[[2]] <- data.frame("target"= c("NIA1", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.2", "CC_a", "CC_b", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d", "1NN", "2NN"), 
                                           "my_counts"= c(323,  656,    15,       22,       0,          23,         160,      20,     90,     320,       135,       50,         0,        30,        0,        20,    60))
  sl$kal[[2]]$bootstrap[[3]] <- data.frame("target"= c("NIA1", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.2", "CC_a", "CC_b", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d", "1NN", "2NN"), 
                                           "my_counts"= c(343,  676,    13,       21,       0,          20,         145,      19,     85,     323,       123,       41,         0,        34,        0,        18,    53))
  
  sl$kal[[3]] <- list()
  sl$kal[[3]]$bootstrap[[1]] <- data.frame("target"= c("NIA1", "1A1N-2", "NIA2", "1D1C:one", "1D1C:two", "1B1C.2", "MIX6.c1", "1A1N-1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.d", "MIX6.nc", "CC_a", "CC_b", "1NN", "2NN"), 
                                           "my_counts"= c(333,  20,       666,    0,          76,         52,       123,       10,       321,       0,         100,       0,        33,        25,     60,     7,     27))
  sl$kal[[3]]$bootstrap[[2]] <- data.frame("target"= c("MIX6.c4", "NIA1", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "1B1C.2", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.nc", "1D1C:two", "MIX6.d", "CC_b", "CC_a", "1NN", "2NN"), 
                                           "my_counts"= c(104,     330,    660,    11,       21,       0,          55,       125,       320,       0,         30,        80,         0,        55,     23,     10,    31))
  sl$kal[[3]]$bootstrap[[3]] <- data.frame("target"= c("NIA1", "NIA2", "1A1N-1", "MIX6.c3", "1A1N-2", "1D1C:one", "MIX6.c4", "1D1C:two", "1B1C.2", "MIX6.c1", "MIX6.c2", "MIX6.nc", "MIX6.d", "CC_a", "CC_b", "1NN", "2NN"), 
                                           "my_counts"= c(340,  670,    7,        0,         18,       0,          99,        72,         50,       120,       325,       32,        0,        22,     50,     9,     32))
  
  sl$kal[[4]] <- list()
  sl$kal[[4]]["bootstrap"] <- list()
  sl$kal[[4]]$bootstrap[[1]] <- data.frame("target"= c("NIA2", "1A1N-1", "NIA1", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.2", "MIX6.c2", "MIX6.c3", "MIX6.c1", "MIX6.c4", "MIX6.nc", "MIX6.d", "CC_b", "CC_a", "1NN", "2NN"), 
                                           "my_counts"= c(666,  12,       333,    23,       0,          25,         150,      125,       40,        325,       0,         33,        0,        93,     22,     22,    61))
  sl$kal[[4]]$bootstrap[[2]] <- data.frame("target"= c("NIA1", "MIX6.d", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "1D1C:two", "MIX6.c1", "1B1C.2", "MIX6.c3", "MIX6.c2", "MIX6.c4", "MIX6.nc", "CC_a", "CC_b", "1NN", "2NN"), 
                                           "my_counts"= c(323,  0,        656,    15,       22,       0,          23,         320,       160,      45,        120,       0,         32,        18,     89,     19,    58))
  sl$kal[[4]]$bootstrap[[3]] <- data.frame("target"= c("NIA1", "NIA2", "1A1N-1", "1A1N-2", "MIX6.nc", "1B1C.2", "MIX6.c1", "MIX6.c2", "MIX6.c3", "1D1C:one", "1D1C:two", "MIX6.c4", "MIX6.d", "CC_a", "CC_b", "1NN", "2NN"), 
                                           "my_counts"= c(343,  676,    13,       21,       34,        145,      323,       123,       41,         0,          20,         0,        0,        23,     80,     17,    50))
  
  names(sl$kal[[1]]$bootstrap[[1]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[1]]$bootstrap[[2]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[1]]$bootstrap[[3]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[2]]$bootstrap[[1]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[2]]$bootstrap[[2]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[2]]$bootstrap[[3]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[3]]$bootstrap[[1]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[3]]$bootstrap[[2]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[3]]$bootstrap[[3]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[4]]$bootstrap[[1]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[4]]$bootstrap[[2]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[4]]$bootstrap[[3]]) <- c(BS_TARGET_COL, COUNTS_COL)
  
  return(list("annot" = tx, "slo" = sl))
}


