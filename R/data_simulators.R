#==============================================================================
#' Generate an artificial dataset of bootstrapped abundance estimates, for code-testing or examples.
#' 
#' @param TARGET_COL Name for annotation column containing transcript identification.
#' @param PARENT_COL Name for the bootdtraps column containing counts.
#' @param errannot_inconsistent Logical. Introduces an inconsistency in the transcript IDs, for testing of sanity checks. (FALSE)
#' @param clean Logical. Remove the intentional inconsistency of IDs between annotation and abundances.
#' @return  A list with 3 elements. First a data.frame with the corresponding annotation. 
#' Then, 2 lists of data.tables. One list per condition. Each data table represents a sample and contains the estimates from the bootstrap iterations.
#' 
#' @import data.table
#' 
#' @export
#'
sim_boot_data <- function(errannot_inconsistent=FALSE, PARENT_COL="parent_id", TARGET_COL="target_id", clean=FALSE) {
  tx <- data.frame(target_id= c("1A1B.a", "1A1B.b", "NIB.1", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.1", "1B1C.2", "CC_a", "CC_b", "1NN", "2NN", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d", "LC1", "LC2", "ALLA1", "ALLB1", "ALLB2"), 
                   parent_id= c("1A1B", "1A1B", "NIB", "1A1N", "1D1C", "1D1C", "1B1C", "1B1C", "CC", "CC", "NN", "NN", "MIX6", "MIX6", "MIX6", "MIX6", "MIX6", "MIX6", "LC", "LC", "ALLA", "ALLB", "ALLB"),
                   stringsAsFactors=FALSE)
  names(tx) <- c(TARGET_COL, PARENT_COL)
  
  a <- list()
  b <- list()
  
  if (clean) {
    a[[1]] <- data.table("target"= c("1A1B.a", "1A1B.b", "LC1", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.2", "CC_a", "CC_b", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d", "1NN", "2NN", "LC2", "ALLA1", "ALLB1", "ALLB2", "1B1C.1", "NIB.1"), 
                         "V1"= c(          69,      0,     3,      20,       0,          76,         52,       20,     50,     103,       321,       0,         100,       90,        0,        10,    30,    4,     50,      0,       0,      0,      0), 
                         "V2"= c(          96,      0,     2,      21,       0,          80,         55,       22,     52,     165,       320,       0,         130,       80,        0,        11,    29,    5,     40,      0,       0,      0,      0), 
                         "V3"= c(          88,      0,     0,      18,       0,          72,         50,       21,     49,     150,       325,       0,         120,       70,        0,         9,    28,    4,     60,      0,       0,      0,      0), 
                         stringsAsFactors=FALSE)
    
    a[[2]] <- data.table("target"= c("1A1B.a", "1A1B.b", "1A1N-2", "LC1", "1D1C:one", "1D1C:two", "1B1C.2", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.d", "MIX6.nc", "CC_a", "CC_b", "1NN", "2NN", "LC2", "ALLA1", "ALLB1", "ALLB2", "1B1C.1", "NIB.1"), 
                         "V1"= c(     121,         0,     20,       3,      0,          76,         52,       100,        360,       0,         100,       0,        180,       25,     60,     7,     27,    13,    35,      0,       0,      0,      0), 
                         "V2"= c(     144,         0,     11,       4,      0,          80,         55,        90,        380,       0,         80,        0,        240,       23,     55,     10,    31,    2,     55,      0,       0,      0,      0),
                         stringsAsFactors=FALSE)
    
    b[[1]] <- data.table("target"= c("1A1B.a", "1A1B.b", "LC1",   "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.2", "CC_a", "CC_b", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d", "1NN", "2NN", "LC2", "ALLA1", "ALLB1", "ALLB2", "1B1C.1", "NIB.1"), 
                         "V1"= c(      0,         91,       6,       23,       0,          25,         150,      15,     80,     325,       105,       40,         0,        200,        0,      15,    45,    12,     0,      80,       200,      0,      0), 
                         "V2"= c(      0,        100,       1,       22,       0,          23,         190,      20,     90,     270,       115,       30,         0,        150,        0,      20,    60,    15,     0,      120,      250,      0,      0), 
                         stringsAsFactors=FALSE)
    
    b[[2]] <- data.table("target"= c("1A1N-2", "1D1C:one", "1D1C:two", "LC1", "1B1C.2", "MIX6.c2", "MIX6.c3", "MIX6.c1", "MIX6.c4", "MIX6.nc", "MIX6.d", "CC_b", "CC_a", "1NN", "2NN", "LC2", "ALLA1", "ALLB1", "ALLB2", "1A1B.b", "1A1B.a", "1B1C.1", "NIB.1"), 
                         "V1"= c(         23,       0,          25,         0,     150,      155,       40,        300,       0,         33,        0,        93,     22,     22,    61,    11,     0,      110,     210,      76,        0,      0,      0), 
                         "V2"= c(         22,       0,          23,         2,     160,      120,       35,        280,       0,         95,        0,        119,    18,     19,    58,    7,     0,       150,     220,      69,        0,      0,      0), 
                         "V3"= c(         21,       0,          20,         3,     145,      133,       30,        270,       0,         55,        0,        80,     23,     17,    50,    14,    0,       130,     200,      36,        0,      0,      0), 
                         stringsAsFactors=FALSE)
  } else {
    a[[1]] <- data.table("target"= c("1A1B.a", "1A1B.b", "LC1", "NIA1", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.2", "CC_a", "CC_b", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d", "1NN", "2NN", "LC2", "ALLA1", "ALLB1", "ALLB2"), 
                         "V1"= c(          69,      0,     3,    333,    666,      10,       20,       0,          76,         52,       20,     50,     103,       321,       0,         100,       90,        0,        10,    30,    4,     50,      0,       0), 
                         "V2"= c(          96,      0,     2,    310,    680,      11,       21,       0,          80,         55,       22,     52,     165,       320,       0,         130,       80,        0,        11,    29,    5,     40,      0,       0), 
                         "V3"= c(          88,      0,     0,    340,    610,      7,        18,       0,          72,         50,       21,     49,     150,       325,       0,         120,       70,        0,         9,    28,    4,     60,      0,       0), 
                         stringsAsFactors=FALSE)
    
    a[[2]] <- data.table("target"= c("NIA1", "1A1B.a", "1A1B.b", "1A1N-2", "LC1", "NIA2", "1D1C:one", "1D1C:two", "1B1C.2", "MIX6.c1", "1A1N-1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.d", "MIX6.nc", "CC_a", "CC_b", "1NN", "2NN", "LC2", "ALLA1", "ALLB1", "ALLB2"), 
                         "V1"= c(     333,    121,         0,     20,       3,     666,    0,          76,         52,       100,       10,       360,       0,         100,       0,        180,       25,     60,     7,     27,    13,    35,      0,       0), 
                         "V2"= c(     330,    144,         0,     11,       4,     560,    0,          80,         55,        90,       11,       380,       0,         80,        0,        240,       23,     55,     10,    31,    2,     55,      0,       0),
                         stringsAsFactors=FALSE)
    
    b[[1]] <- data.table("target"= c("1A1B.a", "1A1B.b", "NIA1", "LC1", "NIA2", "1A1N-1", "1A1N-2", "1D1C:one", "1D1C:two", "1B1C.2", "CC_a", "CC_b", "MIX6.c1", "MIX6.c2", "MIX6.c3", "MIX6.c4", "MIX6.nc", "MIX6.d", "1NN", "2NN", "LC2", "ALLA1", "ALLB1", "ALLB2"), 
                         "V1"= c(      0,         91,     333,    6,     666,    12,       23,       0,          25,         150,      15,     80,     325,       105,       40,         0,        200,        0,        15,    45,    12,     0,      80,       200), 
                         "V2"= c(      0,        100,     323,    1,     606,    15,       22,       0,          23,         190,      20,     90,     270,       115,       30,         0,        150,        0,        20,    60,    15,     0,      120,      250), 
                         stringsAsFactors=FALSE)
    
    b[[2]] <- data.table("target"= c("NIA2", "1A1N-1", "NIA1", "1A1N-2", "1D1C:one", "1D1C:two", "LC1", "1B1C.2", "MIX6.c2", "MIX6.c3", "MIX6.c1", "MIX6.c4", "MIX6.nc", "MIX6.d", "CC_b", "CC_a", "1NN", "2NN", "LC2", "ALLA1", "ALLB1", "ALLB2", "1A1B.b", "1A1B.a"), 
                         "V1"= c(     666,    12,       333,    23,       0,          25,         0,     150,      155,       40,        300,       0,         33,        0,        93,     22,     22,    61,    11,     0,      110,     210,      76,        0), 
                         "V2"= c(     656,    15,       323,    22,       0,          23,         2,     160,      120,       35,        280,       0,         95,        0,        119,    18,     19,    58,    7,     0,       150,     220,      69,        0), 
                         "V3"= c(     676,    13,       343,    21,       0,          20,         3,     145,      133,       30,        270,       0,         55,        0,        80,     23,     17,    50,    14,    0,       130,     200,      36,        0), 
                         stringsAsFactors=FALSE)
  }
  
  if (errannot_inconsistent)
    b[[1]][14, 1] <- "Unexpected-name"
  
  return(list('annot'= tx, 'boots_A'= a, 'boots_B'= b))
}

#==============================================================================
#' Generate an artificial dataset of abundance estimates, for code-testing or examples.
#' 
#' @param PARENT_COL Name for the bootdtraps column containing counts.
#' @param TARGET_COL Name for annotation column containing transcript identification.
#' @param errannot_inconsistent Logical. Introduces an inconsistency in the transcript IDs, for testing of sanity checks. (FALSE)
#' @param clean Logical. Remove the intentional inconsistency of IDs between annotation and abundances.
#' @return A list with 3 elements. First, a data.frame with the corresponding annotation table. 
#' Then, 2 data.tables, one per condition. Each data table contains the abundance estimates for all the samples for the respective condition.
#' 
#' @export
#'
sim_count_data <- function(errannot_inconsistent=FALSE, PARENT_COL="parent_id", TARGET_COL="target_id", clean=FALSE) {
  sim <- sim_boot_data(errannot_inconsistent=errannot_inconsistent, PARENT_COL=PARENT_COL, TARGET_COL=TARGET_COL, clean=clean)
  
  return(list('annot'=sim$annot, 'counts_A'= sim$boots_A[[1]], 'counts_B'= sim$boots_B[[1]]))
}
