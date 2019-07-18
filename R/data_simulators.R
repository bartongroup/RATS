########## ########## ########## ########## ########## ########## ########## ########## ##########
# Helper functions for testing and diagnostics. 
########## ########## ########## ########## ########## ########## ########## ########## ##########


#==============================================================================
#' Generate an artificial dataset of bootstrapped abundance estimates, for code-testing or examples.
#' 
#' Try to cover a range of normal and outlier scenarios that might trip up the filters.
#' 
#' @param TARGET_COL Name for annotation column containing transcript identification.
#' @param PARENT_COL Name for the bootdtraps column containing counts.
#' @param clean Logical. Don;t include annotation inconsistencies.
#' @return  A list with 3 elements. First a data.frame with the corresponding annotation. Then, 2 lists of data.tables. One list per condition, each data.table represents a sample and contains the estimates from the bootstrap iterations.
#' 
#' @import data.table
#' 
#' @export
#'
sim_boot_data <- function(PARENT_COL="parent_id", TARGET_COL="target_id", clean=TRUE) {

  a <- list()
  b <- list()
  
  # Start with clean perfectly consistent annotation and quantification data. Still out of order though, to ensure that tidying up works correctly.
  #
  # SOLO    gene with 1 transcript, expressed in both A and B, not differentially expressed.
  # SAME    gene with 2 transcripts, expressed in both A and B, neither differentially expressed.
  # LONE    gene with 1 transcript, expressed in both A and B, differentially expressed.
  # D2TE    gene with 2 transcripts, expressed in both A and B, differentially expressed, but not differentially used.
  # D2TU    gene with 2 transcripts, expressed in both A and B, a little differentially used, no switching
  # SW      gene with 2 transcripts, expressed in both A and B, differentially used, with switching
  # XSW     gene with 2 transcripts, one expressed in A, one in B, (extreme differential usage).
  # FAKE    gene with 2 transcripts, one expressed in A and B, one not expressed in either.
  # ALLA    gene with 2 transcripts, both expressed only in A.
  # LC      gene with 2 transcripts of low expression, differentially expressed and used.
  # NN      gene with 2 transcripts, neither expressed.
  # MIX     gene with multiple transcripts, mixing the above scenarios all in one.
  # NID     gene annotated but not present in quantifications
  # NIA     gene quantified but not present in the annotation
  
  tx <- data.frame(target_id= c("SW2", "D2TE_b", "SAME_1", "SAME_2", "LONE.a", "D2TE_a", "SOLO.1", "2D2TU", "XSW:one", "XSW:two", "FAKE-1", "LC1", "FAKE-2", "ALLA:1", "ALLA:2", "LC2", "NNa", "NNb", "MIX.n", "MIX.ab", "1D2TU", "MIX.l2", "MIX.l1", "MIX.a", "MIX.b", "SW1"), 
                   parent_id= c("SW", "D2TE", "SAME", "SAME", "LONE", "D2TE", "SOLO", "D2TU", "XSW", "XSW", "FAKE", "LC", "FAKE", "ALLA", "ALLA", "LC", "NN", "NN", "MIX", "MIX", "D2TU", "MIX", "MIX", "MIX", "MIX", "SW"), 
                   stringsAsFactors=FALSE)
  
  # Keep "bootstraps" identical in each sample, to make it easy to mentally predict the outcome and create equivalent unbootstrapped data.
  
  a[[1]] <- data.table("target"= c("SOLO.1", "SAME_1", "SAME_2", "LONE.a", "D2TE_a", "D2TE_b", "1D2TU", "2D2TU", "SW1", "SW2", "XSW:one", "XSW:two", "FAKE-1", "FAKE-2", "ALLA:1", "ALLA:2", "LC1", "LC2", "NNa", "NNb", "MIX.n", "MIX.ab", "MIX.l1", "MIX.l2", "MIX.a", "MIX.b"),
                       "V1"= c(     100,      200,      300,      150,      300,      400,      250,     250,     20,    80,    90,        0,         0,        150,      100,      40,       10,    20,    0,     0,     0,       50,       10,       20,       100,     50    ), 
                       "V2"= c(     100,      200,      300,      150,      300,      400,      250,     250,     20,    80,    90,        0,         0,        150,      100,      40,       10,    20,    0,     0,     0,       50,       10,       20,       100,     50    ), 
                       "V3"= c(     100,      200,      300,      150,      300,      400,      250,     250,     20,    80,    90,        0,         0,        150,      100,      40,       10,    20,    0,     0,     0,       50,       10,       20,       100,     50    ), 
                       stringsAsFactors=FALSE)
  
  a[[2]] <- data.table("target"= c("SOLO.1", "SAME_1", "SAME_2", "LONE.a", "D2TE_a", "D2TE_b", "1D2TU", "2D2TU",  "SW1", "SW2", "XSW:one", "XSW:two", "FAKE-1", "FAKE-2", "ALLA:1", "ALLA:2", "LC1", "LC2", "NNa", "NNb", "MIX.n", "MIX.ab", "MIX.l1", "MIX.l2", "MIX.a", "MIX.b"),
                       "V1"= c(     120,      240,      250,      120,      400,      500,      270,     300,      30,    100,   85,        5,         5,        150,      115,      30,       5,     15,    0,     0,     0,       40,       0,        30,       80,      40    ), 
                       "V2"= c(     120,      240,      250,      120,      400,      500,      270,     300,      30,    100,   85,        5,         5,        150,      115,      30,       5,     15,    0,     0,     0,       40,       0,        30,       80,      40    ), 
                       stringsAsFactors=FALSE)
  # Shuffle order
  a[[2]] <- a[[2]][sample(seq.int(1, nrow(a[[2]]))) ,]
  
  
  b[[1]] <- data.table("target"= c("SOLO.1", "SAME_1", "SAME_2", "LONE.a", "D2TE_a", "D2TE_b", "1D2TU", "2D2TU",  "SW1", "SW2", "XSW:one", "XSW:two", "FAKE-1", "FAKE-2", "ALLA:1", "ALLA:2", "LC1", "LC2", "NNa", "NNb", "MIX.n", "MIX.ab", "MIX.l1", "MIX.l2", "MIX.a", "MIX.b"),
                       "V1"= c(     100,      200,      300,      300,      150,      200,      250,     400,      220,   80,    9,         50,        0,        50,       10,      0,        14,    10,    0,     0,     0,       70,       10,       2,        50,      800   ), 
                       "V2"= c(     100,      200,      300,      300,      150,      200,      250,     400,      220,   80,    9,         50,        0,        50,       10,      0,        14,    10,    0,     0,     0,       70,       10,       2,        50,      800   ), 
                       stringsAsFactors=FALSE)
  
  b[[2]] <- data.table("target"= c("SOLO.1", "SAME_1", "SAME_2", "LONE.a", "D2TE_a", "D2TE_b", "1D2TU", "2D2TU",  "SW1", "SW2", "XSW:one", "XSW:two", "FAKE-1", "FAKE-2", "ALLA:1", "ALLA:2", "LC1", "LC2", "NNa", "NNb", "MIX.n", "MIX.ab", "MIX.l1", "MIX.l2", "MIX.a", "MIX.b"),
                       "V1"= c(     100,      200,      300,      300,      100,      200,      270,     500,      150,   100,   9,         70,        0,        50,       10,       0,        25,    10,    0,     0,     0,       50,       20,       10,       30,      600   ), 
                       "V2"= c(     100,      200,      300,      300,      100,      200,      270,     500,      150,   100,   9,         70,        0,        50,       10,       0,        25,    10,    0,     0,     0,       50,       20,       10,       30,      600   ), 
                       "V3"= c(     100,      200,      300,      300,      100,      200,      270,     500,      150,   100,   9,         70,        0,        50,       10,       0,        25,    10,    0,     0,     0,       50,       20,       10,       30,      600   ), 
                       stringsAsFactors=FALSE)

  
  if (!clean) {
    # An unquantified annotation entry.
    tx <- rbind(tx, data.frame(target_id=c("NID1", "NID2"), parent_id=c("NID", "NID"), stringsAsFactors = FALSE))
    # An unnanotated quantification.
    a[[1]] <- rbind(a[[1]], data.frame("target"=c("NIA1", "NIA2"), "V1"=c(10,10), "V2"=c(20,20), "V3"=c(30,30), stringsAsFactors = FALSE))
    a[[2]] <- rbind(a[[2]], data.frame("target"=c("NIA1", "NIA2"), "V1"=c(10,10), "V2"=c(20,20), stringsAsFactors = FALSE))
    b[[1]] <- rbind(b[[1]], data.frame("target"=c("NIA1", "NIA2"), "V1"=c(10,10), "V2"=c(20,20), stringsAsFactors = FALSE))
    b[[2]] <- rbind(b[[2]], data.frame("target"=c("NIA1", "NIA2"), "V1"=c(10,10), "V2"=c(20,20), "V3"=c(30,30), stringsAsFactors = FALSE))
  }  
  
  names(tx) <- c(TARGET_COL, PARENT_COL)

  return(list('annot'= tx, 'boots_A'= a, 'boots_B'= b))
}

#==============================================================================
#' Generate an artificial dataset of abundance estimates, for code-testing or examples.
#' 
#' @param PARENT_COL Name for the bootdtraps column containing counts.
#' @param TARGET_COL Name for annotation column containing transcript identification.
#' @param clean Logical. Remove the intentional inconsistency of IDs between annotation and abundances.
#' @return A list with 3 elements. First, a data.frame with the corresponding annotation table. 
#' Then, 2 data.tables, one per condition. Each data table contains the abundance estimates for all the samples for the respective condition.
#' 
#' @export
#'
sim_count_data <- function(errannot_inconsistent=FALSE, PARENT_COL="parent_id", TARGET_COL="target_id", clean=TRUE) {
  tx <- data.frame(target_id= c("SW2", "D2TE_b", "SAME_1", "SAME_2", "LONE.a", "D2TE_a", "SOLO.1", "2D2TU", "XSW:one", "XSW:two", "FAKE-1", "LC1", "FAKE-2", "ALLA:1", "ALLA:2", "LC2", "NNa", "NNb", "MIX.n", "MIX.ab", "1D2TU", "MIX.l2", "MIX.l1", "MIX.a", "MIX.b", "SW1"), 
                   parent_id= c("SW", "D2TE", "SAME", "SAME", "LONE", "D2TE", "SOLO", "D2TU", "XSW", "XSW", "FAKE", "LC", "FAKE", "ALLA", "ALLA", "LC", "NN", "NN", "MIX", "MIX", "D2TU", "MIX", "MIX", "MIX", "MIX", "SW"), 
                   stringsAsFactors=FALSE)
  
  # Start with clean perfectly consistent annotation and quantification data. Still out of order though, to ensure that tidying up works correctly.
  #
  # SOLO    gene with 1 transcript, expressed in both A and B, not differentially expressed.
  # SAME    gene with 2 transcripts, expressed in both A and B, neither differentially expressed.
  # LONE    gene with 1 transcript, expressed in both A and B, differentially expressed.
  # D2TE    gene with 2 transcripts, expressed in both A and B, differentially expressed, but not differentially used.
  # D2TU    gene with 2 transcripts, expressed in both A and B, a little differentially used, no switching
  # SW      gene with 2 transcripts, expressed in both A and B, differentially used, with switching
  # XSW     gene with 2 transcripts, one expressed in A, one in B, (extreme differential usage).
  # FAKE    gene with 2 transcripts, one expressed in A and B, one not expressed in either.
  # ALLA    gene with 2 transcripts, both expressed only in A.
  # LC      gene with 2 transcripts of low expression, differentially expressed and used.
  # NN      gene with 2 transcripts, neither expressed.
  # MIX     gene with multiple transcripts, mixing the above scenarios all in one.
  # NID     gene annotated but not present in quantifications
  # NIA     gene quantified but not present in the annotation
  
  a <- data.table("target"= c("SOLO.1", "SAME_1", "SAME_2", "LONE.a", "D2TE_a", "D2TE_b", "1D2TU", "2D2TU", "SW1", "SW2", "XSW:one", "XSW:two", "FAKE-1", "FAKE-2", "ALLA:1", "ALLA:2", "LC1", "LC2", "NNa", "NNb", "MIX.n", "MIX.ab", "MIX.l1", "MIX.l2", "MIX.a", "MIX.b"),
                  "V1"= c(     100,      200,      300,      150,      300,      400,      250,     250,     20,    80,    90,        0,         0,        150,      100,      40,       10,    20,    0,     0,     0,       50,       10,       20,       100,     50    ), 
                  "V2"= c(     120,      240,      250,      120,      400,      500,      270,     300,     30,    100,   85,        5,         5,        150,      115,      30,       5,     15,    0,     0,     0,       40,       0,        30,       80,      40    ), 
                  stringsAsFactors=FALSE)
  
  b <- data.table("target"= c("SOLO.1", "SAME_1", "SAME_2", "LONE.a", "D2TE_a", "D2TE_b", "1D2TU", "2D2TU", "SW1", "SW2", "XSW:one", "XSW:two", "FAKE-1", "FAKE-2", "ALLA:1", "ALLA:2", "LC1", "LC2", "NNa", "NNb", "MIX.n", "MIX.ab", "MIX.l1", "MIX.l2", "MIX.a", "MIX.b"),
                  "V1"= c(     100,      200,      300,      300,      150,      200,      250,     400,      220,   80,    9,         50,        0,        50,       10,      0,        14,    10,    0,     0,     0,       70,       10,       2,        50,      800   ), 
                  "V2"= c(     100,      200,      300,      300,      100,      200,      270,     500,      150,   100,   9,         70,        0,        50,       10,      0,        25,    10,    0,     0,     0,       50,       20,       10,       30,      600   ), 
                  stringsAsFactors=FALSE)
  # Shuffle order
  b <- b[sample(seq.int(1, nrow(b))) ,]
  
  if (!clean) {
    # An unquantified annotation entry.
    tx <- rbind(tx, data.frame(target_id=c("NID1", "NID2"), parent_id=c("NID", "NID"), stringsAsFactors = FALSE))
    # An unnanotated quantification.
    a <- rbind(a, data.frame("target"=c("NIA1", "NIA2"), "V1"=c(10,10), "V2"=c(20,20), stringsAsFactors = FALSE))
    b <- rbind(b, data.frame("target"=c("NIA1", "NIA2"), "V1"=c(10,10), "V2"=c(20,20), stringsAsFactors = FALSE))
  }  
  
  names(tx) <- c(TARGET_COL, PARENT_COL)
  
  return(list('annot'= tx, 'counts_A'= a, 'counts_B'= b))
}
