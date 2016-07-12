#================================================================================
#' Calculate differential transcript usage.
#'
#' @param slo A sleuth object.
#' @param annot A dataframe matching the transcript identifiers to their corresponding gene identifiers.
#' @param name_A The name for one condition, as it appears in the \code{sample_to_covariates} table within the sleuth object.
#' @param name_B The name for the other condition, as it appears in the \code{sample_to_covariates} table within the sleuth object.
#' @param varname The name of the covariate to which the two conditions belong, as it appears in the \code{sample_to_covariates} table within the sleuth object. Default \code{"condition"}.
#' @param verbose Display progress updates, default \code{FALSE}.
#' @param p_thresh The p-value threshold, default 0.05.
#' @param count_thresh Minimum count of fragments per sample, in at least one of the conditions, for transcripts to be considered (default 5).
#' @param dprop_thresh Minimum change in proportion of a transcript for it to be considered (default 0.1).
#' @param testmode One of "G-test", "prop-test", "both" (default both).
#' @param correction The p-value correction to apply, as defined in \code{stats::p.adjust.methods}, default \code{"BH"}.
#' @param boots Bootstrap the p-values of either test. One of "G-test", "prop-test", "both", "none". (default "none").
#' @param bootnum Number of bootstraps (default 10000).
#' @param threads (Not in use currently) Enable parallel processing. Default uses parallel::detectCores().
#' @param COUNTS_COL The name of the counts column to use for the DTU calculation (est_counts or tpm), default \code{"est_counts"}.
#' @param TARGET_COL The name of the transcript identifier column in the annot object, default \code{"target_id"}
#' @param PARENT_COL The name of the parent identifier column in the annot object, default \code{"parent_id"}.
#' @param BS_TARGET_COL The name of the transcript identifier column in the sleuth bootstrap tables, default \code{"target_id"}.
#' @return List of data tables, with gene-level and transcript-level information.
#'
#' @import data.table
#' @import matrixStats
#' @export
call_DTU <- function(slo, annot, name_A, name_B, varname= "condition", 
                          p_thresh= 0.05, count_thresh= 5, dprop_thresh= 0.1, testmode= "both", correction= "BH", 
                          verbose= FALSE, boots= "none", bootnum= 100L, threads= parallel::detectCores(),
                          COUNTS_COL= "est_counts", TARGET_COL= "target_id", PARENT_COL= "parent_id", BS_TARGET_COL= "target_id") {
  #---------- PREP
  
  # Input checks.
  paramcheck <- parameters_good(slo, annot, name_A, name_B, varname, COUNTS_COL,
                                correction, p_thresh, TARGET_COL, PARENT_COL, BS_TARGET_COL, verbose, threads, count_thresh, testmode, 
                                boots, bootnum, dprop_thresh)
  if (paramcheck$error) stop(paramcheck$message)
  threads <- as.integer(threads)  # Plain numbers default to double, unless integer R syntax is explicitly used.
  bootnum <- as.integer(bootnum)
  
  # Set up progress bar
  progress <- init_progress(verbose)
  progress <- update_progress(progress)
  
  #----------- LOOK-UP
  
  progress <- update_progress(progress)
  # Look-up from target_id to parent_id.
  tx_filter <- data.table(target_id = annot[[TARGET_COL]], parent_id = annot[[PARENT_COL]])
  setkey(tx_filter, parent_id, target_id) # Fixates the order of genes and transcripts to be used throughout the rest of this package.
  # Reverse look-up from replicates to covariates.
  samples_by_condition <- group_samples(slo$sample_to_covariates)[[varname]]
  
  #---------- EXTRACT DATA
  
  progress <- update_progress(progress)
  # De-nest, average and index the counts from the bootstraps. 
  # For each condition, a list holds a dataframe for each replicate. 
  # Each dataframe holds the counts from ALL the bootstraps. Target_id is included but NOT used as key so as to ensure the order keeps matching tx_filter.
  data_A <- denest_boots(slo, tx_filter$target_id, samples_by_condition[[name_A]], COUNTS_COL, BS_TARGET_COL )
  data_B <- denest_boots(slo, tx_filter$target_id, samples_by_condition[[name_B]], COUNTS_COL, BS_TARGET_COL )
  bootmeans_A <- as.data.table(lapply(data_A, function(b) { n <- names(b); rowMeans(b[, n[1:length(n)-1], with=FALSE]) }))  # Tables don't have access to column ranges by index so I have todevtools fis out the names?
  bootmeans_B <- as.data.table(lapply(data_B, function(b) { n <- names(b); rowMeans(b[, n[1:length(n)-1], with=FALSE]) }))
  
  #---------- TEST
  
  progress <- update_progress(progress)
  # Do the core work.
  suppressWarnings(resobj <- calculate_DTU(bootmeans_A, bootmeans_B, tx_filter, testmode, "full", count_thresh, dprop_thresh, correction))
  
  #-------- ADD INFO
  
  progress <- update_progress(progress)
  # Fill in run info.
  resobj$Parameters["var_name"] <- varname
  resobj$Parameters["cond_A"] <- name_A
  resobj$Parameters["cond_B"] <- name_B
  resobj$Parameters["p_thresh"] <- p_thresh 
  resobj$Parameters["count_thresh"] <- count_thresh
  resobj$Parameters["dprop_thresh"] <- dprop_thresh
  resobj$Parameters["tests"] <- testmode
  resobj$Parameters["bootstrap"] <- boots
  resobj$Parameters["bootnum"] <- bootnum
  resobj$Parameters["threads"] <- threads
  # Fill in data info.
  resobj$Genes[, known_transc :=  resobj$Transcripts[, length(target_id), by=parent_id][, V1] ]  # V1 is the automatic column name for the lengths in the subsetted data.table
  detected = resobj$Transcripts[, abs((propA + propB))] > 0
  resobj$Genes[, detect_transc :=  resobj$Transcripts[, .(parent_id, ifelse(abs(propA + propB) > 0, 1, 0))][, sum(V2), by = parent_id][, V1] ]
  resobj$Genes[(is.na(detect_transc)), detect_transc := 0]
  resobj$Transcripts[, meanA :=  rowMeans(bootmeans_A) ]
  resobj$Transcripts[, meanB :=  rowMeans(bootmeans_B) ]
  resobj$Transcripts[, stdevA :=  rowSds(as.matrix(bootmeans_A)) ]
  resobj$Transcripts[, stdevB :=  rowSds(as.matrix(bootmeans_B)) ]
  # Process results.
  if (any(testmode == c("prop-test", "both"))) {
    resobj$Transcripts[, Pt_DTU := Pt_pval_corr < p_thresh]
    resobj$Genes[, Pt_DTU :=  resobj$Transcripts[, any(Pt_DTU), by = parent_id][, V1] ]
  }
  if (any(testmode == c("G-test", "g-test", "both"))) {
    resobj$Genes[, Gt_dtuAB := Gt_pvalAB_corr < p_thresh]
    resobj$Genes[, Gt_dtuBA := Gt_pvalBA_corr < p_thresh]
    resobj$Genes[, Gt_DTU := Gt_dtuAB & Gt_dtuBA ]
    resobj$Transcripts[, Gt_DTU := merge(resobj$Genes[, .(parent_id, Gt_DTU)], resobj$Transcripts[, .(parent_id)])[, Gt_DTU] ]
  }
  
  #---------- BOOTSTRAP
  
  progress <- update_progress(progress)
  if (boots != "none") {
    
    #----- Iterations
    
    bootres <- lapply(1:bootnum, function(b) {
      # Grab a bootstrap from each replicate. 
      counts_A <- as.data.table(lapply(data_A, function(smpl) { smpl[[sample( names(smpl)[1:(dim(smpl)[2]-1)], 1)]] }))  # Have to use list syntax to get a vector back. 
                                                                                                                         # Usual table syntax with "with=FALSE" returns a table and I fail to cast it.
                                                                                                                         # Also, the last column is the target_id, so I leave it out.
      counts_B <- as.data.table(lapply(data_B, function(smpl) { smpl[[sample( names(smpl)[1:(dim(smpl)[2]-1)], 1)]] }))
      # Do the work.
      # Ignore warning. Chi-square test generates warnings for counts <5. This is expected behaviour. Transcripts changing between off and on are often culprits.
      suppressWarnings(bout <- calculate_DTU(counts_A, counts_B, tx_filter, testmode, "short", count_thresh, dprop_thresh, correction))
      return(list("pp" = bout$Transcripts[, Pt_pval_corr], 
                  "pgab" = bout$Genes[, Gt_pvalAB_corr],
                  "pgba" = bout$Genes[, Gt_pvalBA_corr] ))
    })
    
    #----- Stats
    
    if (any(boots == c("prop-test", "both"))) {
      # !!! POSSIBLE source of ERRORS if bootstraps * transcripts exceed R's maximum matrix size. (due to number of either) !!!
      pres <- as.matrix(as.data.table(lapply(bootres, function(b) { b[["pp"]] })))
      # Only for eligible transcripts.
      resobj$Transcripts[(eligible), Pt_boot_dtu := apply(pres[resobj$Transcripts[, eligible], ], MARGIN=1, function(row) { mean(ifelse(row > p_thresh | is.na(row), 0, 1)) } )]  # fraction of passes
      resobj$Transcripts[(eligible), Pt_boot_mean := rowMeans(pres[resobj$Transcripts[, eligible], ], na.rm = TRUE)]
      resobj$Transcripts[(eligible), Pt_boot_stdev := rowSds(pres[resobj$Transcripts[, eligible], ], na.rm = TRUE)]
      resobj$Transcripts[(eligible), Pt_boot_min := rowMins(pres[resobj$Transcripts[, eligible], ], na.rm = TRUE)]
      resobj$Transcripts[(eligible), Pt_boot_max := rowMaxs(pres[resobj$Transcripts[, eligible], ], na.rm = TRUE)]
      resobj$Transcripts[(eligible), Pt_boot_na := rowCounts(pres[resobj$Transcripts[, eligible], ], value = NA)]
    }
    if (any(boots == c("G-test", "g-test", "both"))) {
      # !!! POSSIBLE source of ERRORS if bootstraps * genes exceed R's maximum matrix size. (due to number of bootstraps) !!!
      gabres <- as.matrix(as.data.table(lapply(bootres, function(b) { b[["pgab"]] })))
      gbares <- as.matrix(as.data.table(lapply(bootres, function(b) { b[["pgba"]] })))
      resobj$Genes[(eligible), Gt_boot_dtuAB := apply(gabres[resobj$Genes[, eligible], ], MARGIN=1, function(row) { mean(ifelse(row > p_thresh | is.na(row), 0, 1)) } )]
      resobj$Genes[(eligible), Gt_boot_dtuBA := apply(gbares[resobj$Genes[, eligible], ], MARGIN=1, function(row) { mean(ifelse(row > p_thresh | is.na(row), 0, 1)) } )]
      resobj$Genes[(eligible), Gt_boot_meanAB := rowMeans(gabres[resobj$Genes[, eligible], ], na.rm = TRUE)]
      resobj$Genes[(eligible), Gt_boot_meanBA := rowMeans(gbares[resobj$Genes[, eligible], ], na.rm = TRUE)]
      resobj$Genes[(eligible), Gt_boot_stdevAB := rowSds(gabres[resobj$Genes[, eligible], ], na.rm = TRUE)]
      resobj$Genes[(eligible), Gt_boot_stdevBA := rowSds(gbares[resobj$Genes[, eligible], ], na.rm = TRUE)]
      resobj$Genes[(eligible), Gt_boot_minAB := rowMins(gabres[resobj$Genes[, eligible], ], na.rm = TRUE)]
      resobj$Genes[(eligible), Gt_boot_minBA := rowMins(gbares[resobj$Genes[, eligible], ], na.rm = TRUE)]
      resobj$Genes[(eligible), Gt_boot_maxAB := rowMaxs(gabres[resobj$Genes[, eligible], ], na.rm = TRUE)]
      resobj$Genes[(eligible), Gt_boot_maxBA := rowMaxs(gbares[resobj$Genes[, eligible], ], na.rm = TRUE)]
      resobj$Genes[(eligible), Gt_boot_na := rowCounts(gabres[resobj$Genes[, eligible], ], value = NA)]  # It doesn't matter if AB or BA, affected identically by gene eligibility.
    }
  }
  
  #---------- DONE
  
  progress <- update_progress(progress)

  return(resobj)
}









#================================================================================
#' Check input parameters.
#'
#' @return List with a logical value and a message.
#'
parameters_good <- function(slo, annot, name_A, name_B, varname, COUNTS_COL,
                            correction, p_thresh, TARGET_COL, PARENT_COL, BS_TARGET_COL, verbose, 
                            threads, count_thresh, testmode, boots, bootnum, dprop_thresh) 
{
  if (!is.data.frame(annot))
    return(list("error"=TRUE, "message"="The provided annot is not a data.frame!"))
  
  if (any(! c(TARGET_COL, PARENT_COL) %in% names(annot)))
    return(list("error"=TRUE, "message"="The specified target and/or parent IDs field-names do not exist in annot!"))
  if (! BS_TARGET_COL %in% names(slo$kal[[1]]$bootstrap[[1]]))
    return(list("error"=TRUE, "message"="The specified target IDs field-name does not exist in the bootstraps!"))
  
  if (! COUNTS_COL %in% names(slo$kal[[1]]$bootstrap[[1]]))
    return(list("error"=TRUE, "message"="The specified counts field-name does not exist!"))
  
  if (! correction %in% p.adjust.methods)
    return(list("error"=TRUE, "message"="Invalid p-value correction method name. Refer to stats::p.adjust.methods!"))
  
  if ((!is.numeric(p_thresh)) || p_thresh > 1 || p_thresh < 0)
    return(list("error"=TRUE, "message"="Invalid p-value threshold!"))
  
  if (! varname %in% names(slo$sample_to_covariates))
    return(list("error"=TRUE, "message"="The specified covariate name does not exist!"))
  
  if (any(! c(name_A, name_B) %in% slo$sample_to_covariates[[varname]] ))
    return(list("error"=TRUE, "message"="One or both of the specified conditions do not exist!"))
  
  if (!is.logical(verbose))
    return(list("error"=TRUE, "message"="verbose must be a logical value!"))
  
  if ((!is.numeric(threads)) || threads < 1) {
    return(list("error"=TRUE, "message"="Invalid number of threads!"))
  } else if (threads > parallel::detectCores()) {
    return(list("error"=TRUE, "message"=paste("The system does not support that many threads! MAX available: ", parallel::detectCores())))
  }
  if ((!is.numeric(count_thresh)) || count_thresh < 0 )
    return(list("error"=TRUE, "message"="Invalid read-count threshold! Must be zero or a positive number."))
  
  if (! testmode %in% c("G-test", "g-test", "prop-test", "both"))
    return(list("error"=TRUE, "message"="Unrecognized value for testmode!"))
  
  if (! boots %in% c("G-test", "g-test","prop-test", "both", "none"))
    return(list("error"=TRUE, "message"="Unrecognized value for boots!"))
  
  if ((!is.numeric(bootnum)) || bootnum < 1)
    return(list("error"=TRUE, "message"="Invalid number of bootstraps! Must be a positive number."))
  
  tx <- slo$kal[[1]]$bootstrap[[1]][[BS_TARGET_COL]][ order(slo$kal[[1]]$bootstrap[[1]][[BS_TARGET_COL]]) ]
  for (k in 2:length(slo$kal)) {
    if (!all( tx == slo$kal[[k]]$bootstrap[[1]][[BS_TARGET_COL]][ order(slo$kal[[k]]$bootstrap[[1]][[BS_TARGET_COL]]) ] ))
      return(list("error"=TRUE, "message"="Inconsistent set of transcript IDs! Please try again, using the same annotation for all your samples!"))
  }
  
  if ((!is.numeric(dprop_thresh)) || dprop_thresh < 0 || dprop_thresh > 1)
    return(list("error"=TRUE, "message"="Invalid proportion difference threshold! Must be between 0 and 1."))
  
  return(list("error"=FALSE, "message"="All good!"))
}


#================================================================================
#' Initialise progress updates
#'
#' @param on Flag indicating whether updates are on (TRUE) or not (FALSE)
#' @return The progress update object
#'
init_progress <- function(on)
{
  progress_steps <- data.frame(c(1, 2, 3, 10, 11, 15, 30, 31, 100),
                               c("Checked parameters...",
                                 "Creating look-up structures...",
                                 "Extracting counts from bootstraps...",
                                 "Allocating output structure...",
                                 "Calculating counts statistics...",
                                 "Calculating p-values...",
                                 "Filling in more info...",
                                 "Bootstrapping p-values (if applicable)...",
                                 "All done!"),
                               stringsAsFactors = FALSE)
  progress <- TxtProgressUpdate(steps=progress_steps, on=on)
  return(progress)
}


#================================================================================
#' Group sample numbers by factor.
#'
#' @param covariates a dataframe with different factor variables.
#' @return list of lists (per covariate) of vectors (per factor level).
#'
#' Row number corresponds to sample number.
#'
group_samples <- function(covariates) {
  samplesByVariable <- list()
  for(varname in names(covariates)) {
    categories <- levels(as.factor(covariates[, varname]))
    samplesByVariable[[varname]] <- list()
    for (x in categories) {
      samplesByVariable[[varname]][[x]] <- which(covariates[, varname] == x)
    }
  }
  return(samplesByVariable)
}


#================================================================================
#' Extract bootstrap counts into a less nested structure.
#' 
#' NA replaced with 0. Means counts per transcripts are calculated and included as
#' a column per sample.
#' 
#' @param slo A sleuth object.
#' @param tx A vector of transcript ids. The results will be ordered according to this vector.
#' @param samples A numeric vector of samples to extract counts for.
#' @param COUNTS_COL The name of the column with the counts.
#' @param BS_TARGET_COL The name of the column with the transcript IDs.
#' @return a list of data.tables, one per sample, containing all the bootstrap counts.
#'
#' Transcripts in \code{slo} that are missing from \code{tx} will be skipped completely.
#' Transcripts in \code{tx} that are missing from \code{slo} are automatically padded with NA, which we re-assign as 0.
#'
denest_boots <- function(slo, tx, samples, COUNTS_COL, BS_TARGET_COL) {
  lapply(samples, function(smpl) {
    # Extract counts in the order of provided transcript vector, for safety and consistency.
    dt <- as.data.table( lapply(slo$kal[[smpl]]$bootstrap, function(boot) {
      roworder <- match(tx, boot[[BS_TARGET_COL]])
      boot[roworder, COUNTS_COL]
    }))
    # Replace any NAs with 0. Happens when annotation different from that used for DTE.
    dt[is.na(dt)] <- 0
    # Add transcript ID.
    dt[, "target_id" := tx]
  })
}


#================================================================================
#' Create output structure.
#' 
#' @param annot a dataframe with at least 2 columns: target_id and parent_id.
#' @param full Full-sized structure or core fields only. Either "full" or "short".
#' @return a list-based structure of class dtu.
#' 
alloc_out <- function(annot, full){
  if (full == "full") {
    Parameters <- list("var_name"=NA_character_, "cond_A"=NA_character_, "cond_B"=NA_character_,
                       "num_replic_A"=NA_integer_, "num_replic_B"=NA_integer_,
                       "p_thresh"=NA_real_, "count_thresh"=NA_real_, "dprop_thresh"=NA_real_,
                       "tests"=NA_character_, "bootstrap"=NA_character_, "bootnum"=NA_integer_, "threads"=NA_integer_)
    Genes <- data.table("parent_id"=levels(as.factor(annot$parent_id)),
                        "known_transc"=NA_integer_, "detect_transc"=NA_integer_, "elig_transc"=NA_integer_,
                        "eligible"=NA,                              # eligible for testing (reduce number of tests)
                        "Pt_DTU"=NA, "Gt_DTU"=NA, "Gt_dtuAB"=NA, "Gt_dtuBA"=NA,
                        "Gt_pvalAB"=NA_real_, "Gt_pvalBA"=NA_real_, "Gt_pvalAB_corr"=NA_real_, "Gt_pvalBA_corr"=NA_real_,
                        "Gt_boot_dtuAB"=NA_real_, "Gt_boot_dtuBA"=NA_real_,
                        "Gt_boot_meanAB"=NA_real_, "Gt_boot_meanBA"=NA_real_, "Gt_boot_stdevAB"=NA_real_, "Gt_boot_stdevBA"=NA_real_,
                        "Gt_boot_minAB"=NA_real_, "Gt_boot_minBA"=NA_real_, "Gt_boot_maxAB"=NA_real_, "Gt_boot_maxBA"=NA_real_,
                        "Gt_boot_na"=NA_integer_)
    Transcripts <- data.table("target_id"=annot$target_id, "parent_id"=annot$parent_id,
                              "propA"=NA_real_, "propB"=NA_real_, "Dprop"=NA_real_,
                              "eligible"=NA,                        # eligible for testing (reduce number of tests)
                              "Gt_DTU"=NA, "Pt_DTU"=NA,
                              "Pt_pval"=NA_real_,  "Pt_pval_corr"=NA_real_, 
                              "Pt_boot_dtu"=NA_real_, "Pt_boot_mean"=NA_real_, "Pt_boot_stdev"=NA_real_, "Pt_boot_min"=NA_real_,"Pt_boot_max"=NA_real_,
                              "Pt_boot_na"=NA_integer_,
                              "sumA"=NA_real_, "sumB"=NA_real_,      # sum across replicates of means across bootstraps
                              "meanA"=NA_real_, "meanB"=NA_real_,    # mean across replicates of means across bootstraps
                              "stdevA"=NA_real_, "stdevB"=NA_real_,  # standard deviation across replicates of means across bootstraps
                              "totalA"=NA_real_, "totalB"=NA_real_)  # sum of all transcripts for that gene
    
  } else {
    Parameters <- list("num_replic_A"=NA_integer_, "num_replic_B"=NA_integer_)
    Genes <- data.table("parent_id"=levels(as.factor(annot$parent_id)),
                        "elig_transc"=NA_integer_,
                        "eligible"=NA,                              # eligible for testing (reduce number of tests)
                        "Gt_pvalAB"=NA_real_, "Gt_pvalBA"=NA_real_,
                        "Gt_pvalAB_corr"=NA_real_, "Gt_pvalBA_corr"=NA_real_)
    Transcripts <- data.table("target_id"=annot$target_id, "parent_id"=annot$parent_id,
                              "propA"=NA_real_, "propB"=NA_real_, "Dprop"=NA_real_,
                              "eligible"=NA,                        # eligible for testing (reduce number of tests)
                              "Pt_pval"=NA_real_,  "Pt_pval_corr"=NA_real_,
                              "sumA"=NA_real_, "sumB"=NA_real_,      # sum across replicates of means across bootstraps
                              "totalA"=NA_real_, "totalB"=NA_real_)  # sum of all transcripts for that gene
    
  }
  setkey(Genes, parent_id)
  setkey(Transcripts, parent_id, target_id)
  return(list("Parameters"=Parameters, "Genes"=Genes, "Transcripts"=Transcripts))
}


#================================================================================
#' Set up and execute the tests.
#'
#' @param counts_A A data.table of counts for condition A. x: sample, y: transcript.
#' @param counts_B A data.table of counts for condition B. x: sample, y: transcript.
#' @param tx_filter A data.table with target_id and parent_id.
#' @param testmode One of "G-test", "prop-test", "both".
#' @param full Either "full" (for complete output structure) or "short" (for bootstrapping).
#' @param count_thresh Minimum number of counts per replicate.
#' @param dprop_thresh Minimum difference in proportions.
#' @param correction Multiple testing correction type.
#' @return list
#' 
#' @import data.table
#' 
calculate_DTU <- function(counts_A, counts_B, tx_filter, testmode, full, count_thresh, dprop_thresh, correction) {
  
  #---------- PRE-ALLOCATE
  
  # Pre-allocate results object.
  resobj <- alloc_out(tx_filter, full)
  resobj$Parameters["num_replic_A"] <- dim(counts_A)[2]
  resobj$Parameters["num_replic_B"] <- dim(counts_B)[2]
  
  #---------- STATS
  
  # Statistics per transcript across all bootstraps per condition, for filtered targets only.
  resobj$Transcripts[, sumA :=  rowSums(counts_A) ]
  resobj$Transcripts[, sumB :=  rowSums(counts_B) ]
   # Sums and proportions, for filtered targets only.
  resobj$Transcripts[, totalA := sum(sumA), by=parent_id]
  resobj$Transcripts[, totalB := sum(sumB), by=parent_id]
  resobj$Transcripts[, propA := sumA/totalA]
  resobj$Transcripts[, propB := sumB/totalB]
  resobj$Transcripts[, Dprop := propB - propA]
  
  #---------- FILTER
  
  # Filter transcripts and genes to reduce number of tests:
  ctA <- count_thresh * resobj$Parameters[["num_replic_A"]]  # Adjust count threshold for number of replicates.
  ctB <- count_thresh * resobj$Parameters[["num_replic_B"]]
  resobj$Transcripts[, eligible := (abs(Dprop) >= dprop_thresh  &  (sumA >= ctA | sumB >= ctB) )] 
  resobj$Transcripts[(is.na(eligible)), eligible := FALSE]
  resobj$Genes[, elig_transc := resobj$Transcripts[, .(parent_id, ifelse(eligible, 1, 0))][, sum(V2), by = parent_id][, V1] ]
  resobj$Genes[, eligible := elig_transc >= 2]
  
  #---------- TESTS
  
  # Proportion test.
  if ( any(testmode == c("prop-test", "both"))) {
    # The test function is not vectorised, nor easily vectorisable. Looping is needed.
    # Data tables allow values to be changed in place, and the table is pre-allocated. 
    # So access by row should be faster than looking up keys, and there should be no data-copying penalty.
    resobj$Transcripts[(eligible), Pt_pval := as.vector(apply(resobj$Transcripts[(eligible), .(sumA, sumB, totalA, totalB)], MARGIN = 1, 
                                                               FUN = function(row) { prop.test(x = row[c("sumA", "sumB")], n = row[c("totalA", "totalB")], correct = TRUE)[["p.value"]] } )) ]
    resobj$Transcripts[(eligible), Pt_pval_corr := p.adjust(Pt_pval, method=correction)]
  }
  
  # G test.
  if ( any(testmode == c("G-test", "g-test", "both"))) {
    # Set key to parents only. Was previously set to (parent, target) to ensure consistent order of targets for un-keyed manipulations.
    setkey(resobj$Transcripts, parent_id)
    resobj$Genes[(eligible), c("Gt_pvalAB", "Gt_pvalBA") := 
                   as.data.frame( t( as.data.frame( lapply(resobj$Genes[(eligible), parent_id], function(parent) {
                     # Extract all relevant data to avoid repeated look ups in the large table.
                     subdt <- resobj$Transcripts[parent, .(sumA, sumB, propA, propB)]
                     pAB <- g.test(x = subdt[, sumA], p = subdt[, propB])
                     pBA <- g.test(x = subdt[, sumB], p = subdt[, propA])
                     c(pAB, pBA) }) ))) ]
    resobj$Genes[, Gt_pvalAB_corr := p.adjust(Gt_pvalAB, method=correction)]
    resobj$Genes[, Gt_pvalBA_corr := p.adjust(Gt_pvalBA, method=correction)]
  }
  
  return(resobj)
}


#================================================================================
#' Log-likelihood test of goodness of fit.
#'
#' @param x	a numeric vector of positive numbers, with at least one non-zero value.
#' @param p	a vector of probabilities of the same length of x.
#'
#' Sourced and adapted from from:
#' V3.3 Pete Hurd Sept 29 2001. phurd@ualberta.ca
#' http://www.psych.ualberta.ca/~phurd/cruft/g.test.r
#'
g.test <- function(x, p = rep(1/length(x), length(x)))
{
  n <- sum(x)
  E <- n * p
  names(E) <- names(x)
  g <- 0
  for (i in 1:length(x)){
    if (x[i] != 0) g <- g + x[i] * log(x[i]/E[i])
  }
  q <- 1
  STATISTIC <- G <- 2*g/q
  PARAMETER <- length(x) - 1
  PVAL <- pchisq(STATISTIC, PARAMETER, lower = FALSE)
}

