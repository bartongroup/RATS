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
#' @param count_thresh Transcripts with fewer reads per replicate than this will be ignored (default 5).
#' @param testmode One of "G-test", "prop-test", "both" (default both).
#' @param correction The p-value correction to apply, as defined in \code{stats::p.adjust.methods}, default \code{"BH"}.
#' @param boots Bootstrap the p-values. Default TRUE.
#' @param threads Enable parallel processing. Default uses parallel::detectCores(). Try setting to 1 if you are having issues.
#' @param COUNTS_COL The name of the counts column to use for the DTU calculation (est_counts or tpm), default \code{"est_counts"}.
#' @param TARGET_COL The name of the transcript identifier column in the annot object, default \code{"target_id"}
#' @param PARENT_COL The name of the parent identifier column in the annot object, default \code{"parent_id"}.
#' @param BS_TARGET_COL The name of the transcript identifier column in the sleuth bootstrap tables, default \code{"target_id"}.
#' @return List of data tables, with gene-level and transcript-level information.
#'
#' @import data.table
#' @export
calculate_DTU <- function(slo, annot, name_A, name_B, varname="condition", 
                          p_thresh=0.05, count_thresh=5, testmode="both", correction="BH", 
                          verbose=FALSE, boots = TRUE, threads=parallel::detectCores(),
                          COUNTS_COL="est_counts", TARGET_COL="target_id", PARENT_COL="parent_id", BS_TARGET_COL="target_id") {
  #---------- PREP
  
  # Set up progress bar
  progress <- init_progress(verbose)
  
  progress <- update_progress(progress)
  # Input checks.
  threads <- as.integer(threads)  # Plain numbers default to double, unless integer R syntax is explicitly used.
  paramcheck <- parameters_good(slo, annot, name_A, name_B, varname, COUNTS_COL,
                                correction, p_thresh, TARGET_COL, PARENT_COL, BS_TARGET_COL, verbose, threads, count_thresh, testmode, boots)
  if (paramcheck$error) stop(paramcheck$message)
  
  #----------- LOOK-UP
  
  progress <- update_progress(progress)
  # Look-up from parent_id to target_id.
  targets_by_parent <- split(annot[[TARGET_COL]], annot[[PARENT_COL]])
  
  # Look-up from target_id to parent_id.
  tx_filter <- data.table(target_id = annot[[TARGET_COL]], parent_id = annot[[PARENT_COL]])
  setkey(tx_filter, target_id)
  
  # Reverse look-up from replicates to covariates.
  samples_by_condition <- group_samples(slo$sample_to_covariates)[[varname]]
  
  #---------- EXTRACT DATA
  
  progress <- update_progress(progress)
  # De-nest, average and index the counts from the bootstraps. 
  # For each condition, a list holds a dataframe for each replicate. 
  # Each dataframe holds the counts from ALL the bootstraps AND the mean counts across the bootstraps, indexed by transcript.
  # Transcripts found in the bootstraps but missing from the annotation are left out.
  data_A <- denest_boots(slo, tx_filter$target_id, samples_by_condition[[name_A]], COUNTS_COL, BS_TARGET_COL )
  data_B <- denest_boots(slo, tx_filter$target_id, samples_by_condition[[name_B]], COUNTS_COL, BS_TARGET_COL )
  bootmeans_A <- as.data.table(lapply(data_A, function(b) b[, mean_count]))
  bootmeans_B <- as.data.table(lapply(data_B, function(b) b[, mean_count]))
  
  #---------- PRE-ALLOCATE
  
  progress <- update_progress(progress)
  # Pre-allocate results object.
  resobj <- alloc_out(tx_filter)
  # Fill in some stuff I already know.
  resobj$Parameters["num_replic_A"] <- length(data_A)
  resobj$Parameters["num_replic_B"] <- length(data_B)
  resobj$Parameters["p_thresh"] <- p_thresh 
  resobj$Parameters["count_thresh"] <- count_thresh
  resobj$Parameters["tests"] <- testmode
  resobj$Parameters["threads"] <- threads
  resobj$Genes[, known_transc := resobj$Transcripts[, length(target_id), by=parent_id][, V1] ]  # V1 is the automatic column name for the lengths in the subsetted data.table
  
  #---------- STATS
  
  progress <- update_progress(progress)
  # Statistics per transcript across all bootstraps per condition, for filtered targets only.
  resobj$Transcripts[, sumA :=  rowSums(bootmeans_A) ]
  resobj$Transcripts[, sumB :=  rowSums(bootmeans_B) ]
  resobj$Transcripts[, meanA :=  rowMeans(bootmeans_A) ]
  resobj$Transcripts[, meanB :=  rowMeans(bootmeans_B) ]
  resobj$Transcripts[, stdevA :=  sqrt(matrixStats::rowVars(as.matrix(bootmeans_A))) ]
  resobj$Transcripts[, stdevB :=  sqrt(matrixStats::rowVars(as.matrix(bootmeans_B))) ]
  
  # Sums and proportions, for filtered targets only.
  resobj$Transcripts[, totalA := sum(sumA), by=parent_id]
  resobj$Transcripts[, totalB := sum(sumB), by=parent_id]
  resobj$Transcripts[, propA := sumA/totalA]
  resobj$Transcripts[, propB := sumB/totalB]
  resobj$Transcripts[, Dprop := propB - propA]
  
  # Filter transcripts and genes to reduce number of tests:
  resobj$Transcripts[, test_elig := (propA + propB > 0  &  
                                     propA + propB < 2  &  
                                     sumA > resobj$Parameters[["num_replic_A"]] * count_thresh  &  
                                     sumB > resobj$Parameters[["num_replic_B"]] * count_thresh)]
  resobj$Genes[, usable_transc := resobj$Transcripts[, .(parent_id, ifelse(test_elig, 1, 0))][, sum(V2), by = parent_id][, V1] ]  # Assign 1 for eligible transcripts, 0 for non, then sum to get the number of eligible trascripts per gene.
  resobj$Gene[, test_elig := usable_transc >= 2]

  #---------- TESTS
    
  progress <- update_progress(progress)
  if (testmode %in% c("prop-test", "both")) {
    #   # G test, only for parents and targets that survived filtering.
    #   # Compare B counts to A ratios:
    #   results$Genes[actual_parents, Gt_pvalAB := parallel::parSapply(wcl, actual_targets_by_parent, function(targets)
    #     g.test(results$Transcripts[targets, sumB],
    #            p=results$Transcripts[targets, propA])[["p.value"]])]
    #   # Compare A counts to B ratios:
    #   results$Genes[actual_parents, Gt_pvalBA := parallel::parSapply(wcl, actual_targets_by_parent, function(targets)
    #     g.test(results$Transcripts[targets, sumA],
    #            p=results$Transcripts[targets, propB])[["p.value"]])]
    #   # Correct p-values and apply threshold.
    #   results$Genes[, Gt_pvalAB_corr := p.adjust(Gt_pvalAB, method=correction)]
    #   results$Genes[, Gt_dtuAB := Gt_pvalAB_corr < p_thresh]
    #   results$Genes[, Gt_pvalBA_corr := p.adjust(Gt_pvalBA, method=correction)]
    #   results$Genes[, Gt_dtuBA := Gt_pvalBA_corr < p_thresh]
    #   # Find the agreements.
    #   results$Genes[, Gt_DTU := Gt_dtuAB & Gt_dtuBA ]
  }

  progress <- update_progress(progress)
  if (testmode %in% c("prop-test", "both")) {
    # Proportion test.
    # The test function is not vectorised, nor easily vectorisable. Looping is needed.
    # Data tables allow values to be changed in place, and the table is pre-allocated. 
    # So access by row should be faster than looking up keys, and there should be no data-copying penalty.
    resobj$Transcripts[(test_elig), Pt_pval := as.vector(apply(resobj$Transcripts[(test_elig), .(sumA, sumB, totalA, totalB)], MARGIN = 1, 
        FUN = function(row) { prop.test(x = row[c("sumA", "sumB")], n = row[c("totalA", "totalB")], correct = TRUE)[["p.value"]] } )) ]
    resobj$Transcripts[(test_elig), Pt_pval_corr := p.adjust(Pt_pval, method=correction)]
    resobj$Transcripts[, Pt_DTU := Pt_pval_corr < p_thresh]
    tmp <- resobj$Transcripts[, any(Pt_DTU), by = parent_id][, V1]
    
    #resobj$Genes[, Pt_DTU := tmp]
    
    
  }
  
  #---------- DONE
  
  progress <- update_progress(progress)
  return(resobj)
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
#' @param tx A vector of transcript ids.
#' @param samples A numeric vector of samples to extract counts for.
#' @param COUNTS_COL The name of the column with the counts.
#' @param BS_TARGET_COL The name of the column with the transcript IDs.
#' @return a list of data.tables, one per sample, containing all the bootstrap counts.
#'
denest_boots <- function(slo, tx, samples, COUNTS_COL, BS_TARGET_COL) {
  lapply(samples, function(smpl) {
    # Extract counts in the order of provided transcript vector, for safety and consistency.
    dt <- as.data.table( lapply(slo$kal[[smpl]]$bootstrap, function(boot) {
      roworder <- match(tx, boot[[BS_TARGET_COL]])
      boot[roworder, COUNTS_COL]
    }))
    # Do something about the ugly huge default names.
    nm <- names(dt)
    setnames(dt, nm, as.character(seq(1, length(nm), 1)))
    # Replace any NAs with 0. Happens when annotation different from that used for DTE.
    dt[is.na(dt)] <- 0
    # Add mean count and transcript ID.
    means <- rowMeans(dt)  # Compute it while the table is purely counts.
    dt[, c("mean_count", "target_id") := list(means, tx)]
  })
}


#================================================================================
#' Create output structure.
#' 
#' @param annot a dataframe with at least 2 columns: target_id and parent_id.
#' @return a list-based structure of class dtu.
#' 
alloc_out <- function(annot){
  Parameters <- list("var_name"=NA_character_, "cond_A"=NA_character_, "cond_B"=NA_character_,
                     "num_replic_A"=NA_integer_, "num_replic_B"=NA_integer_,
                     "p_thresh"=NA_real_, "count_thresh"=NA_real_, "tests"=NA_character_, threads=NA_integer_)
  Genes <- data.table("parent_id"=levels(as.factor(annot$parent_id)),
                      "known_transc"=NA_integer_, "usable_transc"=NA_integer_,
                      "test_elig"=NA,                              # eligible for testing (reduce number of tests)
                      "Pt_DTU"=NA, "Gt_DTU"=NA,
                      "Gt_pvalAB"=NA_real_, "Gt_pvalBA"=NA_real_,
                      "Gt_pvalAB_stdev"=NA_real_, "Gt_pvalBA_stdev"=NA_real_,
                      "Gt_pvalAB_corr"=NA_real_, "Gt_pvalBA_corr"=NA_real_,
                      "Gt_dtuAB"=NA, "Gt_dtuBA"=NA)
  Transcripts <- data.table("target_id"=annot$target_id, "parent_id"=annot$parent_id,
                            "propA"=NA_real_, "propB"=NA_real_, "Dprop"=NA_real_,
                            "test_elig"=NA,                        # eligible for testing (reduce number of tests)
                            "Gt_DTU"=NA, "Pt_DTU"=NA,
                            "Pt_pval"=NA_real_, "Pt_pval_stdev"=NA_real_,  "Pt_pval_corr"=NA_real_, 
                            "sumA"=NA_real_, "sumB"=NA_real_,      # sum across replicates of means across bootstraps
                            "meanA"=NA_real_, "meanB"=NA_real_,    # mean across replicates of means across bootstraps
                            "stdevA"=NA_real_, "stdevB"=NA_real_,  # standard deviation across replicates of means across bootstraps
                            "totalA"=NA_real_, "totalB"=NA_real_)  # sum of all transcripts for that gene
  setkey(Genes, parent_id)
  setkey(Transcripts, parent_id, target_id)
  return(structure(list("Parameters"=Parameters, "Genes"=Genes, "Transcripts"=Transcripts), class="dtu"))
}


#================================================================================
#' Check input parameters.
#'
#' @return List with a logical value and a message.
#'
parameters_good <- function(slo, annot, ref_name, comp_name, varname, COUNTS_COL,
                            correction, p_thresh, TARGET_COL, PARENT_COL, BS_TARGET_COL, verbose, 
                            threads, count_thresh, testmode, boots) {
  if ( ! is.data.frame(annot))
    return(list("error"=TRUE, "message"="transcripts is not a data.frame!"))
  if (any( ! c(TARGET_COL, PARENT_COL) %in% names(annot)))
    return(list("error"=TRUE, "message"="The specified target and parent IDs field-names do not exist in transcripts!"))
  if ( ! BS_TARGET_COL %in% names(slo$kal[[1]]$bootstrap[[1]]))
    return(list("error"=TRUE, "message"="The specified target IDs field-name does not exist in the bootstraps!"))
  if ( ! COUNTS_COL %in% names(slo$kal[[1]]$bootstrap[[1]]))
    return(list("error"=TRUE, "message"="The specified counts field-name does not exist!"))
  if ( ! correction %in% p.adjust.methods)
    return(list("error"=TRUE, "message"="Invalid p-value correction method name. Refer to stats::p.adjust.methods!"))
  if ( ( ! is.numeric(p_thresh)) || p_thresh > 1 || p_thresh < 0)
    return(list("error"=TRUE, "message"="Invalid p-value threshold!"))
  if ( ! varname %in% names(slo$sample_to_covariates))
    return(list("error"=TRUE, "message"="The specified covariate name does not exist!"))
  if ( any( ! c(ref_name, comp_name) %in% slo$sample_to_covariates[[varname]] ))
    return(list("error"=TRUE, "message"="One or both of the specified conditions do not exist!"))
  if ( ! is.logical(verbose))
    return(list("error"=TRUE, "message"="verbose must be a logical value!"))
  if ( ( ! is.numeric(threads)) || threads < 1) {
    return(list("error"=TRUE, "message"="Invalid number of threads!"))
  } else if (threads > parallel::detectCores()) {
    return(list("error"=TRUE, "message"=paste("The system does not support that many threads! MAX available: ", parallel::detectCores())))
  }
  if ( (! is.numeric(count_thresh)) || count_thresh < 0 )
    return(list("error"=TRUE, "message"="Invalid read-count threshold! Must be zero or a positive number."))
  if ( ! testmode %in% c("G-test", "proportion-test", "both"))
    return(list("error"=TRUE, "message"="Unrecognized value for testmode!"))
  if ( ! is.logical(boots))
    return(list("error"=TRUE, "message"="Boots must be a logical value."))
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
  progress_steps <- data.frame(c(1, 3, 7, 14, 15, 20, 60, 100),
                               c("Checking parameters...",
                                 "Creating look-up structures...",
                                 "Extracting counts from bootstraps...",
                                 "Allocating output structure...",
                                 "Calculating counts statistics...",
                                 "Calculating G-test p-values...",
                                 "Calculating proportions-test p-values...",
                                 "All done!"),
                               stringsAsFactors = FALSE)
  progress <- TxtProgressUpdate(steps=progress_steps, on=on)
  return(progress)
}

