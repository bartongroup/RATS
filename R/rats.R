#================================================================================
#' Calculate differential transcript usage.
#'
#' @param sleuth_data A sleuth object.
#' @param transcripts A dataframe matching the transcript identifiers to their corresponding gene identifiers.
#' @param name_A The name for one condition, as it appears in the \code{sample_to_covariates} table within the sleuth object.
#' @param name_B The name for the other condition, as it appears in the \code{sample_to_covariates} table within the sleuth object.
#' @param varname The name of the covariate to which the two conditions belong, as it appears in the \code{sample_to_covariates} table within the sleuth object. Default \code{"condition"}.
#' @param verbose Display progress updates, default \code{FALSE}.
#' @param p_thresh The p-value threshold, default 0.05.
#' @param count_thresh Transcripts with fewer reads will be ignored (default 5).
#' @param testmode One of "G-test", "proportion-test", "both" (default both).
#' @param correction The p-value correction to apply, as defined in \code{stats::p.adjust.methods}, default \code{"BH"}.
#' @param boots Bootstrap the p-values. Default TRUE.
#' @param threads Enable parallel processing. Default uses parallel::detectCores(). Try setting to 1 if you are having issues.
#' @param COUNTS_COL The name of the counts column to use for the DTU calculation (est_counts or tpm), default \code{"est_counts"}.
#' @param TARGET_COL The name of the transcript identifier column in the transcripts object, default \code{"target_id"}
#' @param PARENT_COL The name of the parent identifier column in the transcripts object, default \code{"parent_id"}.
#' @param BS_TARGET_COL The name of the transcript identifier column in the sleuth bootstrap tables, default \code{"target_id"}.
#' @return List of data tables, with gene-level and transcript-level information.
#'
#' @import data.table
#' @export
calculate_DTU <- function(sleuth_data, transcripts, name_A, name_B, varname="condition", 
                          p_thresh=0.05, count_thresh=5, testmode="both", correction="BH", 
                          verbose=FALSE, boots = TRUE, threads=parallel::detectCores(),
                          COUNTS_COL="est_counts", TARGET_COL="target_id", PARENT_COL="parent_id", BS_TARGET_COL="target_id") {
  #---------- PREP
  # Set up progress bar
  progress <- init_progress(verbose)
  
  progress <- update_progress(progress)
  # Input checks.
  threads <- as.integer(threads)  # Plain numbers default to double, unless integer R syntax is explicitly used.
  paramcheck <- parameters_good(sleuth_data, transcripts, name_A, name_B, varname, COUNTS_COL,
                                correction, p_thresh, TARGET_COL, PARENT_COL, BS_TARGET_COL, verbose, threads, count_thresh, testmode, boots)
  if (paramcheck$error) stop(paramcheck$message)
  
  progress <- update_progress(progress)
  # Initialize workers cluster.
  wcl <- parallel::makeCluster(threads, type="PSOCK")  # PSOCK is the most portable option but may get blocked by over-zealous security software.
  
  #----------- LOOK UPS
  progress <- update_progress(progress)
  # Look-up from parent_id to target_id
  targets_by_parent <- split(as.matrix(transcripts[TARGET_COL]), transcripts[[PARENT_COL]])
  
  # Look-up from target_id to parent_id
  tx_filter <- data.table(target_id = transcripts[[TARGET_COL]], parent_id = transcripts[[PARENT_COL]])
  setkey(tx_filter, target_id)
  
  # Reverse look-up from replicates to covariates.
  samples_by_condition <- group_samples(sleuth_data$sample_to_covariates)[[varname]]
  
  #---------- EXTRACT DATA
  progress <- update_progress(progress)
  # De-nest and index the counts from the bootstraps.
  data_A <- denest_boots(wcl, sleuth_data, tx_filter$target_id, samples_by_condition[[name_A]], COUNTS_COL, BS_TARGET_COL )
  data_B <- denest_boots(wcl, sleuth_data, tx_filter$target_id, samples_by_condition[[name_B]], COUNTS_COL, BS_TARGET_COL )
  
  #----------
  progress <- update_progress(progress)
  # Pre-allocate output structure.
  results <- list("Parameters"=list("var_name"=varname, "cond_A"=name_A, "cond_B"=name_B,
                                    "num_replic_A"=length(data_A), "num_replic_B"=length(data_B),
                                    "p_thresh"=p_thresh, "count_thresh"=count_thresh, "tests"=testmode),
                  "Genes"=data.table("parent_id"=levels(as.factor(tx_filter[[PARENT_COL]])),
                                     "known_transc"=NA_integer_, "usable_transc"=NA_integer_,
                                     "Gt_DTU"=NA, "Pt_DTU"=NA,
                                     "Gt_pvalAB"=NA_real_, "Gt_pvalBA"=NA_real_,
                                     "Gt_pvalAB_stdev"=NA_real_, "Gt_pvalBA_stdev"=NA_real_,
                                     "Gt_pvalAB_corr"=NA_real_, "Gt_pvalBA_corr"=NA_real_,
                                     "Gt_dtuAB"=NA, "Gt_dtuBA"=NA),
                  "Transcripts"=data.table("target_id"=tx_filter[[TARGET_COL]], "parent_id"=tx_filter[[PARENT_COL]],
                                           "propA"=NA_real_, "propB"=NA_real_,   # proportion of sums across replicates
                                           "Dprop"=NA_real_,                       # propB - propA
                                           "Gt_DTU"=NA, "Pt_DTU"=NA,
                                           "Pt_pval"=NA_real_, "Pt_pval_stdev"=NA_real_,  # proportion test
                                           "Pt_pval_corr"=NA_real_, 
                                           "sumA"=NA_real_, "sumB"=NA_real_,     # sum across replicates of means across bootstraps
                                           "meanA"=NA_real_, "meanB"=NA_real_,   # mean across replicates of means across bootstraps
                                           "stdevA"=NA_real_, "stdevB"=NA_real_,  # standard deviation across replicates of means across bootstraps
                                           "totalA"=NA_real_, "totalB"=NA_real_))  # sum of all transcripts for that gene
  
  
  setkey(results$Genes, parent_id)
  setkey(results$Transcripts, target_id)
  results$Genes[, known_transc := parallel::parSapply(wcl, results$Genes[[PARENT_COL]], function(p) length(targets_by_parent[[p]]))]
  
#   progress <- update_progress(progress)
#   # Statistics per transcript across all bootstraps per condition, for filtered targets only.
#   results$Transcripts[actual_targets, sumA :=  rowSums(count_data[[name_A]])]
#   results$Transcripts[actual_targets, sumB :=  rowSums(count_data[[name_B]])]
#   results$Transcripts[actual_targets, meanA :=  rowMeans(count_data[[name_A]])]
#   results$Transcripts[actual_targets, meanB :=  rowMeans(count_data[[name_B]])]
#   results$Transcripts[actual_targets, stdevA :=  sqrt(matrixStats::rowVars(as.matrix(count_data[[name_A]]))) ]
#   results$Transcripts[actual_targets, stdevB :=  sqrt(matrixStats::rowVars(as.matrix(count_data[[name_B]]))) ]
#   
#   # Sums and proportions, for filtered targets only.
#   results$Transcripts[actual_targets, totalA := sum(sumA), by=parent_id]
#   results$Transcripts[actual_targets, totalB := sum(sumB), by=parent_id]
#   results$Transcripts[actual_targets, propA := sumA/totalA]
#   results$Transcripts[actual_targets, propB := sumB/totalB]
#   
#   progress <- update_progress(progress)
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
#   
#   progress <- update_progress(progress)
#   # Proportion test.
#   results$Transcripts[actual_targets, Pt_pval:= unlist(parallel::parLapply(wcl, actual_targets, function(target) 
#     prop.test(x=c(results$Transcripts[target, sumA], results$Transcripts[target, sumB]), 
#               n=c(results$Transcripts[target, totalA], results$Transcripts[target, totalB]), 
#               correct=TRUE
#     )[["p.value"]] )) ]
#   results$Transcripts[actual_targets, Pt_pval_corr:= p.adjust(results$Transcripts[actual_targets, Pt_pval], method=correction)]
#   results$Transcripts[, Pt_DTU := Pt_pval_corr < p_thresh]
  
  # Dismiss workers.
  parallel::stopCluster(wcl)
  
  progress <- update_progress(progress)
  return(results)
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
#' NA replaced with 0.
#' 
#' @param wcl A cluster of workers (threads).
#' @param slo A sleuth object.
#' @param tx A vector of transcript ids.
#' @param samples A numeric vector of samples to extract counts for.
#' @param COUNTS_COL The name of the column with the counts.
#' @param BS_TARGET_COL The name of the column with the transcript IDs.
#' @return a list of data.tables, one per sample, containing all the bootstrap counts.
#'
denest_boots <- function(wcl, slo, tx, samples, COUNTS_COL, BS_TARGET_COL) {
  lapply(samples, function(smpl) {
    dt <- as.data.frame( parallel::parLapply(wcl, slo$kal[[smpl]]$bootstrap, function(boot) {
      roworder <- match(tx, boot[[BS_TARGET_COL]])
      boot[roworder, COUNTS_COL]
    }))
    # Do something about the ugly huge default names.
    nm <- names(dt)
    names(dt) <- seq(1, length(nm), 1)
    # replace any NAs with 0. Happens when annotation different from that used for DTE.
    dt[is.na(dt),] <- 0
    # Add target_id as index.
    dt <- data.table(dt)
    dt[, target_id := tx]
    setkey(dt, target_id)
  })
}


#================================================================================
#' Check input parameters.
#'
#' @return List with a logical value and a message.
#'
parameters_good <- function(sleuth_data, transcripts, ref_name, comp_name, varname, COUNTS_COL,
                            correction, p_thresh, TARGET_COL, PARENT_COL, BS_TARGET_COL, verbose, 
                            threads, count_thresh, testmode, boots) {
  if ( ! is.data.frame(transcripts))
    return(list("error"=TRUE, "message"="transcripts is not a data.frame!"))
  if (any( ! c(TARGET_COL, PARENT_COL) %in% names(transcripts)))
    return(list("error"=TRUE, "message"="The specified target and parent IDs field-names do not exist in transcripts!"))
  if ( ! BS_TARGET_COL %in% names(sleuth_data$kal[[1]]$bootstrap[[1]]))
    return(list("error"=TRUE, "message"="The specified target IDs field-name does not exist in the bootstraps!"))
  if ( ! COUNTS_COL %in% names(sleuth_data$kal[[1]]$bootstrap[[1]]))
    return(list("error"=TRUE, "message"="The specified counts field-name does not exist!"))
  if ( ! correction %in% p.adjust.methods)
    return(list("error"=TRUE, "message"="Invalid p-value correction method name. Refer to stats::p.adjust.methods!"))
  if ( ( ! is.numeric(p_thresh)) || p_thresh > 1 || p_thresh < 0)
    return(list("error"=TRUE, "message"="Invalid p-value threshold!"))
  if ( ! varname %in% names(sleuth_data$sample_to_covariates))
    return(list("error"=TRUE, "message"="The specified covariate name does not exist!"))
  if ( any( ! c(ref_name, comp_name) %in% sleuth_data$sample_to_covariates[[varname]] ))
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
  progress_steps <- data.frame(c(1, 3, 7, 10, 20, 25, 35, 70, 100),
                               c("Checking parameters...",
                                 "Initializing threads...",
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

