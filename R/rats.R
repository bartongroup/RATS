#================================================================================
#' Calculate differential transcript usage.
#'
#' @param sleuth_data A sleuth object.
#' @param transcripts A dataframe matching the transcript identifiers to their corresponding gene identifiers.
#' @param name_A The name for one condition, as it appears in the \code{sample_to_covariates} table within the sleuth object.
#' @param name_B The name for the other condition, as it appears in the \code{sample_to_covariates} table within the sleuth object.
#' @param varname The name of the covariate to which the two conditions belong, as it appears in the \code{sample_to_covariates} table within the sleuth object. Default \code{"condition"}.
#' @param counts_col The name of the counts column to use for the DTU calculation (est_counts or tpm), default \code{"est_counts"}.
#' @param verbose Display progress updates, default \code{FALSE}.
#' @param p_thresh The p-value threshold, default 0.05.
#' @param count_thresh Transcripts with fewer reads will be ignored (default 5).
#' @param testmode One of "G-test", "proportion-test", "both" (default both).
#' @param correction The p-value correction to apply, as defined in \code{stats::p.adjust.methods}, default \code{"BH"}.
#' @param threads Enable parallel processing. Default uses parallel::detectCores(). Try setting to 1 if you are having issues.
#' @param TARGET_ID The name of the transcript identifier column in the transcripts object, default \code{"target_id"}
#' @param PARENT_ID The name of the parent identifier column in the transcripts object, default \code{"parent_id"}.
#' @param BS_TARGET_ID The name of the transcript identifier column in the sleuth bootstrap tables, default \code{"target_id"}.
#' @return List of data tables, with gene-level and transcript-level information.
#'
#' @export
#' @import data.table
calculate_DTU <- function(sleuth_data, transcripts, name_A, name_B,
                          varname="condition", counts_col="est_counts",
                          p_thresh=0.05, count_thresh=5, testmode="both", correction="BH", 
                          verbose=FALSE, threads=parallel::detectCores(),
                          TARGET_ID="target_id", PARENT_ID="parent_id", BS_TARGET_ID="target_id") {
  # Set up progress bar
  progress <- init_progress(verbose)
  
  progress <- update_progress(progress)
  # Input checks.
  threads <- as.integer(threads)  # Plain numbers default to double, unless integer R syntax is explicitly used.
  paramcheck <- parameters_good(sleuth_data, transcripts, name_A, name_B, varname, counts_col,
                                correction, p_thresh, TARGET_ID, PARENT_ID, BS_TARGET_ID, verbose, threads, count_thresh, testmode)
  if (paramcheck$error) stop(paramcheck$message)
  
  progress <- update_progress(progress)
  # Initialize workers cluster.
  wcl <- parallel::makeCluster(threads, type="PSOCK")  # PSOCK is the most portable option but may get blocked by over-zealous security software.
  
  progress <- update_progress(progress)
  # Look-up from parent_id to target_id
  targets_by_parent <- split(as.matrix(transcripts[TARGET_ID]), transcripts[[PARENT_ID]])
  
  progress <- update_progress(progress)
  # Identify genes with a single transcript. Order by gene ID and transcript ID.
  tx_filter <- transcripts[order(transcripts[[PARENT_ID]], transcripts[[TARGET_ID]]), ]
  tx_filter["has_siblings"] <- TRUE
  
  progress <- update_progress(progress)
  # Reverse look-up from replicates to covariates.
  samples_by_condition <- group_samples(sleuth_data$sample_to_covariates)[[varname]]
  
  progress <- update_progress(progress)
  # Build list of dataframes, one for each condition.
  # Each dataframe contains filtered and correctly ordered mean counts per sample from the bootstraps
  count_data <- parallel::parLapply(wcl, samples_by_condition, function(condition) make_filtered_bootstraps(sleuth_data, condition, tx_filter, counts_col, TARGET_ID, BS_TARGET_ID))
  
  progress <- update_progress(progress)
  # Remove entries which are entirely 0 across all conditions.
  nonzero <-  parallel::parLapply(wcl, count_data, function(condition) apply(condition, 1, function(row) !all(row == 0 )))
  count_data <- parallel::parLapply(wcl, count_data, function(condition) condition[Reduce("&", nonzero),, drop=FALSE])
  
  progress <- update_progress(progress)
  # Which IDs am I actually working with after the filters?
  actual_targets <- rownames(count_data[[name_A]])
  actual_parents <- levels(as.factor(tx_filter[[PARENT_ID]][match(actual_targets, tx_filter[[TARGET_ID]])]))
  actual_txs <- transcripts[transcripts[[TARGET_ID]] %in% actual_targets,]
  actual_targets_by_parent <- (split(as.matrix(actual_txs[TARGET_ID]), actual_txs[[PARENT_ID]]))[actual_parents]
  # Reject parents that now are left with a single child, as g.test() won't accept them.
  actual_targets_by_parent <- actual_targets_by_parent[parallel::parSapply(wcl, actual_targets_by_parent, function(targets) length(targets) > 1)]
  actual_parents <- names(actual_targets_by_parent)
  
  progress <- update_progress(progress)
  # Pre-allocate output structure.
  results <- list("Parameters"=list("var_name"=varname, "cond_A"=name_A, "cond_B"=name_B,
                                    "replicates_A"=dim(count_data[[name_A]])[2], "replicates_B"=dim(count_data[[name_B]])[2],
                                    "p_thresh"=p_thresh, "count_thresh"=count_thresh, "tests"=testmode),
                  "Genes"=data.table("parent_id"=levels(as.factor(tx_filter[[PARENT_ID]])),
                                     "known_transc"=NA_integer_, "usable_transc"=NA_integer_,
                                     "Gt_DTU"=NA, "Pt_DTU"=NA,
                                     "Gt_pvalAB"=NA_real_, "Gt_pvalBA"=NA_real_,
                                     "Gt_pvalAB_stdev"=NA_real_, "Gt_pvalBA_stdev"=NA_real_,
                                     "Gt_pvalAB_corr"=NA_real_, "Gt_pvalBA_corr"=NA_real_,
                                     "Gt_dtuAB"=NA, "Gt_dtuBA"=NA),
                  "Transcripts"=data.table("target_id"=tx_filter[[TARGET_ID]], "parent_id"=tx_filter[[PARENT_ID]],
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
  results$Genes[, known_transc := parallel::parSapply(wcl, results$Genes[[PARENT_ID]], function(p) length(targets_by_parent[[p]]))]
  results$Genes[, usable_transc := parallel::parSapply(wcl, results$Genes[[PARENT_ID]], function(p) ifelse(any(actual_parents == p), length(actual_targets_by_parent[[p]]), 0))]
  
  progress <- update_progress(progress)
  # Statistics per transcript across all bootstraps per condition, for filtered targets only.
  results$Transcripts[actual_targets, sumA :=  rowSums(count_data[[name_A]])]
  results$Transcripts[actual_targets, sumB :=  rowSums(count_data[[name_B]])]
  results$Transcripts[actual_targets, meanA :=  rowMeans(count_data[[name_A]])]
  results$Transcripts[actual_targets, meanB :=  rowMeans(count_data[[name_B]])]
  results$Transcripts[actual_targets, stdevA :=  sqrt(matrixStats::rowVars(as.matrix(count_data[[name_A]]))) ]
  results$Transcripts[actual_targets, stdevB :=  sqrt(matrixStats::rowVars(as.matrix(count_data[[name_B]]))) ]
  
  # Sums and proportions, for filtered targets only.
  results$Transcripts[actual_targets, totalA := sum(sumA), by=parent_id]
  results$Transcripts[actual_targets, totalB := sum(sumB), by=parent_id]
  results$Transcripts[actual_targets, propA := sumA/totalA]
  results$Transcripts[actual_targets, propB := sumB/totalB]
  
  progress <- update_progress(progress)
  # G test, only for parents and targets that survived filtering.
  # Compare B counts to A ratios:
  results$Genes[actual_parents, Gt_pvalAB := parallel::parSapply(wcl, actual_targets_by_parent, function(targets)
    g.test(results$Transcripts[targets, sumB],
           p=results$Transcripts[targets, propA])[["p.value"]])]
  # Compare A counts to B ratios:
  results$Genes[actual_parents, Gt_pvalBA := parallel::parSapply(wcl, actual_targets_by_parent, function(targets)
    g.test(results$Transcripts[targets, sumA],
           p=results$Transcripts[targets, propB])[["p.value"]])]
  # Correct p-values and apply threshold.
  results$Genes[, Gt_pvalAB_corr := p.adjust(Gt_pvalAB, method=correction)]
  results$Genes[, Gt_dtuAB := Gt_pvalAB_corr < p_thresh]
  results$Genes[, Gt_pvalBA_corr := p.adjust(Gt_pvalBA, method=correction)]
  results$Genes[, Gt_dtuBA := Gt_pvalBA_corr < p_thresh]
  # Find the agreements.
  results$Genes[, Gt_DTU := Gt_dtuAB & Gt_dtuBA ]
  
  progress <- update_progress(progress)
  # Proportion test.
  results$Transcripts[actual_targets, Pt_pval:= unlist(parallel::parLapply(wcl, actual_targets, function(target) 
    prop.test(x=c(results$Transcripts[target, sumA], results$Transcripts[target, sumB]), 
              n=c(results$Transcripts[target, totalA], results$Transcripts[target, totalB]), 
              correct=TRUE
          )[["p.value"]] )) ]
  results$Transcripts[actual_targets, Pt_pval_corr:= p.adjust(results$Transcripts[actual_targets, Pt_pval], method=correction)]
  results$Transcripts[, Pt_DTU := Pt_pval_corr < p_thresh]
  
  # Dismiss workers.
  parallel::stopCluster(wcl)
  
  progress <- update_progress(progress)
  return(results)
}






#================================================================================
#' Compute a logical vector marking as FALSE the single-target parents in a data frame.
#'
#' @param ids a data frame with at least two variables, \code{target_id} & \code{parent_id}.
#' @param p2t a list of vectors, listing the \code{target_id}s per \code{parent_id}.
#' @param TARGET_ID The name of transcript id column in transcripts object.
#' @param PARENT_ID The name of parent id column in transcripts object.
#' @return data.frame An updated version of the input ids.
#'
mark_sibling_targets <- function(ids, p2t, TARGET_ID, PARENT_ID) {
  rownames(ids) <- ids[[TARGET_ID]]
  
  # function testing for length > 1, for use in aggregate
  f <- function(x) { length(x) > 1 }
  # build has_siblings column by PARENT_ID by aggregating count of target ids
  # and testing if count > 1
  has_siblings <- aggregate(ids[TARGET_ID], by=ids[PARENT_ID], FUN=f)
  
  colnames(has_siblings)[2] <- "has_siblings"
  
  # inner join has_siblings to ids to give has_sibllings value for each target_id
  ids <- merge(ids, has_siblings, by=PARENT_ID)
  
  return(ids[order(ids[[PARENT_ID]], ids[[TARGET_ID]]), ])
}

#--------------------------------------------------------------------------------
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
#' For a given condition in the sleuth object, construct a dataframe containing counts from each bootstrap,
#' filtered according to tx_filter, and ordered according to target_ids
#'
#' @param sleuth_data A sleuth object.
#' @param condition A vector of sample numbers.
#' @param tx_filter A dataframe containing \code{target_id} and \code{has_siblings}.
#' @param counts_col The sleuth column name for the type of counts to use.
#' @param TARGET_ID The name of transcript id column in transcripts object.
#' @param BS_TARGET_ID The name of transcript id column in sleuth bootstrap tables.
#' @return A dataframe containing the counts from all bootstraps of all the samples for the condition.
#'
make_filtered_bootstraps <- function(sleuth_data, condition, tx_filter, counts_col, TARGET_ID, BS_TARGET_ID) {
  
  # make a list of dataframes, one df for each condition, containing the counts from its bootstraps
  count_data <- as.data.frame(lapply(condition, function(sample)
    rowMeans(sapply(sleuth_data$kal[[sample]]$bootstrap, function(e)
      filter_and_match(e, tx_filter, counts_col, TARGET_ID, BS_TARGET_ID) )) ))
  
  # now set the filtered target ids as rownames - previous call returns target ids in this order
  rownames(count_data) <- tx_filter[[TARGET_ID]][tx_filter$has_siblings]
  
  # replace any NAs with 0: some bootstrap did not have an entry for the transcript for the row
  count_data[is.na(count_data)] <- 0
  
  return(count_data)
}

#--------------------------------------------------------------------------------
#' Extract a filtered column from a bootstrap table
#'
#' @param bootstrap A bootstrap dataframe
#' @param tx_filter A boolean vector for filtering the bootstrap, with matching transcript ids
#' @param counts_col The column in the bootstrap table to extract
#' @param TARGET_ID The name of transcript id column in transcripts object.
#' @param BS_TARGET_ID The name of transcript id column in sleuth bootstrap tables.
#' @return The column from the bootstrap table
#'
filter_and_match <- function(bootstrap, tx_filter, counts_col, TARGET_ID, BS_TARGET_ID)
{
  # create map from bootstrap to filter(i.e. main annotation) target ids
  b_to_f_rows <- match(tx_filter[[TARGET_ID]], bootstrap[[BS_TARGET_ID]])
  
  # map the bootstrap to the filter target ids and then apply the filter
  result <- (bootstrap[b_to_f_rows, counts_col]) [tx_filter$has_siblings]
  
  return(result)
}

#================================================================================
#' Check input parameters.
#'
#' @return List with a logical value and a message.
#'
parameters_good <- function(sleuth_data, transcripts, ref_name, comp_name, varname, counts_col,
                            correction, p_thresh, TARGET_ID, PARENT_ID, BS_TARGET_ID, verbose, 
                            threads, count_thresh, testmode) {
  if ( ! is.data.frame(transcripts))
    return(list("error"=TRUE, "message"="transcripts is not a data.frame!"))
  if (any( ! c(TARGET_ID, PARENT_ID) %in% names(transcripts)))
    return(list("error"=TRUE, "message"="The specified target and parent IDs field-names do not exist in transcripts!"))
  if ( ! BS_TARGET_ID %in% names(sleuth_data$kal[[1]]$bootstrap[[1]]))
    return(list("error"=TRUE, "message"="The specified target IDs field-name does not exist in the bootstraps!"))
  if ( ! counts_col %in% names(sleuth_data$kal[[1]]$bootstrap[[1]]))
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
  if ( (! is.numeric(count_thresh)) || count_thresh < 0 ) {
    return(list("error"=TRUE, "message"="Invalid read-count threshold! Must be zero or a positive number."))
  }
  if ( ! testmode %in% c("G-test", "proportion-test", "both")){
    return(list("error"=TRUE, "message"="Unrecognized value for testmode!"))
  }
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
  progress_steps <- data.frame(c(1, 2, 4, 6, 8, 10, 20, 25, 27, 30, 35, 70, 100),
                               c("Checking parameters...",
                                 "Initializing threads...",
                                 "Mapping genes to transcripts...",
                                 "Identifying genes with multiple transcripts...",
                                 "Grouping samples by condition...",
                                 "Extracting counts from bootstraps...",
                                 "Removing all-zero count cases...",
                                 "Allocating output structure...",
                                 "Identifying genes and transcripts represented in the data...",
                                 "Calculating counts statistics...",
                                 "Calculating G-test p-values...",
                                 "Calculating proportions-test p-values...",
                                 "All done!"),
                               stringsAsFactors = FALSE)
  progress <- TxtProgressUpdate(steps=progress_steps, on=on)
  return(progress)
}

