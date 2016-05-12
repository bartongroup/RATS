#================================================================================
#' Calculate differential transcript usage.
#'
#' @param sleuth_data A sleuth object
#' @param transcripts A dataframe matching the transcript IDs to their corresponding gene IDs.
#' @param name_A The sleuth name for one condition.
#' @param name_B The sleuth name for the other condition.
#' @param varname The sleuth name of the covariate to which the two conditions belong, default \code{"condition"}.
#' @param counts_col The sleuth column to use for the calculation (est_counts or tpm), default \code{"est_counts"}.
#' @param correction p-value correction to apply, as defined in \code{stats::p.adjust.methods}, default \code{"BH"}.
#' @param p_thresh p-value threshold, default 0.05.
#' @param TARGET_ID The name of the transcript id column in transcripts object, default \code{"target_id"}
#' @param PARENT_ID The name of the parent id column in transcripts object, default \code{"parent_id"}.
#' @param BS_TARGET_ID The name of the transcript id column in sleuth bootstrap tables, default \code{"target_id"}.
#' @param verbose Whether to update progress updates, default \code{FALSE}.
#' @return List of data frames, with gene-level and transcript-level information.
#'
#' @export
#' @import data.table
calculate_DTU <- function(sleuth_data, transcripts, name_A, name_B,
                          varname="condition", counts_col="est_counts", correction="BH", p_thresh=0.05,
                          TARGET_ID="target_id", PARENT_ID="parent_id", BS_TARGET_ID="target_id",
                          verbose=FALSE)
{
  # Input checks.
  paramcheck <- parameters_good(sleuth_data, transcripts, name_A, name_B, varname, counts_col,
                                correction, p_thresh, TARGET_ID, PARENT_ID, BS_TARGET_ID, verbose)
  if (paramcheck$error) stop(paramcheck$message)

  # Set up progress bar
  progress <- init_progress(verbose)

  # Look-up from parent_id to target_id
  targets_by_parent <- split(as.matrix(transcripts[TARGET_ID]), transcripts[[PARENT_ID]])
  progress <- update(progress)

  # Identify genes with a single transcript. Order by gene ID and transcript ID.
#  tx_filter <- mark_sibling_targets(transcripts, targets_by_parent, TARGET_ID, PARENT_ID)
  tx_filter <- transcripts[order(transcripts[[PARENT_ID]], transcripts[[TARGET_ID]]), ]
  tx_filter["has_siblings"] <- TRUE
  progress <- update(progress)

  # Reverse look-up from replicates to covariates.
  samples_by_condition <- group_samples(sleuth_data$sample_to_covariates)[[varname]]
  progress <- update(progress)

  # build list of dataframes, one for each condition
  # each dataframe contains filtered and correctly ordered counts from all the bootstraps for the condition
  count_data <- lapply(samples_by_condition, function(condition) make_filtered_bootstraps(sleuth_data, condition, tx_filter, counts_col, TARGET_ID, BS_TARGET_ID))
  progress <- update(progress)

  # remove entries which are entirely 0 across all conditions
  nonzero <-  lapply(count_data, function(condition) apply(condition, 1, function(row) !all(row == 0 )))
  count_data <- lapply(count_data, function(condition) condition[Reduce("&", nonzero),, drop=FALSE])
  progress <- update(progress)

  # Which IDs am I actually working with after the filters?
  actual_targets <- rownames(count_data[[name_A]])
  actual_parents <- levels(as.factor(tx_filter[[PARENT_ID]][match(actual_targets, tx_filter[[TARGET_ID]])]))
  actual_targets_by_parent <- lapply(actual_parents, function(p) {                                                                # BOTTLENECK
    targets_by_parent[[p]][targets_by_parent[[p]] %in% actual_targets]  # the transcripts for which we have non-zero counts.      # BOTTLENECK
  })
  names(actual_targets_by_parent) <- actual_parents
  # Reject parents that now are left with a single child, as g.test() won't accept them.
  actual_targets_by_parent <- actual_targets_by_parent[sapply(actual_targets_by_parent, function(targets) length(targets) > 1)]
  actual_parents <- names(actual_targets_by_parent)
  progress <- update(progress)

  # Pre-allocate output structure.
  results <- list("Parameters"=list("var_name"=varname, "cond_A"=name_A, "cond_B"=name_B,
                                    "replicates_A"=dim(count_data[[name_A]])[2], "replicates_B"=dim(count_data[[name_B]])[2],
                                    "p_thresh"=p_thresh),
                  "Genes"=data.table("parent_id"=levels(as.factor(tx_filter[[PARENT_ID]])),
                                     "known_transc"=NA_integer_, "usable_transc"=NA_integer_,
                                     "pval_AB"=NA_real_, "pval_BA"=NA_real_,
                                     "pval_AB_corr"=NA_real_, "pval_BA_corr"=NA_real_,
                                     "dtu_AB"=NA, "dtu_BA"=NA, "dtu"=NA),
                  "Transcripts"=data.table("target_id"=tx_filter[[TARGET_ID]], "parent_id"=tx_filter[[PARENT_ID]],
                                           "prop_A"=NA_real_, "prop_B"=NA_real_,   # proportion of sums across replicates
                                           "sum_A"=NA_real_, "sum_B"=NA_real_,     # sum across replicates of means across bootstraps
                                           "mean_A"=NA_real_, "mean_B"=NA_real_,   # mean across replicates of means across bootstraps
                                           "var_A"=NA_real_, "var_B"=NA_real_))    # var across replicates of means across bootstraps
  setkey(results$Genes, parent_id)
  setkey(results$Transcripts, target_id)
  results$Genes[, known_transc := sapply(results$Genes[[PARENT_ID]], function(p) length(targets_by_parent[[p]]))]
  results$Genes[, usable_transc := sapply(results$Genes[[PARENT_ID]], function(p) ifelse(any(actual_parents == p), length(actual_targets_by_parent[[p]]), 0))]
  progress <- update(progress)

  # Statistics per transcript across all bootstraps per condition, for filtered targets only.
  results$Transcripts[actual_targets, sum_A :=  rowSums(count_data[[name_A]])]                         # TODO change to hybrid
  results$Transcripts[actual_targets, sum_B :=  rowSums(count_data[[name_B]])]                         # TODO change to hybrid
  results$Transcripts[actual_targets, mean_A :=  rowMeans(count_data[[name_A]])]                       # TODO change to hybrid
  results$Transcripts[actual_targets, mean_B :=  rowMeans(count_data[[name_B]])]                       # TODO change to hybrid
  results$Transcripts[actual_targets, var_A :=  matrixStats::rowVars(as.matrix(count_data[[name_A]]))] # TODO change to hybrid
  results$Transcripts[actual_targets, var_B :=  matrixStats::rowVars(as.matrix(count_data[[name_B]]))] # TODO change to hybrid
  progress <- update(progress)

  # Proportions = sum of tx / sum(sums of all related txs), for filtered targets only.
  results$Transcripts[, prop_A := sum_A/sum(sum_A), by=parent_id]
  results$Transcripts[, prop_B := sum_B/sum(sum_B), by=parent_id]
  progress <- update(progress)

  # P values, only for parents and targets that survived filtering.
  # Compare B counts to A ratios:
  results$Genes[actual_parents, pval_AB := sapply(actual_targets_by_parent, function(targets)
                                                  g.test(results$Transcripts[targets, sum_B],
                                                         p=results$Transcripts[targets, prop_A])[["p.value"]])]
  results$Genes[, pval_AB_corr := p.adjust(pval_AB, method=correction)]
  results$Genes[, dtu_AB := pval_AB < p_thresh]
  # Compare A counts to B ratios:
  results$Genes[actual_parents, pval_BA := sapply(actual_targets_by_parent, function(targets)
                                                  g.test(results$Transcripts[targets, sum_A],
                                                         p=results$Transcripts[targets, prop_B])[["p.value"]])]
  results$Genes[, pval_BA_corr := p.adjust(pval_BA, method=correction)]
  results$Genes[, dtu_BA := pval_BA < p_thresh]
  # Find the agreements.
  results$Genes[, dtu := dtu_AB & dtu_BA ]
  progress <- update(progress)

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
                            correction, p_thresh, TARGET_ID, PARENT_ID, BS_TARGET_ID, verbose) {
  if ( ! is.data.frame(transcripts))
    return(list("error"=TRUE, "message"="transcripts is not a data.frame."))
  if (any( ! c(TARGET_ID, PARENT_ID) %in% names(transcripts)))
    return(list("error"=TRUE, "message"="The specified target and parent IDs field-names do not exist in transcripts."))
  if ( ! BS_TARGET_ID %in% names(sleuth_data$kal[[1]]$bootstrap[[1]]))
    return(list("error"=TRUE, "message"="The specified target IDs field-name does not exist in the bootstraps."))
  if ( ! counts_col %in% names(sleuth_data$kal[[1]]$bootstrap[[1]]))
    return(list("error"=TRUE, "message"="The specified counts field-name does not exist."))
  if ( ! correction %in% p.adjust.methods)
    return(list("error"=TRUE, "message"="Invalid p-value correction method name. Refer to stats::p.adjust.methods."))
  if ( ( ! is.numeric(p_thresh)) || p_thresh > 1 || p_thresh < 0)
    return(list("error"=TRUE, "message"="Invalid p-value threshold."))
  if ( ! varname %in% names(sleuth_data$sample_to_covariates))
    return(list("error"=TRUE, "message"="The specified covariate name does not exist."))
  if ( any( ! c(ref_name, comp_name) %in% sleuth_data$sample_to_covariates[[varname]] ))
    return(list("error"=TRUE, "message"="One or both of the specified conditions do not exist."))
  if ( ! is.logical(verbose))
    return(list("error"=TRUE, "message"="verbose must be a logical value."))
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
  progress_steps <- data.frame(c(10,20,30,40,50,60,70,80,90,100),
                               c("Mapped genes to target targets",
                                 "Identified genes with multiple transcripts",
                                 "Grouped samples by condition",
                                 "Extracted counts from bootstraps",
                                 "Removed all-zero count cases",
                                 "Allocated output structure",
                                 "Identified genes and transcripts represented in the data",
                                 "Calculated counts statistics",
                                 "Calculated transcript proportions",
                                 "Calculated p-values"),
                               stringsAsFactors = FALSE)
  progress <- TxtProgressUpdate(steps=progress_steps, on=on)
  return(progress)
}

