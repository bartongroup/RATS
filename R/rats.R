#================================================================================
#' Calculate differential transcript usage.
#'
#' @param sleuth_data A sleuth object
#' @param transcripts A dataframe matching the transcript IDs to their corresponding gene IDs.
#' @param ref_name The sleuth name of the reference condition.
#' @param comp_name The sleuth name of the condition to compare.
#' @param varname The sleuth name of the covariate to which the two conditions belong, default \code{"condition"}.
#' @param counts_col The sleuth column to use for the calculation (est_counts or tpm), default \code{"est_counts"}.
#' @param correction p-value correction to apply, as defined in \code{stats::p.adjust.methods}, default \code{"BH"}.
#' @param TARGET_ID The name of the transcript id column in transcripts object, default \code{"target_id"}
#' @param PARENT_ID The name of the parent id column in transcripts object, default \code{"parent_id"}.
#' @param BS_TARGET_ID The name of the transcript id column in sleuth bootstrap tables, default \code{"target_id"}.
#' @param verbose Whether to update progress updates, default \code{FALSE}.
#' @return List of data frames, with gene-level and transcript-level information.
#'
#' @export
calculate_DTU <- function(sleuth_data, transcripts, ref_name, comp_name,
                          varname="condition", counts_col="est_counts", correction="BH",
                          TARGET_ID="target_id", PARENT_ID="parent_id", BS_TARGET_ID="target_id",
                          verbose=FALSE)
{
  # Input checks.
  paramcheck <- parameters_good(sleuth_data, transcripts, ref_name, comp_name, varname, counts_col,
                  correction, TARGET_ID, PARENT_ID, BS_TARGET_ID, verbose)
  if (paramcheck$error) stop(paramcheck$message)

  # Set up progress bar
  progress_steps <- data.frame(c(10,20,30,40,50,60,70,80,90,95,100),
                               c("Built parent ids to target map","Built multiple transcript filter",
                                 "Grouped samples by condition","Extracted counts from bootstraps","Removed 0 count cases",
                                 "Filtered parent ids", "Allocated output structure", "Calculated statistics",
                                 "Calculated proportions", "Calculated p-values", "Finished!"),
                               stringsAsFactors = FALSE)
  progress <- TxtProgressUpdate(steps=progress_steps, on=verbose)

  # Look-up from parent_id to target_id
  targets_by_parent <- split(as.matrix(transcripts[TARGET_ID]), transcripts[[PARENT_ID]])

  progress <- update(progress)

  # Identify genes with a single transcript. Order by gene ID and transcript ID.
  tx_filter <- mark_sibling_targets3(transcripts, targets_by_parent, TARGET_ID, PARENT_ID)

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
  count_data <- lapply(count_data, function(condition) condition[Reduce("&", nonzero),])

  progress <- update(progress)

  # Pre-allocate output structure.
  results <- list("Comparison"=c("variable_name"=varname, "reference"=ref_name, "compared"=comp_name),
                  "Genes"=data.frame("considered"=c(FALSE),
                                     "parent_id"=levels(as.factor(tx_filter[[PARENT_ID]])),
                                     "dtu"=NA, "p_value"=NA_real_, "corrected_p"=NA_real_,
                                     "num_known_transc"=NA_integer_, "num_applicable_transc"=NA_integer_),
                  "Transcripts"=data.frame("considered"=c(FALSE), "target_id"=tx_filter[[TARGET_ID]], "parent_id"=tx_filter[[PARENT_ID]],
                                           "ref_proportion"=NA_real_, "comp_proportion"=NA_real_,
                                           "ref_sum"=NA_real_, "comp_sum"=NA_real_,
                                           "ref_mean"=NA_real_, "ref_variance"=NA_real_,
                                           "comp_mean"=NA_real_, "comp_variance"=NA_real_))
  rownames(results$Genes) <- results$Genes$parent_id
  rownames(results$Transcripts) <- results$Transcripts$target_id

  progress <- update(progress)

  # Which IDs am I actually working with after the filters?
  actual_targets <- rownames(count_data[[ref_name]])
  actual_parents <- levels(as.factor(tx_filter[[PARENT_ID]][match(actual_targets, tx_filter[[TARGET_ID]])]))

  progress <- update(progress)

  # Statistics per transcript across all bootstraps per condition, for filtered targets only.
  results$Transcripts[actual_targets, c("ref_sum", "ref_mean", "ref_variance", "comp_sum", "comp_mean", "comp_variance")] <-
    c(rowSums(count_data[[ref_name]]), rowMeans(count_data[[ref_name]]), matrixStats::rowVars(as.matrix(count_data[[ref_name]])),
      rowSums(count_data[[comp_name]]), rowMeans(count_data[[comp_name]]), matrixStats::rowVars(as.matrix(count_data[[comp_name]])))

  progress <- update(progress)

  # Proportions = sum of tx / sum(sums of all related txs), for filtered targets only.
  # The need to access multiple rows and different objects, means this cannot be done in a simple vector-style assignment. Hence the loop.
  for (target in results$Transcripts$target_id) {
    results$Transcripts[target, c("ref_proportion", "comp_proportion")] <-
      c(results$Transcripts[target, "ref_sum"] / sum(results$Transcripts[ targets_by_parent[[results$Transcripts[target, "parent_id"]]], "ref_sum"]),
        results$Transcripts[target, "comp_sum"] / sum(results$Transcripts[ targets_by_parent[[results$Transcripts[target, "parent_id"]]], "comp_sum"]))
  }

  progress <- update(progress)

  # P values, only for parents and targets that survived filtering.
  for (p in actual_parents) {
    targets <- targets_by_parent[[p]]  # all the annotated transcripts for the gene.
    results$Genes[p,"num_known_transc"] <- length(targets)
    targets <- targets[targets %in% actual_targets]  # the transcripts for which we have non-zero counts.
    results$Genes[p,"num_applicable_transc"] <- length(targets)

    if (length(targets) < 2) next  # Can't do DTU with less than two transcripts.

    gt <- g.test(results$Transcripts[targets, "comp_sum"], p=results$Transcripts[targets, "ref_proportion"])
    results$Genes[p, "p_value"] <- gt$p.value
    results$Genes[p, c("considered", "dtu")] <- c(TRUE, gt$p.value < 0.05)
    results$Transcripts[targets, "considered"] <- TRUE
  }

  # Multiple testing correction.
  results$Genes[["corrected_p"]] <- p.adjust(results$Genes[["p_value"]], method=correction)

  progress <- update(progress)

  # Tidy up some information gaps.
  # Proportion for singletons is 1.
  results$Transcripts[ !tx_filter$has_siblings, "ref_proportion"] <- 1
  results$Transcripts[ !tx_filter$has_siblings, "comp_proportion"] <- 1
  # Parents of singletons have 1 known transcript and 0 applicable transcripts.
  for (p in results$Genes$parent_id) {
    if (length(targets_by_parent[[p]]) == 1 ) {
      results$Genes[p, c("num_known_transc", "num_applicable_transc")] <- c(1,0)
    }
  }

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
mark_sibling_targets3 <- function(ids, p2t, TARGET_ID, PARENT_ID) {
  rownames(ids) <- ids[[TARGET_ID]]
  for (target in ids[[TARGET_ID]]) {
    ids[target, "has_siblings"] <- length(p2t[[ ids[target, PARENT_ID] ]]) > 1
  }

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
#' For each condition in the sleuth object, construct a dataframe containing counts from each bootstrap,
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
  count_data <- as.data.frame(lapply (condition, function(sample)
      sapply(sleuth_data$kal[[sample]]$bootstrap, function(e) filter_and_match(e, tx_filter, counts_col, TARGET_ID, BS_TARGET_ID) )))

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
  result <- (bootstrap[b_to_f_rows,][,counts_col])[tx_filter$has_siblings]

  return(result)
}

#================================================================================
#' Check input parameters.
#'
#' @return Logical value
#'
parameters_good <- function(sleuth_data, transcripts, ref_name, comp_name, varname, counts_col,
                            correction, TARGET_ID, PARENT_ID, BS_TARGET_ID, verbose) {
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
  if ( ! varname %in% names(sleuth_data$sample_to_covariates))
    return(list("error"=TRUE, "message"="The specified covariate name does not exist."))
  if ( any( ! c(ref_name, comp_name) %in% sleuth_data$sample_to_covariates[[varname]] ))
    return(list("error"=TRUE, "message"="One or both of the specified conditions do not exist."))
  if ( ! is.logical(verbose))
    return(list("error"=TRUE, "message"="verbose must be a logical value."))
  return(list("error"=FALSE, "message"="All good!"))
}
