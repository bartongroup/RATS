# constant declarations - do not export in NAMESPACE
TARGET_ID = "target_id"     # name of transcript id column in transcripts object
PARENT_ID = "parent_id"     # name of parent id column in transcripts object
BS_TARGET_ID = "target_id"  # name of transcript id column in sleuth bootstrap tables

#================================================================================
#================================================================================
#' Calculate differential transcript usage.
#'
#' @param sleuth_data A sleuth object
#' @param transcripts A dataframe matching the transcript IDs (\code{target_id}) to their corresponding gene IDs (\code{parent_id}).
#' @param ref_name The sleuth name of the reference condition.
#' @param comp_name The sleuth name of the condition to compare.
#' @param varname The sleuth name of the covariate to which the two conditions belong.
#' @param counts_col The sleuth column to use for the calculation (est_counts or tpm), default est_counts
#' @return List of data frames, with gene-level and transcript-level information.
#'
#' @export
calculate_DTU <- function(sleuth_data, transcripts, ref_name, comp_name, varname="condition", counts_col="est_counts") {

  # set up progress bar
  progress_steps <- c(10,20,30,40,50,60,70,80,90,100)
  pb = init_progress(progress_steps)
  pb = update_progress(pb, progress_steps)

  # Look-up from parent_id to target_id (slow).
  targets_by_parent <- parent_to_targets(transcripts)

  pb = update_progress(pb, progress_steps)

  # Identify genes with a single transcript. Order by gene ID and transcript ID.
  tx_filter <- mark_sibling_targets3(transcripts, targets_by_parent)

  pb = update_progress(pb, progress_steps)

  # Reverse look-up from relicates to covariates.
  samples_by_condition <- group_samples(sleuth_data$sample_to_covariates)[[varname]]

  pb = update_progress(pb, progress_steps)

  # build list of dataframes, one for each condition
  # each dataframe contains filtered and correctly ordered counts from all the bootstraps for the condition
  count_data <- lapply(samples_by_condition, function(condition) make_filtered_bootstraps(sleuth_data, condition, tx_filter, counts_col))

  pb = update_progress(pb, progress_steps)

  # remove entries which are entirely 0 across all conditions
  subsets <-  lapply(count_data, function(condition) apply(condition, 1, function(row) !all(row == 0 )))
  count_data <- lapply(count_data, function(condition) condition[Reduce("&", subsets),])

  pb = update_progress(pb, progress_steps)

  # Which IDs am I actually working with after the filters?
  actual_targets <- rownames(count_data[[ref_name]])
  actual_parents <- levels(as.factor(tx_filter[[PARENT_ID]][match(actual_targets, tx_filter[[TARGET_ID]])]))

  pb = update_progress(pb, progress_steps)

  # Pre-allocate output structure.
  results <- list("Comparison"=c("variable_name"=varname, "reference"=ref_name, "compared"=comp_name),
                  "Genes"=data.frame("considered"=c(FALSE),
                                     "parent_id"=levels(as.factor(tx_filter[[PARENT_ID]])),
                                     "dtu"=NA, "p_value"=NA_real_),
                  "Transcripts"=data.frame("considered"=c(FALSE), "target_id"=tx_filter[[TARGET_ID]], "parent_id"=tx_filter[[PARENT_ID]],
                                           "ref_proportion"=NA_real_, "comp_proportion"=NA_real_,
                                           "ref_sum"=NA_real_, "comp_sum"=NA_real_,
                                           "ref_mean"=NA_real_, "ref_variance"=NA_real_,
                                           "comp_mean"=NA_real_, "comp_variance"=NA_real_))
  rownames(results$Genes) <- results$Genes$parent_id
  rownames(results$Transcripts) <- results$Transcripts$target_id

  pb = update_progress(pb, progress_steps)

  # Statistics per transcript across all bootstraps per condition.
  results$Transcripts[actual_targets, c("ref_sum", "ref_mean", "ref_variance", "comp_sum", "comp_mean", "comp_variance")] <-
    c(rowSums(count_data[[ref_name]]), rowMeans(count_data[[ref_name]]), matrixStats::rowVars(as.matrix(count_data[[ref_name]])),
      rowSums(count_data[[comp_name]]), rowMeans(count_data[[comp_name]]), matrixStats::rowVars(as.matrix(count_data[[comp_name]])))

  pb = update_progress(pb, progress_steps)

  # Proportions = sum of tx / sum(sums of all related txs).
  for (target in results$Transcripts$target_id) {
    results$Transcripts[target, c("ref_proportion", "comp_proportion")] <-
      c(results$Transcripts[target, "ref_sum"] / sum(results$Transcripts[ targets_by_parent[[results$Transcripts[target, "parent_id"]]], "ref_sum"]),
        results$Transcripts[target, "comp_sum"] / sum(results$Transcripts[ targets_by_parent[[results$Transcripts[target, "parent_id"]]], "comp_sum"]))
  }

  pb = update_progress(pb, progress_steps)

  # P values.
  for (p in actual_parents) {
    targets <- targets_by_parent[[p]]
    gt <- g.test(results$Transcripts[targets, "comp_sum"], p=results$Transcripts[targets, "ref_proportion"])
    results$Genes[p, "p_value"] <- gt$p.value
    results$Genes[p, c("considered", "dtu")] <- c(TRUE, gt$p.value < 0.05)
    results$Transcripts[targets, "considered"] <- TRUE
  }

  # TODO : Multiple testing correction.

  pb = update_progress(pb, progress_steps)
  close(pb)

  return(results)
}
#================================================================================
#================================================================================

#================================================================================
#' Compute a logical vector marking as FALSE the single-target parents in a data frame.
#'
#' @param ids a data frame with at least two variables, \code{target_id} & \code{parent_id}.
#' @param p2t a list of vectors, listing the \code{target_id}s per \code{parent_id}.
#' @return data.frame An updated version of the input ids.
#'
mark_sibling_targets3 <- function(ids, p2t) {
  rownames(ids) <- ids[[TARGET_ID]]
  for (target in ids[[TARGET_ID]]) {
    ids[target, "has_siblings"] <- length(p2t[[ ids[target, PARENT_ID] ]]) > 1
  }

  return(ids[order(ids[[PARENT_ID]], ids[[TARGET_ID]]), ])
}

#--------------------------------------------------------------------------------
#' Group targets by parent.
#'
#' @param ids A dataframe with \code{target_id} and \code{parent_id}.
#' @return A list of vectors: one vector for each parent, containing the repsective targets.
#'
parent_to_targets <- function(ids) {
#  p <- "AT1G01020"
#  ids <- transcripts

  parents <- levels(as.factor(ids[[PARENT_ID]]))
  p2t <- lapply(parents, function(p) ids[[TARGET_ID]][ids[[PARENT_ID]]==p])
  names(p2t) <- parents

  return(p2t)
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
#' @return A dataframe containing the counts from all bootstraps of all the samples for the condition.
#'
make_filtered_bootstraps <- function(sleuth_data, condition, tx_filter, counts_col) {

  # make a list of dataframes, one df for each condition, containing the counts from its bootstraps
  count_data <- as.data.frame(lapply (condition, function(sample)
      sapply(sleuth_data$kal[[sample]]$bootstrap, function(e) filter_and_match(e, tx_filter, counts_col) )))

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
#' @return The column from the bootstrap table
#'
filter_and_match <- function(bootstrap, tx_filter, counts_col)
{
  # create map from bootstrap to filter(i.e. main annotation) target ids
  b_to_f_rows <- match(tx_filter[[TARGET_ID]], bootstrap[[BS_TARGET_ID]])

  # map the bootstrap to the filter target ids and then apply the filter
  result <- (bootstrap[b_to_f_rows,][,counts_col])[tx_filter$has_siblings]

  return(result)
}
