# constant declarations - do not export in NAMESPACE
CONDITION_COL = "condition" # name of condition column in sleuth sample_to_covariates object
TARGET_ID = "target_id"     # name of transcript id column in transcripts object
BS_TARGET_ID = "target_id"  # name of transcript id column in sleuth bootstrap tables

#' Calculate differential transcript usage
#' @param sleuth_data A sleuth object
#' @param transcripts A dataframe listing the transcripts to process, and their parent genes,
#' with at least column \code{target_id} & \code{parent_id}
#' @param counts_col The sleuth column to use for the calculation (est_counts or tpm), default est_counts
#' @return (TODO) dataframe with a row for each transcript, indicating if it is DTU and giving p-value etc.
calculate_DTU <- function(sleuth_data, transcripts, counts_col="est_counts") {

  # get full set of target_id filters
  filter <- mark_sibling_targets2(transcripts)

  # get which samples correspond to which condition
  samples_by_condition <- group_samples(sleuth_data$sample_to_covariates)[[CONDITION_COL]]

  # build list of dataframes, one for each condition
  # each dataframe contains filtered and correctly ordered counts from all the bootstraps for the condition
  count_data <- make_filtered_bootstraps(sleuth_data, samples_by_condition, filter, counts_col)

  # calculate the relative proportion of expression of each transcript
  proportions <- calculate_tx_proportions(count_data)

  # TODO DTU calc and error of proportion test
  return(proportions) # for now just return proportions
}

#' Compute a logical vector filter marking single-target parents in a data frame.
#'
#' @param df a data frame with at least two variables, \code{target_id} & \code{parent_id}
#'
#' @export
mark_sibling_targets <- function(df){
  parent_group <- dplyr::group_by(df, parent_id)
  unique_parent_count<-dplyr::summarise(parent_group,n())
  singles_filter <- unique_parent_count[2]>1
  multi_parent_ids=unique_parent_count[singles_filter,]
  df_filter <- apply(df, 1, function(r,s) any(r["parent_id"] %in% s ), s=multi_parent_ids[[1]])
  if (sum(df_filter)!=sum(multi_parent_ids[2])) {
    warning("Something went wrong with the identification. The number targets of identified
            multiple target parent ids does not match the number in the final filter!")
  }
  df$has_siblings <- df_filter
  return(df)
}


#' Compute a logical vector marking as FALSE the single-target parents in a data frame.
#'
#' @param ids a data frame with at least two variables, \code{target_id} & \code{parent_id}
#' @param duptx a boolean switch indicating whether to account for duplicate target_ids (default FALSE)
#'
#' This is a second implementation of this routine. Faster, but less general.
#' Target IDs are assumed to have the format \code{parentID.childextension} .
#'
#' @export
mark_sibling_targets2 <- function(ids, duptx=FALSE) {
  if (duptx) {
    # EITHER count for duplicate target_ids. Then I can use the hash keys as a non-redundant list of target_ids.
    branch_count <- as.data.frame(table(ids$target_id))
    root_ids <- sapply(branch_count$Var1, function(x) strsplit(as.character(x), "\\.")[[1]][1])
  } else {
    # OR assume the target_ids are already non-redundant
    root_ids <- sapply(ids$target_id, function(x) strsplit(as.character(x), "\\.")[[1]][1])
  }
  # Count parents and mark as TRUE targets that share their parent with other targets.
  root_id_counts <- table(root_ids)
  ids$has_siblings <- as.vector(root_id_counts[ids$parent_id] > 1)
  return(ids)
}


#' Group sample numbers by factor.
#'
#' @param covariates a dataframe with different factor variables.
#'
#' Row number corresponds to sample number. Does not assume proximity of same-factor samples.
#' Assumes that a factor's value is not also another factor's name.
#' Returns nested lists of vectors. Dataframe inappropriate as the vectors may differ in length.
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

#' For each condition in the sleuth object, construct a dataframe containing counts from each bootstrap,
#' filtered according to filter, and ordered according to target_ids
#' @param sleuth_data A sleuth object
#' @param samples_by_condition The samples which correspond to each condition
#' @param filter A filter for the bootstraps, with corresponding target_ids
#' @param counts_col The sleuth column to use for the calculation
#' @return List of dataframes, where each dataframe contains the counts from all bootstraps for one condition
make_filtered_bootstraps <- function(sleuth_data, samples_by_condition, filter, counts_col) {

  # make a list of dataframes, one df for each condition, containing the counts from its bootstraps
  count_data <- lapply(samples_by_condition, function(condition)
    as.data.frame(lapply (condition, function(sample)
      sapply(sleuth_data$kal[[sample]]$bootstrap, function(e) filter_and_match(e, filter, counts_col) ))))

  # now set the filtered target ids as rownames on each condition - previous call returns target ids in this order
  count_data <- lapply(count_data, function(condition) {rownames(condition) <- filter[[TARGET_ID]][filter$has_siblings]; condition})

  # remove any row containing NAs: some bootstrap did not have an entry for the transcript for the row
  # TODO this may not be required behaviour
  count_data <- lapply(count_data, function(condition) condition[complete.cases(condition),])

  return(count_data)
}

#' Calculate the proportion of counts which are assigned to each transcript in a gene
#'
#' @param count_data List of dataframes, where each dataframe contains the counts from all bootstraps for one condition
#' @return Proportion of counts which are assigned to each transcript in a gene (but currently just mean and variance per transcript)
calculate_tx_proportions <- function(count_data) {

  # TODO
  # try to do mean/var calculations in place rather than creating a new count_data dataframe
  # calculate the proportions

  # calculate mean and variance across all samples of the same condition
  metrics <- lapply(count_data, function(condition) calculate_stats(condition))
  return(metrics)
}

#' Extract a filtered column from a bootstrap table
#'
#' @param bootstrap A bootstrap dataframe
#' @param filter A boolean vector for filtering the bootstrap, with matching transcript ids
#' @param counts_col The column in the bootstrap table to extract
#' @return The column from the bootstrap table
filter_and_match <- function(bootstrap, filter, counts_col)
{
  # create map from bootstrap to filter(i.e. main annotation) target ids
  b_to_f_rows <- match(filter[[TARGET_ID]], bootstrap[[BS_TARGET_ID]])

  # map the bootstrap to the filter target ids and then apply the filter
  result <- (bootstrap[b_to_f_rows,][,counts_col])[filter$has_siblings]

  return(result)
}

#' Calculate statistics across the columns of a dataframe
#'
#' @param df dataframe
#' @return A dataframe containing the calculated statistics, one in each column
calculate_stats <- function(df)
{
  metric <- as.data.frame(t(apply(df, 1, mean_and_var)))
  return(metric)
}

#' Calculate mean and variance across the columns of a row
#'
#' @param row the row
#' @return A vector containing the mean and variance
mean_and_var <- function(row)
{
  return(c(mean = mean(row), variance = var(row)))
}

#' Calculate mean and variance.
#'
#' @param kal A list of kallisto objects.
#' @param targets A logic vector indicating which rows in the bootstraps to calculate for.
#' @param counttype "est_counts" or "tmp".
#'
#' Returns a dataframe with mean and variance per target_id
#'
#' @export
do_count_stats <- function(kal, targets=c(TRUE), counttype="tpm") {
#  targets <- filtered_ids
#  reps <- covariates_to_samples[["condition"]][["Col"]]
#  counttype <- "est_counts"

  # Initialize dataframe
  metrics <- data.frame("target_id"=kal[[1]]$bootstrap[[1]]$target_id[targets], "mean"=NA_real_, "variance"=NA_real_)

  # Extact counts across bootstraps across replicates.
  counts <- as.data.frame(lapply(kal, function(k) sapply(k$bootstrap, function(b) b[targets, counttype])))

  # Calculate stats
  metrics["mean"] <- rowMeans(counts)
  metrics["variance"] <- rowVars(as.matrix(counts))

  return(metrics)
}


