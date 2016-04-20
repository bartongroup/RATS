# constant declaration - do not export in NAMESPACE
CONDITION_COL = "condition" # name of condition column in sleuth sample_to_covariates object

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
  ids$has_siblings <- root_id_counts[ids$parent_id] > 1
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

#' Calculate the proportion of counts which are assigned to each transcript in a gene
#'
#' @param sleuth_data a sleuth object
#' @param transcripts a dataframe listing the transcripts to process, and their parent genes,
#' with at least column \code{target_id} & \code{parent_id}
#' @param counts_col the sleuth column to use for the calculation (est_counts or tpm), default est_counts
#'
#' @export
calculate_tx_proportions <- function(sleuth_data, transcripts, counts_col="est_counts") {

  # TODO
  # add transcript ids as index
  # filter by transcript ids
  # try to do mean/var calculations in place rather than creating a new count_data dataframe
  # calculate the proportions

  # get full set of target_id filters
  filter <- mark_sibling_targets(At_tair10_t2g)

  # reduce filter to match entries in bootstraps (assumes all bootstraps have same entries)
  filter <- filter[ filter$target_id %in% sleuth_data$kal[[1]]$bootstrap[[1]]$target_id, "has_siblings" ]

  # make a list of dataframes, one df for each condition, containing the counts from its bootstraps
  samples_by_condition <- group_samples(sleuth_data$sample_to_covariates)[[CONDITION_COL]]
  count_data <- lapply(samples_by_condition, function(condition)
    as.data.frame(lapply (condition, function(sample) sapply(sleuth_data$kal[[sample]]$bootstrap, function(e) e[filter,][, counts_col]))))

  # calculate mean and variance across all samples of the same condition
  means <- lapply(count_data, function(condition) apply(condition, 1, mean))
  vars <- lapply(count_data, function(condition) apply(condition, 1, var))

}

