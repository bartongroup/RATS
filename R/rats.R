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
#' Row number corresponds to smaple number. Does not assume proximity of same-factor samples.
#' Assumes that a factor's value is not also another factor's name.
#' Returns list of vectors. Dataframe inappropriate as the vectors may differ in length.
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

  max_b <- 100

  # make a list of dataframes, one df for each sample, containing the counts from its bootstraps
  samples <- 1:nrow(sleuth_data$sample_to_covariates)
  count_data <- lapply(samples, function(x) as.data.frame(sapply(sleuth_data$kal[[x]]$bootstrap, function(e) e[, counts_col])))

  # add in transcript ids
  tx_ids <- sleuth_data$kal[[1]]$bootstrap[[1]]["target_id"] #assume target ids in same order in all dataframes
  count_data <- Map(cbind, count_data, target_id = tx_ids)

  # calculate mean/variance across samples of same condition
  # count_data[[1]]$testMean <- rowMeans(count_data[[1]][,1:max_b], na.rm=TRUE)
  # count_data[[1]]$testVar <- apply(count_data[[1]][,1:max_b],1,var)

  # merge with gene ids no not yet - use merge(data, transcripts)
  #count_stats <- aggregate(count_data[[1]][, 1:100], list(count_data[[1]]$target_id), mean)

}
