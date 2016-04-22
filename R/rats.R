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
  # try to do mean/var calculations in place rather than creating a new count_data dataframe
  # calculate the proportions

  # get full set of target_id filters
  filter <- mark_sibling_targets2(transcripts)

  samples_by_condition <- group_samples(sleuth_data$sample_to_covariates)[[CONDITION_COL]]

  # reduce filter to match entries in bootstraps (assumes all bootstraps have same entries)
  f_to_b_rows <- match(sleuth_data$kal[[1]]$bootstrap[[1]]$target_id, filter$target_id)
  filter <- filter[f_to_b_rows,]$has_siblings
  target_ids <- transcripts [f_to_b_rows,]$target_id

  # make a list of dataframes, one df for each condition, containing the counts from its bootstraps

  count_data <- lapply(samples_by_condition, function(condition)
    as.data.frame(lapply (condition, function(sample)
    sapply(sleuth_data$kal[[sample]]$bootstrap, function(e) e[filter,][, counts_col]))))

  # calculate mean and variance across all samples of the same condition
  metrics <- lapply(count_data, function(condition) calculate_stats(condition))
  metrics <- lapply(metrics, function(m) cbind.data.frame(target_id = target_ids[filter],m))
  return(metrics)
}

calculate_stats <- function(x)
{
  return (t(apply(x, 1, mean_and_var)))
}

mean_and_var <- function(x)
{
  return(c(mean = mean(x), variance = var(x)))
}

bs_filter <- function(bootstrap, global_filter)
{
  f_to_b_rows <- match(bootstrap$target_id, global_filter$target_id)
  filter <- filter[f_to_b_rows,]$has_siblings
  return(filter)
}


#' Calculate mean and variance.
#'
#' @param targets A vector of target_id's.
#' @param sl A sleuth object.
#' @param reps A vector indicating which replicates to use.
#' @param counttype "est_counts" or "tmp".
#'
#' Returns a dataframe with mean and variance per target_id
#'
#' @export
do_count_stats <- function(targets, sl, reps, counttype="tpm") {
#  targets <- filtered_ids
#  reps <- covariates_to_samples[["condition"]][["Col"]]
#  counttype <- "est_counts"

  # Initialize dataframe
  metrics <- data.frame("target_id"=targets, "mean"=NA, "variance"=NA)
  rownames(metrics) <- metrics$target_id

  # I think as.dataframe() will copy the values into a new location of RAM in the requested format, so I might as well name it and re-use it.
  counts <- as.data.frame(lapply(sl$kal[reps], function(k) sapply(k$bootstrap, function(b) b[metrics$target_id, counttype])))
  rownames(counts) <- metrics$target_id   # counts was created based on metrics$target_id so I am able to do this with confidence that it is indeed so.

  # Calculate means
  metrics["mean"] <- rowMeans(counts)
#  metrics["mean"] <- rowMeans( as.data.frame(lapply(sl$kal[reps], function(k) sapply(k$bootstrap, function(b) b[metrics$target_id, counttype])))  )

  # Calculate variances (from package matrixStats)
  metrics["variance"] <- rowVars(as.matrix(counts))

  return(metrics)
}


#' Main
#'
#' @param sl A sleuth object.
#' @param ids A dataframe with \code{target_id}, corresponding \code{parent_id}, and other corresponding  id's that are ignored.
#'
#' Does not produce any visible output. For workflow demonstration purpose.
#'
#' @export
main <- function(sl=At_250genes_50singletx_200multitx_513totaltx, ids=At_tair10_t2g) {
#  devtools::load_all()
#  sl <- At_250genes_50singletx_200multitx_513totaltx   # sleuth object
#  ids <- At_tair10_t2g

  # Create reverse look-up of covariates to samples.
  covariates_to_samples <- group_samples(sl$sample_to_covariates)

  # Mark target_id's that have alternative transcript siblings.
  ids <- mark_sibling_targets2(At_tair10_t2g)

  # Assign index to ids
  rownames(ids) <- ids$target_id      # assign transcript_id as row index, without dropping

  # Assign index to bootstraps
  for (k in 1:length(sl$kal)) {
    for (b in 1:length(sl$kal[[k]]$bootstrap)) {
      rownames(sl$kal[[k]]$bootstrap[[b]]) <- sl$kal[[k]]$bootstrap[[b]]$target_id
    }
  }

  # Reduce the ids to the ones actually present in the sleuth object and select those that have alternative splicing siblings.
  exist_ids <- ids[ ids$target_id[ ids$target_id %in% sl$kal[[1]]$bootstrap[[1]]$target_id ], ]   # still a dataframe
  filtered_ids <- exist_ids$target_id[exist_ids$has_siblings]                                     # vector

  # Calculate means and variances per transcript (only for those that have alternative transcript siblings).
  metrics <- do_count_stats(filtered_ids,
                            sl,
                            reps=covariates_to_samples[["condition"]][["Col"]],
                            counttype="est_counts")
}





