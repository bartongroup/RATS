#================================================================================
#' Compute a logical vector filter marking single-target parents in a data frame.
#'
#' @param df a data frame with at least two variables, \code{target_id} & \code{parent_id}
#'
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
  rownames(df) <- df$target_id
  return(df[order(df$parent_id, df$target_id),])
  }

#--------------------------------------------------------------------------------
#' Compute a logical vector marking as FALSE the single-target parents in a data frame.
#' (alternative implementation)
#'
#' @param ids a data frame with at least two variables, \code{target_id} & \code{parent_id}
#' @param duptx a boolean switch indicating whether to account for duplicate target_ids (default FALSE)
#' @return data.frame An updated version of the input ids.
#'
#' Target IDs are assumed to have the format \code{parentID.childextension}.
#'
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
  rownames(ids) <- ids$target_id
  return(ids[order(ids$parent_id, ids$target_id),])
}

#--------------------------------------------------------------------------------
#' Extract filtered counts, assuming the bootstraps and replicates are identical in their
#' order, number and presence of transcripts.
#'
#' @param kal A list of kallisto objects.
#' @param tx_filter A logic vector to select bootstap rows (the same ones across all bootstraps).
#' @param counttype "tpm" or "est_counts".
#'
filter_and_match2 <- function(kal, tx_filter=c(TRUE), counttype="tpm") {
  counts <- as.data.frame(lapply(kal, function(k) sapply(k$bootstrap, function(b) b[tx_filter, counttype])))
  return(counts)
}

#================================================================================
#' Calculate statistics across the columns of a dataframe
#'
#' @param df dataframe
#' @return A dataframe containing the calculated statistics, one in each column
#'
calculate_stats <- function(df)
{
  metric <- as.data.frame(t(apply(df, 1, mean_and_var)))
  return(metric)
}

#--------------------------------------------------------------------------------
#' Calculate mean and variance across the columns of a row
#'
#' @param row the row
#' @return A vector containing the mean and variance
#'
mean_and_var <- function(row)
{
  return(c(sum = sum(row), mean = mean(row), variance = var(row)))
}

#--------------------------------------------------------------------------------
#' Calculate row-wise statistics in a dataframe.
#'
#' @param df A dataframe.
#' @return A dataframe containing the sum, mean, variance per row.
#'
#'It expects and preserves rownames.
#'
calculate_stats2 <- function(df) {
  # Prepare empty dataframe.
  txid = rownames(df)
  metrics <- data.frame("sum"=rep(NA_real_, length(txid)), "mean"=NA_real_, "variance"=NA_real_)
  rownames(metrics) <- txid
  metrics["sum"] <- rowSums(df)
  metrics["mean"] <- rowMeans(df)
  metrics["variance"] <- matrixStats::rowVars(as.matrix(df))
  return(metrics)
}


#================================================================================
#' Calculate the proportion of counts which are assigned to each transcript in a gene
#'
#' @param count_data A dataframe containing the counts from all bootstraps for one condition
#' @param transcripts A dataframe listing the transcripts to process, and their parent genes
#' @return Proportion of counts which are assigned to each transcript in a gene (but currently just mean and variance per transcript)
#'
#' @import data.table
calculate_tx_proportions <- function(count_data, transcripts) {

  # calculate mean and variance across all samples of the same condition
  metrics <- calculate_stats(count_data)

  # add gene ids to dataframe
  metrics$parent_id <- transcripts[[PARENT_ID]][match(rownames(metrics), transcripts[[TARGET_ID]])];

  # calculate reads proportion per transcript by gene
  dt = data.table(target_id=rownames(metrics), metrics)
  dt[,proportion := mean/sum(mean), by=parent_id] # import data.table needed to make := work

  return(dt)
}
