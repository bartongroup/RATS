#' Compute a logical vector filter marking single-target parents in a data frame.
#'
#' @param df a data frame with at least two variables, \code{target_id} & \code{parent_id}
#'
#' @export
mark_singleparent_ids <- function(df){
  parent_group <- group_by(df, parent_id)
  unique_parent_count<-summarise(parent_group,n())
  singles_filter <- unique_parent_count[2]>1
  multi_parent_ids=unique_parent_count[singles_filter,]
  df_filter <- apply(df, 1, function(r,s) any(r["parent_id"] %in% s ), s=multi_parent_ids[[1]])
  if (sum(df_filter)!=sum(multi_parent_ids[2])) {
    warning("Something went wrong with the identification. The number targets of identified
            multiple target parent ids does not match the number in the final filter!")
  }
  df$has_multi <- df_filter
  return(df)
}

#' Compute a logical vector filter marking multi-child parent IDs in a data frame.
#'
#' @param ids a data frame with at least two variables, \code{target_id} & \code{parent_id}
#'
#' This is a second implementation of this routine. Faster, but less general.
#'
#' @export
mark_singleparent_ids2 <- function(ids) {
  # Look for duplicate target_ids.
  # branch_count is basically a hash table with Var1 contining the target_ids exactly once.
  branch_count <- as.data.frame(table(ids$target_id))
  # Use the non redundant list of target_ids to find out how many times the parent_ids are used.
  root_ids <- sapply(branch_count$Var1, function(x) strsplit(as.character(x), "\\.")[[1]][1])
  root_id_counts <- table(root_ids)
  # Add boolean column
  ids$has_multi <- root_id_counts[ids$parent_id] > 1
  return(ids)
}
