#' Compute a logical vector filter marking multi-child parent IDs in a data frame
#'
#' @param ids a data frame with at least two variables, \code{target_id} & \code{parent_id}
#'
#'
#'
#' @export

mark_branched_IDs <- function(ids) {
  # Look for duplicate target_ids.
  branch_count <- as.data.frame(table(ids$target_id))  # basically a hash table with Var1 contining the target_ids exactly once.
  # Use the non redundant list of target_ids to find out how many times the parent_ids are used.
  root_ids <- sapply(branch_count$Var1, function(x) strsplit(as.character(x), "\\.")[[1]][1])
  root_id_counts <- table(root_ids)
  # Add boolean column
  ids$has_multi <- root_id_counts[ids$parent_id] > 1
  return ids
}
