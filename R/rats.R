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
