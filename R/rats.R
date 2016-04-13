#' Compute a logical vector filter marking single-target parents in a data frame
#'
#' @param df a data frame with at least two variables, \code{target_id} & \code{parent_id}
#'
#' @export

markSingles <- function(df){
  parent_group <- group_by(df, parent_id)
  unique_parent_count<-summarise(parent_group,n())
  singles_filter <- unique_parent_count[2]>1
  multi_parent_ids=unique_parent_count[singles_filter,]
  df_filter <- apply(df, 1, function(r,s) any(r["parent_id"] %in% s ), s=multi_parent_ids[[1]])
  if (sum(df_filter)!=sum(multi_parent_ids[2])) {
    warning("Something went wrong with the identification. The number targets of identified
            multiple target parent ids does not match the number in the final filter!")
  }
  df$singles_filter <- df_filter
  return(df)
}
