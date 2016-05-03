#================================================================================
#' Initialise progress reporting
#' @param progress_steps Vector of values to use as progress bar updates
#' (should match number of calls being made to update_progress)
#'
init_progress <- function(progress_steps)
{
  max <- progress_steps[length(progress_steps)]
  pb <- txtProgressBar(min = 0, max, style = 3)
  return(pb)
}

#================================================================================
#' Update progress reporting
#'
#' @param pb A progress bar
#' @param progress_steps Vector of values to use as progress bar updates
#' (should match number of calls being made to update_progress)
#'
update_progress <- function(pb, progress_steps)
{
  # get current value of progress bar
  current <- getTxtProgressBar(pb)

  # work out index of this value in progress_steps
  if (current == 0)
  {
    thisstep <- 0
  }
  else
  {
    thisstep <- match(current, progress_steps)
  }

  # set progress bar to next value
  setTxtProgressBar(pb, progress_steps[thisstep+1])
  return(pb)
}
