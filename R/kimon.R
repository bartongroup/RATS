#' Compute mean, variance, etc per transcript across bootstraps for a single replicate
#'
#' @param boots a list of bootstrap dataframes for a single replicate.
#' @param targets a vector of target_id (pre-filtered to exclude id without AS siblings).
#' @param count_type either 'est_counts' or 'tpm'
#' @export


# Faking the stuff the function expects to receive.
targets <- ids$targets <- ids$target_id[ids$has_siblings]
sl <- At_250genes_50singletx_200multitx_513totaltx   # sleuth object
boots <- sl$kal[[1]]$bootstrap
count_type <- "tpm"

one_rep_stats <- function(boots, targets, count_type="tpm") {
  # Gather counts into one dataframe
  gather <- data.frame(row.names=targets)
  for (b in 1:length(boots)) {
    gather[paste("boot",b)] <- boots[[b]][targets,count_type]
  }
  # Stats by row
  metrics <- stat.desc(t(gather))
}



# Testing how to access things...
sl <- At_250genes_50singletx_200multitx_513totaltx   # sleuth object
r <- sl$kal[[1]]                                     # kallisto object for R-eplicate 1
rb <- r$bootstrap[[1]]                               # dataframe of transcript results for B-ootstrap 1 of replicate 1
rbt <- rb$target_id[1]                               # id value for T-ranscript 1 in bootstrap 1 of replicate 1.







