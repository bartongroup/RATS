#================================================================================
#' Import abundances directly from salmon and kallisto output.
#' 
#' Rather than re-inventing the wheel, this function is a wrapper for the wasabi converter,
#' the sleuth data loader and RATs' sleuth data extractor. This also has the advantage that the
#' sleuth loader scales TPMs to estimated transcript counts, which is the appropriate abundance
#' type for RATs. TPMs directly are not suitable, neither are read counts.
#' 
#' Converting, normalising and importing multiple bootstrapped abundance files takes a bit of time.
#' 
#' \code{wasabi} automatically skips format conversion if a folder already contains an \code{abundance.h5} file.
#' 
#' @param A_paths (character) A vector of strings, listing the directory paths to the quantifications for the first condition. One directory per replicate. The directory name should be a unique identifier for the sample.
#' @param B_paths (character) A vector of strings, listing the directory paths to the quantifications for the second condition. One directory per replicate. The directory name should be a unique identifier for the sample.
#' @param annot (data.frame) A table matching transcript identifiers to gene identifiers. This should be the same that you will use with \code{call_DTU()}.
#' @param TARGET_COL (character) The name of the column in \code{annot} that contains the transcript IDs. (Default "target_id")
#' @param half_cooked (logical) If TRUE, input is already in \code{Kallisto} h5 format and \code{wasabi} conversion will be skipped. (Default FALSE)
#' @return A list of two, representing the abundances per condition. These will be formatted in the RATs generic bootstrapped data input format.
#' 
#' @import wasabi
#' @import sleuth
#' @import data.table
#' @export
#' 
fish4rodents <- function(A_paths, B_paths, half_cooked=FALSE, TARGET_COL="target_id") {
  # Wasabi?
  if (!half_cooked) {
    prepare_fish_for_sleuth(c(A_paths, B_paths))
  }
  
  # Sleuth loader.
  smpl2cond <- data.frame(sample=c(sapply(A_paths, basename), sapply(B_paths, basename)),
                          condition=c(rep.int("A", length(A_paths)), rep.int("B", length(B_paths))),
                          path=c(A_paths, B_paths), 
                          stringsAsFactors=FALSE)
  slo <- sleuth_prep(smpl2cond, ~ condition)
  
  # Extract.
  samples_by_condition <- group_samples(slo$sample_to_covariates)
  
  boots_data_A <- denest_sleuth_boots(slo, annot, samples_by_condition$condition[[name_A]], TARGET_COL= TARGET_COL)
  boots_data_B <- denest_sleuth_boots(slo, annot, samples_by_condition$condition[[name_B]], TARGET_COL= TARGET_COL)
    
  return(list("boot_data_A"= boots_data_A, "boot_data_B"= boots_data_B))
}
