#' {Arabidopsis thaliana} TAIR10 annotation information
#'
#' Annotation information extracted from ensembl with \code{\link{biomaRt}} that links
#' transcript level annotations with the parent gene level annotations and an external
#' annotation.
#'
#' The data were generated with: \preformatted{
#' library("biomaRt")
#' mart <- biomaRt::useMart(biomart = "plants_mart", host="plants.ensembl.org", dataset = "athaliana_eg_gene")
#' At_tair10_t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_id"), mart = mart)
#' At_tair10_t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, parent_id = ensembl_gene_id, other_id_1 = external_gene_id)
#' devtools::use_data(At_tair10_t2g)
#' }
#'
#' @format A data frame with 41,671 observations of 3 variables
#' \describe{
#'   \item{target_id}{ENSEMBL transcript id.}
#'   \item{parent_id}{ENSEMBL gene id.}
#'   \item{other_id_1}{ENSEMBL external id.}
#'
#' }
"At_tair10_t2g"

#' Example sleuth object containing transcript abundance estimates and bootstrap data
#'
#' The transcript abundance estimates and bootstraps were calculated with
#' \href{http://combine-lab.github.io/salmon/}{salmon} for 20 biological replicates in 3
#' conditions:
#' \itemize{
#'   \item{Col: 6 replicates}
#'   \item{Vir: 7 replicates}
#'   \item{L: replicates}
#' }
#' The \href{http://combine-lab.github.io/salmon/}{salmon} abundance information was
#' converted to the correct hd5 format using \code{\link[wasabi]{prepare_fish_for_sleuth}}.
#' The converted files are read by \code{\link[sleuth]{sleuth_prep}} and combined with the
#' experiment metadata (including a potental batch effect). to produce the output
#' \code{\link{sleuth}} object. A sleuth object contains:
#' \describe{
#'   \item{sample_to_covariates}{A data frame with 20 observations of 3 variables relating the
#'   sample information to the 'covariates' int he model (i.e. the conditions and, in this case,
#'   the batch.)}
#'   \item{full_formula}{The formula that will be used by sleuth for doing DTE.}
#'   \item{design_matrix}{The experiemntal design mapping covariates to samples.}
#'   \item{target_mapping}{Currently NA - probably there target_ids eventually get mapped to
#'   genes.}
#'   \item{kal}{A list of 20 kallisto objects (1/sample).}
#'   \item{kal_versions}{The program that make the kallisto objects.}
#'   \item{obs_raw}{A dataframe with all the raw data in (707,720 observations of 6 variables).}
#'   \item{obs_norm}{A dataframe with normalized data in (707,720 observations of 6 variables).}
#'   \item{filter_df}{A data frame listing those transcripts that pass the filtering
#'   criteria. The filtering applied here requires that a transcript has >5 counts in at
#'   least 47\% of the samples.}
#'   \item{filter_bool}{Logical array marking those transcripts that pass the filtering
#'   criteria.}
#'   \item{obs_norm_filt}{A dataframe with filtered normalized data in (552,220 observations
#'   of 6 variables).}
#'   \item{est_counts_sf}{Not quite sure what this is!}
#'   \item{tpm_sf}{Not quite sure what this is!}
#'   \item{bootstrap_summary}{Currently NA.}
#'   \item{bs_summary}{A list of two objects. More about this below!}
#' }
#' For our purposes the realy interesing bits here are the kallisto objects, the
#' obs_norm_filt data and the summary_bs data. The kallisto objects are essentially lists of
#' five objects:
#' \describe{
#'   \item{abundance}{A data frame with 35,386 observations of 6 variables that contains the
#'   data for this sample. This data includes the target ID (\code{target_id}), transcript
#'   length (\code{len}), the transcript effective length (\code{eff_len}) the estimated
#'   transcript counts (\code{est_counts}), and the transcripts-per-million values
#'   (\code{tpm}).}
#'   \item{bias_normalized}{Logical NA}
#'   \item{bias_observed}{Logical NA}
#'   \item{bootstrap}{A list of 100 objects, each of which is a data frame with
#'   35,386 observations of 3 variables that contains the target ID (\code{target_id}),
#'   estimated transcript counts (\code{est_counts}), and the transcripts-per-million
#'   values (\code{tpm}) from a bootstrap.}
#' }
#' \code{summary_bs} is a list containing two objects:
#' \describe{
#'   \item{obs_counts}{A 27,611 x 20 matrix of (I guess) the variance of each transcripts
#'   estimated counts across the bootstraps in each sample.}
#'   \item{sigma_q_sq}{A vector with the sigma-squared values for each transcript across
#'   the replicate bootraps.}
#' }
#' @format A sleuth object.
"At_RNAMeth_sleuth"
