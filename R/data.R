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

##' Example sleuth object containing transcript abundance estimates and bootstrap data
##'
##' The transcript abundance estimates and bootstraps were calculated with
##' \href{http://combine-lab.github.io/salmon/}{salmon} for 20 biological replicates in 3
##' conditions:
##' \itemize{
##'   \item{Col: 6 replicates}
##'   \item{Vir: 7 replicates}
##'   \item{L: replicates}
##' }
##' The \href{http://combine-lab.github.io/salmon/}{salmon} abundance information was
##' converted to the correct hd5 format using \code{\link[wasabi]{prepare_fish_for_sleuth}}.
##' The converted files are read by \code{\link[sleuth]{sleuth_prep}} and combined with the
##' experiment metadata (including a potental batch effect). to produce the output
##' \code{\link{sleuth}} object. A sleuth object contains:
##' \describe{
##'   \item{sample_to_covariates}{A data frame with 20 observations of 3 variables relating the
##'   sample information to the 'covariates' int he model (i.e. the conditions and, in this case,
##'   the batch.)}
##'   \item{full_formula}{The formula that will be used by sleuth for doing DTE.}
##'   \item{design_matrix}{The experiemntal design mapping covariates to samples.}
##'   \item{target_mapping}{Currently NA - probably there target_ids eventually get mapped to
##'   genes.}
##'   \item{kal}{A list of 20 kallisto objects (1/sample).}
##'   \item{kal_versions}{The program that make the kallisto objects.}
##'   \item{obs_raw}{A dataframe with all the raw data in (707,720 observations of 6 variables).}
##'   \item{obs_norm}{A dataframe with normalized data in (707,720 observations of 6 variables).}
##'   \item{filter_df}{A data frame listing those transcripts that pass the filtering
##'   criteria. The filtering applied here requires that a transcript has >5 counts in at
##'   least 47\% of the samples.}
##'   \item{filter_bool}{Logical array marking those transcripts that pass the filtering
##'   criteria.}
##'   \item{obs_norm_filt}{A dataframe with filtered normalized data in (552,220 observations
##'   of 6 variables).}
##'   \item{est_counts_sf}{Not quite sure what this is!}
##'   \item{tpm_sf}{Not quite sure what this is!}
##'   \item{bootstrap_summary}{Currently NA.}
##'   \item{bs_summary}{A list of two objects. More about this below!}
##' }
##' For our purposes the realy interesing bits here are the kallisto objects, the
##' obs_norm_filt data and the summary_bs data. The kallisto objects are essentially lists of
##' five objects:
##' \describe{
##'   \item{abundance}{A data frame with 35,386 observations of 6 variables that contains the
##'   data for this sample. This data includes the target ID (\code{target_id}), transcript
##'   length (\code{len}), the transcript effective length (\code{eff_len}) the estimated
##'   transcript counts (\code{est_counts}), and the transcripts-per-million values
##'   (\code{tpm}).}
##'   \item{bias_normalized}{Logical NA}
##'   \item{bias_observed}{Logical NA}
##'   \item{bootstrap}{A list of 100 objects, each of which is a data frame with
##'   35,386 observations of 3 variables that contains the target ID (\code{target_id}),
##'   estimated transcript counts (\code{est_counts}), and the transcripts-per-million
##'   values (\code{tpm}) from a bootstrap.}
##' }
##' \code{summary_bs} is a list containing two objects:
##' \describe{
##'   \item{obs_counts}{A 27,611 x 20 matrix of (I guess) the variance of each transcripts
##'   estimated counts across the bootstraps in each sample.}
##'   \item{sigma_q_sq}{A vector with the sigma-squared values for each transcript across
##'   the replicate bootraps.}
##' }
##' @format A sleuth object.
#"At_RNAMeth_sleuth"

#' Small example sleuth object containing transcript abundance estimates and bootstrap data
#'
#' \href{http://combine-lab.github.io/salmon/}{Salmon} abundance estimates and bootstraps were
#' calculated for 513 Arabidopsis thaliana transcripts from the TAIR10 annotation. These
#' transcripts map to 250 genes - 50 ofwhich are single-transcript genes, 200 of which have
#' multiple transcripts (463 in total).
#'
#' The data were computed for 3 biological replicates in 2 conditions:
#' \itemize{
#'   \item{Col: 3 replicates}
#'   \item{Vir: 3 replicates}
#' }
#'
#' The \href{http://combine-lab.github.io/salmon/}{salmon} abundance information is then
#' processed with the routine \code{\link{gen_example_sleuth}} (including a convertion to
#' the correct hd5 format if required) to generate theoutput sleuth object that is saved here.
#'
#' See the help page for \code{\link{At_RNAMeth_sleuth}} for a description of the sleuth
#' object structure.
#'
#' @format A sleuth object.
"At_250genes_50singletx_200multitx_513totaltx"

#' A function that generates the example sleuth object storred in
#' At_250genes_50singletx_200multitx_513totaltx.rda from the raw salmon 0.5.1. data
#'
#' @param overwrite Logical flag for triggering an overwrite of At_250genes_50singletx_200multitx_513totaltx.rda.
#' @param annot name of the annotation used for the data - assumes this is part of the path to the files.
#' @param path root path to the raw salmon data - combined with annot to get the data for all replicates.
#'
#' @export
gen_example_sleuth <- function(overwrite=FALSE, annot="TAIR10",
                               path="/cluster/gjb_lab/nschurch/Projects/Arabidopsis_RNAMeth/salmon/rats_testdata_0.5.1_bootstraps"){

  # construct paths
  files <- Sys.glob(file.path(path,"*-*", annot))

  # convert the sailfish dirs to kallisto h5 format
  wasabi::prepare_fish_for_sleuth(files)

  # make the table for the data
  colreps <- sapply(1:3, function(repname) paste("Col-",repname, sep=""))
  virreps <- sapply(1:3, function(repname) paste("Vir-",repname, sep=""))
  allreps <- c(colreps, virreps)

  # construct condition metadata
  colcond <- sapply(1:3, function(repname) paste("Col"))
  vircond <- sapply(1:3, function(repname) paste("Vir"))
  allcond <- c(colcond, vircond)

  # construct batch effect metadata
  batch1 <- sapply(1:4, function(repname) paste("batch1"))
  batch2 <- sapply(1:2, function(repname) paste("batch2"))
  batch <- c(batch1, batch2)

  # initiate data frame (remembering to skip col4 cos its bad)
  exp <- matrix(, ncol=4, nrow=(length(allreps)))
  colnames(exp) <- c("sample","condition","path","batch")
  exp[,1] <- allreps
  exp[,2] <- allcond
  exp[,3] <- files
  exp[,4] <- batch
  exp <- as.data.frame(exp, stringsAsFactors=FALSE)

  # construct sleuth object
  so <- sleuth::sleuth_prep(exp, ~ condition+batch)

  # fit linear model
  so <- sleuth::sleuth_fit(so)
  sleuth::models(so)
  so <- sleuth::sleuth_wt(so, 'conditionVir')

  if (overwrite) {
    # save the sleuth object includind the sleuth results for this data.
    At_250genes_50singletx_200multitx_513totaltx <- so
    devtools::use_data(At_250genes_50singletx_200multitx_513totaltx, overwrite=TRUE)
  }

  return(so)
}
