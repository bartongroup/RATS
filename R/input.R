#================================================================================
#' Extract matching transcript and gene IDs from a GTF file.
#' 
#' This function performs no file format validation. Please use the appropriate type.
#' It is NOT compatible GFF3, as there is no rigid specification for ID attributes.
#' 
#' @param annotfile A GTF file.
#' @param transc_header The tile for the transcripts column in the output. (target_id)
#' @param gene_header The tile for the genes column in the output. (parent_id)
#' @return A data.table with two columns, matching transcript IDs to gene IDs.
#' 
#' @export
#'
annot2ids <- function(annotfile, transc_header= "target_id", gene_header= "parent_id")
{
  annot <- as.data.table(read.delim(annotfile, comment.char= "#", header= FALSE))
  annot <- annot[(grepl("transcript_id", annot[[9]]) & grepl("gene_id", annot[[9]])), ]
  t2g <- data.table("transc_id" = gsub(".*transcript_id \"?(.+?)\"?;.*", "\\1", annot[[9]]),
                    "gene_id" = gsub(".*gene_id \"?(.+?)\"?;.*", "\\1", annot[[9]]) )
  # Clean up.
  t2g <- unique(t2g)
  names(t2g) <- c(transc_header, gene_header)
  
  return(t2g)
}


#================================================================================
#' Import abundances directly from salmon and kallisto output.
#' 
#' Uses \code{wasabi} to convert Salmon read counts format to Kallisto RHDF5 format,
#' then applies TPM normalisation using the info available from the abundance.h5 files.
#' 
#' Converting, normalising and importing multiple bootstrapped abundance files takes a bit of time.
#' 
#' \code{wasabi} automatically skips format conversion if a folder already contains an \code{abundance.h5} file.
#' 
#' @param A_paths (character) A vector of strings, listing the directory paths to the quantifications for the first condition. One directory per replicate. The directory name should be a unique identifier for the sample.
#' @param B_paths (character) A vector of strings, listing the directory paths to the quantifications for the second condition. One directory per replicate. The directory name should be a unique identifier for the sample.
#' @param annot (data.frame) A table matching transcript identifiers to gene identifiers. This should be the same that you used for quantification and that you will use with \code{call_DTU()}. It is used to order the transcripts consistently throughout RATs.
#' @param TARGET_COL The name of the column for the transcript identifiers in \code{annot}. (Default \code{"target_id"})
#' @param PARENT_COL The name of the column for the gene identifiers in \code{annot}. (Default \code{"parent_id"})
#' @param half_cooked (logical) If TRUE, input is already in \code{Kallisto} h5 format and \code{wasabi} conversion will be skipped. Wasabi automatically skips conversion if abundance.h5 is present, so this parameter is redundant, unless wasabi is not installed. (Default FALSE)
#' @param threads (integer) For parallel processing. (Default 1)
#' @return A list of two, representing the TPM abundances per condition. These will be formatted in the RATs generic bootstrapped data input format.
#' 
#' @import wasabi
#' @import rhdf5
#' @import data.table
#' @import parallel
#'
#' @export
#' 
fish4rodents <- function(A_paths, B_paths, annot, TARGET_COL="target_id", PARENT_COL="parent_id", half_cooked=FALSE, threads= 1) 
{
  threads <- as.integer(threads)  # Can't be decimal.
  # Ensure data.table complies.
  # if (packageVersion("data.table") >= "1.9.8")
    setDTthreads(1)
  
  # Don't assume the annotation is ordered properly.
  annot <- tidy_annot(annot, TARGET_COL, PARENT_COL)
  
  # Wasabi?
  if (!half_cooked) {
    prepare_fish_for_sleuth(c(A_paths, B_paths))
  }
  
  # Sleuth loader.
  # smpl2cond <- data.frame(sample=c(sapply(A_paths, basename), sapply(B_paths, basename)),
  #                         condition=c(rep.int("A", length(A_paths)), rep.int("B", length(B_paths))),
  #                         path=c(A_paths, B_paths), 
  #                         stringsAsFactors=FALSE)
  # slo <- sleuth_prep(smpl2cond, ~ condition)
  # Extract.
  # samples_by_condition <- group_samples(slo$sample_to_covariates)
  # 
  # boots_data_A <- denest_sleuth_boots(slo, annot, samples_by_condition$condition[["A"]], TARGET_COL= TARGET_COL, PARENT_COL=PARENT_COL, BS_TARGET_COL=BS_TARGET_COL, COUNTS_COL=COUNTS_COL)
  # boots_data_B <- denest_sleuth_boots(slo, annot, samples_by_condition$condition[["B"]], TARGET_COL= TARGET_COL, PARENT_COL=PARENT_COL, BS_TARGET_COL=BS_TARGET_COL, COUNTS_COL=COUNTS_COL)
  
  # Load and convert manually.
  boots_A <- mclapply(A_paths, function(x) {
    ids <- as.data.table( h5read(file.path(x, "abundance.h5"), "/aux/ids") )
    counts <- as.data.table( h5read(file.path(x, "abundance.h5"), "/bootstrap") )
    effl <- h5read(file.path(x, "abundance.h5"), "/aux/eff_lengths")
    tpm <- as.data.table( lapply(counts, function (y) {
      cpb <- y / effl
      tcpb <- 1000000 / sum(cpb)
      return(cpb * tcpb)
    }) )
    dt <- as.data.table( cbind(ids, tpm) )
    with(dt, setkey(dt, V1) )
    names(dt)[1] <- TARGET_COL
    # Order transcripts to match annotation.
    ord <- match(annot$target_id, dt$target_id)
    return (dt[ord, ])
  }, mc.cores = threads, mc.preschedule = TRUE, mc.allow.recursive = FALSE)
  boots_B <- mclapply(B_paths, function(x) {
    ids <- as.data.table( h5read(file.path(x, "abundance.h5"), "/aux/ids") )
    counts <- as.data.table( h5read(file.path(x, "abundance.h5"), "/bootstrap") )
    effl <- h5read(file.path(x, "abundance.h5"), "/aux/eff_lengths")
    tpm <- as.data.table( lapply(counts, function (y) {
      cpb <- y / effl
      tcpb <- 1000000 / sum(cpb)
      return(cpb * tcpb)
    }) )
    dt <- as.data.table( cbind(ids, tpm) )
    with(dt, setkey(dt, V1) )
    names(dt)[1] <- "target_id"
    # Order transcripts to match annotation.
    ord <- match(annot$target_id, dt$target_id)
    return (dt[ord, ])
  }, mc.cores = threads, mc.preschedule = TRUE, mc.allow.recursive = FALSE)
    
  return(list("boot_data_A"= boots_A, "boot_data_B"= boots_B))
}











# annot <- read.delim("~/Homo_sapiens.GRCh37.60.t2g.tsv")
# 
# basedir <- "/Volumes/kfroussios/PROJECTS/rats"
# quantdir <- "salmon_quant"
# 
# A_paths <- file.path(basedir, quantdir, c("SRR393010", "SRR393011", "SRR393012"))
# B_paths <- file.path(basedir, quantdir, c("SRR393013", "SRR393014", "SRR393015"))
# 
# h5ls(file.path(A_paths[[1]], "abundance.h5"))
# 
# h5read(file.path(A_paths[[1]], "abundance.h5"), "/bootstrap")
# h5read(file.path(A_paths[[1]], "abundance.h5"), "/aux/ids")
# h5read(file.path(A_paths[[1]], "abundance.h5"), "/aux/fld")
# 
# H5close()
# 
# col <-  h5read(file.path(A_paths[[1]], "abundance.h5"), "/bootstrap/bs0")
# effl <- h5read(file.path(A_paths[[1]], "abundance.h5"), "/aux/eff_lengths")
# 
# cpb <- col / effl
# tcpb  <- 1000000 / sum(cpb)
# tpm <- cpb * tcpb
# 
# head(col)
# head(tpm)
# 
# 
