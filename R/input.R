########## ########## ########## ########## ########## ########## ########## ########## ##########
# Helper functions for setting up the input, if necessary.
########## ########## ########## ########## ########## ########## ########## ########## ##########


#================================================================================
#' Extract matching transcript and gene IDs from a GRanges object.
#'
#' The GRanges must have at least the following metadata columns: `gene_id` and `transcript_id` (such as GRanges imported from GTF).
#' GFF3-style metadata columns are currently not supported.
#'
#' @param annot A GRanges object with gene_id and transcript_id metadata columns.
#' @param transc_header The title for the transcripts column in the output. (target_id)
#' @param gene_header The title for the genes column in the output. (parent_id)
#' @return A data.table with two columns, matching transcript IDs to gene IDs.
#'
#' @importFrom GenomicRanges mcols
#' @import data.table
#' @export
#'
granges2ids <- function(annot, transc_header= "target_id", gene_header= "parent_id")
{
  # If GTF style...
  if ('transcript_id' %in% names(mcols(annot))) {
    t2g <- unique(data.table(transcript_id=annot$transcript_id, gene_id=annot$gene_id))
  # If GFF3 style...
  } else {
    stop('It seems you supplied a GFF3 file. This is currenlty not implemented. Please convert to GTF using available tools and try again.')
  }
  # Don't want any unpaired IDs (genes without transcripts are irrelevant, transcripts without gene are unusable).
  t2g <- t2g[!is.na(t2g$gene_id) & !is.na(t2g$transcript_id), ]    # R Check really does not like data.table syntax... :(
  
  names(t2g) <- c(transc_header, gene_header)
  return(t2g)
}


#================================================================================
#' Extract matching transcript and gene IDs from a GTF file.
#'
#'Previously named annot2ids(). The old name will be discontinued.
#'
#' GFF3 not supported.
#'
#' @param annotfile A GTF file.
#' @param transc_header The title for the transcripts column in the output. (target_id)
#' @param gene_header The title for the genes column in the output. (parent_id)
#' @return A data.table with two columns, matching transcript IDs to gene IDs.
#'
#' @importFrom rtracklayer import
#' @export
#'
gtf2ids <- function(annotfile, transc_header= "target_id", gene_header= "parent_id")
{
  annot <- rtracklayer::import(annotfile)
  granges2ids(annot, transc_header, gene_header)
}


#================================================================================
#' (deprecated) Extract matching transcript and gene IDs from a GTF file.
#'
#' This function has been renamed to gtf2ids(). The old name will be discontinued.
#' 
#' GFF3 not supported.
#'
#' @param annotfile A GTF file.
#' @param transc_header The title for the transcripts column in the output. (target_id)
#' @param gene_header The title for the genes column in the output. (parent_id)
#' @return A data.table with two columns, matching transcript IDs to gene IDs.
#'
#' @export
#'
annot2ids <- function(annotfile, transc_header= "target_id", gene_header= "parent_id")
{
  warning("This function has been renamed to 'gtf2ids()' for clarity. The old name 'annot2ids()' will be discontinued.")
  gtf2ids(annotfile, transc_header, gene_header)
}


#================================================================================
#' Prepare GTF for gene model plotting.
#' 
#' Create GRanges objects for each transcript and bind them into a GRangeList for each gene.
#' Then one can easily \code{ggbio::autoplot()} a gene's models by its gene ID.
#' 
#' @param annotfile A GTF file.
#' @param threads Number of processing threads (Default 1).
#' @return A list (by gene) of GRangeLists (by transcript).
#'
#' @importFrom rtracklayer import
#' @import parallel
#' @importFrom GenomicRanges makeGRangesListFromDataFrame
#' @export
#'
annot2models <- function(annotfile, threads= 1L)
{
  message('This will take a few minutes...')
  gr <- rtracklayer::import(annotfile)
  
  # Split entries by gene and then by transcript.
  lg <- list()
  
  # If GTF style...
  if ('transcript_id' %in% names(mcols(gr))) {
    genes <- unique(gr$gene_id)
    na <- is.na(genes)
    if (any(na)){
      warning("The annotation contains entries without gene_id assignment. These will be ommitted.")
      genes <- genes[!na]
      gr <- gr[!is.na(gr$gene_id)]
    }
    lg <- mclapply(genes, function(g){
      features <- gr[gr$gene_id == g]
      isoforms <- makeGRangesListFromDataFrame(as.data.frame(features), split.field='transcript_id', keep.extra.columns=TRUE)
      
      return(isoforms)
    }, mc.cores = threads)
    names(lg) <- genes
  # If GFF3 style...
  } else {
    stop("It looks like you supplied a GFF3 fileinstead of GTF. This is currently disabled as the result does not plot properly. Please convert to GTF first, using available tools, and try again.")
  }
  
  return(lg)
}


#================================================================================
#' Import abundances directly from kallisto output.
#'
#' Apply TPM normalisation using the info available from the abundance.h5 files.
#'
#' Converting, normalising and importing multiple bootstrapped abundance files takes a bit of time.
#' IMPORTANT: This function is currently not intended to be used to import non-bootstrapped quantifications.
#'
#' @param A_paths (character) A vector of strings, listing the directory paths to the quantifications for the first condition. One directory per replicate, without trailing path dividers. The directory name should be a unique identifier for the sample.
#' @param B_paths (character) A vector of strings, listing the directory paths to the quantifications for the second condition. One directory per replicate, without trailing path dividers.. The directory name should be a unique identifier for the sample.
#' @param annot (data.frame) A table matching transcript identifiers to gene identifiers. This should be the same that you used for quantification and that you will use with \code{call_DTU()}. It is used to order the transcripts consistently throughout RATs.
#' @param scaleto (double) Scaling factor for normalised abundances. (Default 1000000 gives TPM). If a numeric vector is supplied instead, its length must match the total number of samples. The value order should correspond to the samples in group A followed by group B. This allows each sample to be scaled to its own actual library size, allowing higher-throughput samples to carry more weight in deciding DTU.
#' @param beartext (logical) Instead of importing bootstrap data from the \code{abundance.h5} file of each sample, import it from plaintext files in a \code{bootstrap} subdirectory created by running \code{kallisto}'s \code{h5dump} subcommand (Default FALSE). This workaround circumvents some mysterious .h5 parsing issues on certain systems.
#' @param TARGET_COL The name of the column for the transcript identifiers in \code{annot}. (Default \code{"target_id"})
#' @param PARENT_COL The name of the column for the gene identifiers in \code{annot}. (Default \code{"parent_id"})
#' @param threads (integer) For parallel processing. (Default 1)
#' @return A list of two, representing the TPM abundances per condition. These will be formatted in the RATs generic data input format, preferably for bootstrapped estimates (if bootstraps are available) or otherwise for plain count estimates.
#'
#' @import parallel
#' @import data.table
#' @import rhdf5
#'
#' @export

fish4rodents <- function(A_paths, B_paths, annot, TARGET_COL="target_id", PARENT_COL="parent_id", beartext=FALSE, threads= 1L, scaleto= 1000000)
{
  threads <- as.integer(threads)  # Can't be decimal.
  setDTthreads(threads, restore_after_fork = TRUE)

  # Don't assume the annotation is ordered properly.
  annot <- tidy_annot(annot, TARGET_COL, PARENT_COL)

  # Sort out scaling factor per sample.
  lgth <- c("A"=length(A_paths), "B"=length(B_paths))
  sfac <- list("A"=NA_real_, "B"=NA_real_)
  if (length(scaleto) == 1) {     # uniform scaling
    sfac$A <- rep(scaleto, lgth["A"])
    sfac$B <- rep(scaleto, lgth["B"])
  } else {                        # individual scaling
    sfac$A <- scaleto[1:lgth$A]
    sfac$B <- scaleto[(1 + lgth$A):(lgth$A + lgth$B)]
  }

  # Load and convert.
  res <- lapply(c('A', 'B'), function(cond) {
    boots <- mclapply(1:lgth[cond], function(x) {
      # Get the respective file and scaling factor.
      sf <- sfac[[cond]][x]
      if (cond=='A') {
        fil <- A_paths[x]
      } else {
        fil <- B_paths[x]
      }
      
      # If data from Kallisto plaintext...
      if (beartext) {
        # list the bootstrap files
        bfils <- list.files(path=file.path(fil, "bootstraps"), full.names=TRUE, no..=TRUE)
        bfils <- bfils[ grepl('bs_abundance', bfils) ]
        # parse info.
        meta <- fread(bfils[[1]], header=TRUE)[, c('target_id', 'eff_length'), with=FALSE]  
            # the IDs should all come out in the same order in every iteration file of a given sample, 
            # and the transcript lengths should not change either.
        counts <- as.data.table( lapply(bfils, function(bf){ fread(bf, header=TRUE)[['est_counts']] }) )  # plaintext already has TPMs computed, but I stick with the raw counts
                                                                                                          # for consistency with the .h5 mode and to allow normalisation to other target sizes 
      # If data from Salmon/Wasabi or Kallisto abundance.h5...
      } else {
        content <- rhdf5::h5ls(file.path(fil, "abundance.h5"))
        meta <- as.data.table(list( target_id=as.character(rhdf5::h5read(file.path(fil, "abundance.h5"), "/aux/ids")), 
                                    eff_length=as.numeric(rhdf5::h5read(file.path(fil, "abundance.h5"), "/aux/eff_lengths")) ))
        
        if ('bootstrap' %in% content$name) {
          counts <- as.data.table( rhdf5::h5read(file.path(fil, "abundance.h5"), "/bootstrap") )
          # Normalise raw counts.
          tpm <- as.data.table( lapply(counts, function (y) {
            cpb <- y / meta$eff_length
            tcpb <- sf / sum(cpb)
            return(cpb * tcpb)
          }) )
          # Structure.
          dt <- as.data.table( cbind(meta$target_id, tpm) )
          with(dt, setkey(dt, V1) )
          names(dt)[1] <- TARGET_COL
          # Order transcripts to match annotation.
          dt <- merge(annot[, c(TARGET_COL), with=FALSE], dt, by=TARGET_COL, all=TRUE)
          
        } else {
          counts <- rhdf5::h5read(file.path(fil, "abundance.h5"), "/est_counts")
          # Normalize.
          cpb <- counts / meta$eff_length
          tcpb <- sf / sum(cpb)
          tpm <- cpb * tcpb
          # Structure
          dt <- as.data.table( cbind(meta$target_id, tpm) )
          with(dt, setkey(dt, V1) )
          names(dt)[1] <- TARGET_COL
          names(dt)[2] <- basename(fil)
          # Order transcripts to match annotation.
          dt <- merge(annot[, c(TARGET_COL), with=FALSE], dt, by=TARGET_COL, all=TRUE)
        }
        
        return (dt)
      }
      return (dt)
    }, mc.cores = min(threads,lgth[cond]), mc.preschedule = TRUE, mc.allow.recursive = FALSE)
    
    # If single measurement for all samples, ie. est_counts instead of bootsraps, merge and unnest to meet un-bootsrapped format.
    if( sum(vapply(boots, function(b){length(b)-1}, numeric(1))) == length(boots) ) {
      boots <- Reduce(merge, boots)
    }
    
    return (boots)
  })
  
  if (is.data.table(res[[1]])) {
    names(res) <- c("count_data_A", "count_data_B")
  } else {
    names(res) <- c("boot_data_A", "boot_data_B")
  }
  
  return(res)
}


#================================================================================
#' Extract bootstrap counts into a less nested structure.
#'
#' *Legacy function*
#'
#' It extracts the bootstrap data from the older-style \code{sleuth} object.
#' As of sleuth version 0.29, the bootstrap data is no longer kept in the object.
#'
#' @param slo A sleuth object.
#' @param annot A data.frame matching transcript identfier to gene identifiers.
#' @param samples A numeric vector of samples to extract counts for.
#' @param COUNTS_COL The name of the column with the counts. (Default "tpm")
#' @param BS_TARGET_COL The name of the column with the transcript IDs. (Default "target_id")
#' @param TARGET_COL The name of the column for the transcript identifiers in \code{annot}. (Default \code{"target_id"})
#' @param PARENT_COL The name of the column for the gene identifiers in \code{annot}. (Default \code{"parent_id"})
#' @param threads Number of threads. (Default 1)
#' @return A list of data.tables, one per sample, containing all the bootstrap counts of the smaple. First column contains the transcript IDs.
#'
#' NA replaced with 0.
#'
#' Transcripts in \code{slo} that are missing from \code{annot} will be skipped completely.
#' Transcripts in \code{annot} that are missing from \code{slo} are automatically padded with NA, which we re-assign as 0.
#'
#' @import parallel
#' @import data.table
#' @export
#'
denest_sleuth_boots <- function(slo, annot, samples, COUNTS_COL="tpm", BS_TARGET_COL="target_id",
                                TARGET_COL="target_id", PARENT_COL="parent_id", threads= 1) {
  
  warning("DEPRECATED: denest_sleuth_boots() no longer serves a purpose in RATs and will be eventually removed.")
  
  # Ensure data.table complies.
  setDTthreads(threads, restore_after_fork = TRUE)

  # Compared to the size of other structures involved in calling DTU, this should make negligible difference.
  tx <- tidy_annot(annot, TARGET_COL, PARENT_COL)$target_id

  mclapply(samples, function(smpl) {
    # Extract counts in the order of provided transcript vector, for safety and consistency.
    dt <- as.data.table( lapply(slo$kal[[smpl]]$bootstrap, function(boot) {
      roworder <- match(tx, boot[[BS_TARGET_COL]])
      boot[roworder, COUNTS_COL]
    }))
    # Replace any NAs with 0. Happens when annotation different from that used for DTE.
    dt[is.na(dt)] <- 0
    # Add transcript ID.
    with(dt, dt[, target_id := tx])
    nn <- names(dt)
    ll <- length(nn)
    # Return reordered so that IDs are in first column.
    return(dt[, c(nn[ll], nn[seq.int(1, ll-1)]), with=FALSE])
  }, mc.cores= threads, mc.allow.recursive= FALSE, mc.preschedule= TRUE)
}

