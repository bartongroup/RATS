########## ########## ########## ########## ########## ########## ########## ########## ##########
# The functions in this file are mainly helper functions, not intended to be invoked directly.
########## ########## ########## ########## ########## ########## ########## ########## ##########


#================================================================================
#' Check input parameters.
#'
#' @param annot Annotation dataframe.
#' @param correction P-value correction method.
#' @param p_thresh Significance level.
#' @param TARGET_COL Name of transcript id column in annotation.
#' @param PARENT_COL Name of gene id column in annotation.
#' @param abund_thresh Minimum transcript abundance per sample.
#' @param testmode Which test to run.
#' @param qboot Whether to bootstrap against quantifications.
#' @param qbootnum Number of bootstrap iterations.
#' @param qrep_thresh Confidence threshold.
#' @param dprop_thresh Minimum change in proportion.
#' @param count_data_A A dataframe of estimated counts.
#' @param count_data_B A dataframe of estimated counts.
#' @param boot_data_A A list of dataframes, one per sample, each with all the bootstrapped estimates for the sample.
#' @param boot_data_B A list of dataframes, one per sample, each with all the bootstrapped estimates for the sample.
#' @param rboot Whether to bootstrap against samples.
#' @param rrep_thresh Confidence threshold.
#' @param threads Number of threads.
#' @param seed Seed for random engine.
#' @param scaling Abundance scaling factor.
#' @param reckless whether to ignore detected annotation discrepancies.
#' @param lean whether to report descriptive statistics for Dprop and pval.
#' @param verbose Verbose status.
#' @param use_sums Whether to use sums or means.
#'
#' @return List: \itemize{
#'  \item{"error"}{logical}
#'  \item{"message"}{string}
#'  \item{"maxboots"}{Infered rule-of-thumbs number of iterations to do.}
#'  \item{"warn"}{logical}
#'  \item{"warnings}{list of strings}
#' }
#'
#' @import data.table
#' @importFrom parallel detectCores
#' @importFrom stats p.adjust.methods
#'
parameters_are_good <- function(annot, count_data_A, count_data_B, boot_data_A, boot_data_B, TARGET_COL, PARENT_COL,
                                correction, testmode, scaling, threads, seed, p_thresh, abund_thresh, dprop_thresh, 
                                qboot, qbootnum, qrep_thresh, rboot, rrep_thresh, reckless, lean, verbose, use_sums) {
  warnmsg <- list()

  # Input format.
  if(any(is.null(annot),
         all( any(is.null(count_data_A), is.null(count_data_B)),
              any(is.null(boot_data_A), is.null(boot_data_B)) ) ))
    return(list("error"=TRUE, "message"="Insufficient data parameters!"))
  if(any(is.null(count_data_A) != is.null(count_data_B), is.null(boot_data_A) != is.null(boot_data_B)))
    return(list("error"=TRUE, "message"="You must specify (the same type of) data for both conditions!"))
  # Booted counts.
  if (!is.null(boot_data_A)) {
    if (any(!is.list(boot_data_A), !is.list(boot_data_A), !is.data.table(boot_data_A[[1]]), !is.data.table(boot_data_B[[1]]) ))
      return(list("error"=TRUE, "message"="The bootstrap data are not lists of data.tables!"))
  }
  # Counts.
  if (!is.null(count_data_A)) {  # If boot_data available, ignore count_data.
    if(is.null(boot_data_A)) {
      if (!is.data.table(count_data_A) | !is.data.table(count_data_B))
        return(list("error"=TRUE, "message"="The counts data are not data.tables!"))
    } else {
      warnmsg["cnts&boots"] <- "Received multiple input formats! Only the bootstrapped data will be used."
    }
  }
  
  # Annotation.
  if (!is.logical(reckless) || is.na(reckless))
    return(list("error"=TRUE, "message"="Invalid value to reckless!"))
  if (!is.data.frame(annot))
    return(list("error"=TRUE, "message"="The provided annot is not a data.frame!"))
  if (any(!c(TARGET_COL, PARENT_COL) %in% names(annot)))
    return(list("error"=TRUE, "message"="The specified target and/or parent IDs field names do not exist in annot!"))
  if (length(annot[[TARGET_COL]]) != length(unique(annot[[TARGET_COL]])))
    return(list("error"=TRUE, "message"="Some transcript identifiers are not unique in `annot`!"))
  
  # Simple parameters.
  if (!correction %in% stats::p.adjust.methods)
    return(list("error"=TRUE, "message"="Invalid p-value correction method name. Refer to stats::p.adjust.methods."))
  if (!testmode %in% c("genes", "transc", "both"))
    return(list("error"=TRUE, "message"="Unrecognized value for testmode!"))
  if ((!is.numeric(p_thresh)) || p_thresh > 1 || p_thresh < 0)
    return(list("error"=TRUE, "message"="Invalid p-value threshold!"))
  if ((!is.numeric(dprop_thresh)) || dprop_thresh < 0 || dprop_thresh > 1)
    return(list("error"=TRUE, "message"="Invalid proportion difference threshold! Must be between 0 and 1."))
  if ((!is.numeric(abund_thresh)) || abund_thresh < 0)
    return(list("error"=TRUE, "message"="Invalid abundance threshold! Must be a count >= 0."))
  if ((!is.numeric(qbootnum)) || qbootnum < 0)
    return(list("error"=TRUE, "message"="Invalid number of bootstraps! Must be a positive integer number."))
  if (!is.logical(qboot))
    return(list("error"=TRUE, "message"="Unrecognized value for qboot! Must be TRUE/FALSE."))
  if (!is.logical(rboot))
    return(list("error"=TRUE, "message"="Unrecognized value for rboot! Must be TRUE/FALSE."))
  if ((!is.numeric(qrep_thresh)) || qrep_thresh < 0 || qrep_thresh > 1)
    return(list("error"=TRUE, "message"="Invalid quantification reproducibility threshold! Must be between 0 and 1."))
  if ((!is.numeric(rrep_thresh)) || rrep_thresh < 0 || rrep_thresh > 1)
    return(list("error"=TRUE, "message"="Invalid cross-replicate reproducibility threshold! Must be between 0 and 1."))
  if ((!is.numeric(threads)) || threads <= 0)
    return(list("error"=TRUE, "message"="Invalid number of threads! Must be positive integer."))
  if (threads > parallel::detectCores(logical= TRUE))
    return(list("error"=TRUE, "message"="Number of threads exceeds system's reported capacity."))
  if (!is.logical(lean))
    return(list("error"=TRUE, "message"="Invalid value for lean! Must be TRUE/FALSE."))
  if (!is.logical(verbose))
    return(list("error"=TRUE, "message"="Invalid value for verbose! Must be TRUE/FALSE."))
  if (!is.logical(use_sums))
    return(list("error"=TRUE, "message"="Invalid value for use_sums! Must be TRUE/FALSE."))
  
  # Scaling
  nsmpl <- NULL
  if (!is.null(count_data_A)){
    nsmpl <- dim(count_data_A)[2] + dim(count_data_B)[2]
  } else {
    nsmpl <- length(boot_data_A) + length(boot_data_B)
  }
  if (!is.numeric(scaling) || any(scaling < 1) || (length(scaling) != 1 && length(scaling) != nsmpl))
    return(list("error"=TRUE, "message"="Invalid scaling factor(s)! Must be vector of non-zero numbers. Vector length must be either 1 or equal to the combined number of samples in the two conditions."))
  if (!is.na(seed) && (!is.numeric(seed) || as.integer(seed) != seed) )
    return(list("error"=TRUE, "message"="Invalid seed! Must be integer or NA_integer_."))
  
  # Bootstrap.
  minboots <- NA_integer_
  samples_by_condition <- NULL
  numsamples <- NA_integer_
  if (!is.null(boot_data_A)) {
  	numsamples <- length(boot_data_A) + length(boot_data_B)
  }
  maxmatrix <- 2^31 - 1
  if (qboot) {
    if (!is.null(boot_data_A)) {
      # Compared to annotation.
      if ( !all(boot_data_A[[1]][[1]] %in% annot[[TARGET_COL]]) && !all(annot[[TARGET_COL]] %in% boot_data_A[[1]][[1]]) ) {
        if(reckless){
          warnmsg["annotation-mismatch"] <- "The transcript IDs in the quantifications and the annotation do not match completely! Did you use different annotations? Reckless mode enabled, continuing at your own risk..."
        } else {
          return(list("error"=TRUE, "message"="The transcript IDs in the quantifications and the annotation do not match completely! Did you use different annotations?"))
        }
      }
      # Among samples
      tx <- boot_data_A[[1]][[1]][order(boot_data_A[[1]][[1]])]
      for (k in 2:length(boot_data_A)){
        if (!all( tx == boot_data_A[[k]][[1]][order(boot_data_A[[k]][[1]])] )) {
          if (reckless) {
            warnmsg[paste0("bootIDs-A", k)] <- paste0("Inconsistent set of transcript IDs across samples (A", k, ")! Did you use different annotations? Reckless mode enabled, continuing at your own risk...")
          } else {
            return(list("error"=TRUE, "message"="Inconsistent set of transcript IDs across samples! Did you use different annotations?"))
          }
        }
      }
      for (k in 1:length(boot_data_B)){
        if (!all( tx == boot_data_B[[k]][[1]][order(boot_data_B[[k]][[1]])] ))
          if (reckless) {
            warnmsg[paste0("bootIDs-B", k)] <- paste0("Inconsistent set of transcript IDs across samples (B", k, ")! Did you use different annotations? Reckless mode enabled, continuing at your own risk...")
          } else {
            return(list("error"=TRUE, "message"="Inconsistent set of transcript IDs across samples! Did you use different annotations?"))
          }
      }
    
      # Number of iterations.
      minboots <- infer_bootnum(boot_data_A, boot_data_B)
      if (!is.na(minboots)) {
        if (minboots <= 1)
          return(list("error"=TRUE, "message"="It appears some of your samples have no bootstraps!"))
        if (minboots < 100) {
          warnmsg["toofewboots"] <- "Your quantifications have few bootstrap iterations, which reduces reproducibility of the calls."
        }
        if (qbootnum < 100 && qbootnum != 0) {
          warnmsg["toolowbootnum"] <- "The requested qbootnum is low, which reduces reproducibility of the calls."
        }
        bootcombos <- minboots^numsamples  # Conservative estimate.
        if (qbootnum >= bootcombos/100)
          warnmsg["toomanyboots"] <- "The requested number of quantification bootstraps is very high, relatively to the supplied data. Over 1% chance of duplicate iterations."
        if (qbootnum >= maxmatrix/dim(annot)[1])
          return(list("error"=TRUE,"message"="The requested number of quantification bootstraps would exceed the maximum capacity of an R matrix."))
      } # else it is probably count data and qboot will be auto-set to FALSE
    } else {
      warnmsg["noboots"] <- "qboot is TRUE but no bootstrapped estimates were provided! Continuing without bootstrapping."
    }
  }
  if (rboot & any( length(samples_by_condition[[1]]) * length(samples_by_condition[[2]]) > maxmatrix/dim(annot)[1],
  				         length(boot_data_A)               * length(boot_data_B)               > maxmatrix/dim(annot)[1] ) )
    warnmsg["toomanyreplicates"] <- "The number of replicates is too high. Exhaustive 1vs1 would exceed maximum capacity of an R matrix."
  
  return(list("error"=FALSE, "message"="All good!", "maxboots"=minboots, "warn"=(length(warnmsg) > 0), "warnings"= warnmsg))
}


#================================================================================
#' Feature IDs match between 
#'
#' This is applied after bootstrapped counts are checked and combined down to a single vector of values per sample.
#' So it covers both bootstrapped and non-bootstrapped input.
#' 
#' @param annot Annotation dataframe.
#' @param tids List of transcript ID vectors.
#' @return bool
#' 
annotation_inconsistent <- function(annot, tids) {
  return( any( vapply(tids, function(v) { return( !all(annot$target_id %in% v) | !all(v %in% annot$target_id) ) }, logical(1)) ) )
}

#================================================================================
#' Rule-of-thumb number of iterations.
#'
#' @param boot_data_A List of tables of bootstrapped counts.
#' @param boot_data_B List of tables of bootstrapped counts.
#' @return The least number of iterations seen in the input.
#'
infer_bootnum <- function(boot_data_A, boot_data_B){
  minboot <- NA_integer_
  minboot <- length(boot_data_A[[1]]) - 1   # ID column
  for (k in 2:length(boot_data_A)){
    minboot <- min(minboot, length(boot_data_A[[k]]) - 1)
  }
  for (k in 1:length(boot_data_B)){
    minboot <- min(minboot, length(boot_data_B[[k]]) - 1)
  }
  return(minboot)
}


#================================================================================
#' Structure input data.
#' 
#' @param tx_filter Pre-processed annotation.
#' @param steps Infered input type.
#' @param threads Threads.
#' @param count_data_A Non-bootstrapped data for one condition.
#' @param count_data_B Non-bootstrapped data for other condition.
#' @param boot_data_A Boostrapped data for one condition.
#' @param bootsB Bootstrapped data for other condition.
#' @return List.
#'   
#' @import data.table
#' @importFrom parallel mclapply
#'
structure_data <- function(tx_filter, steps, threads, count_data_A, count_data_B, boot_data_A, boot_data_B) {
  # Preprocess bootstrapped data.
  tids <- list()
  if (steps == 2) {
    # Re-order rows to match the annotation.
    boot_data_A <- lapply(boot_data_A, function(x) { x[match(tx_filter$target_id, x[[1]]), ] })
    boot_data_B <- lapply(boot_data_B, function(x) { x[match(tx_filter$target_id, x[[1]]), ] })
    # Remove ID columns so I don't have to always subset for math operations.
    for (smpl in c(boot_data_A, boot_data_B))
      smpl[, 1 := NULL]
    # Calculate mean count across bootstraps.
    count_data_A <- as.data.table(mclapply(boot_data_A, function(b) { return(rowMeans(b)) },
                                           mc.cores = threads, mc.preschedule = TRUE, mc.allow.recursive = FALSE))
    count_data_B <- as.data.table(mclapply(boot_data_B, function(b) { return(rowMeans(b)) },
                                           mc.cores = threads, mc.preschedule = TRUE, mc.allow.recursive = FALSE))
    # Preprocess unbootstrapped/debootstrapped data.
  } else {
    # Just re-order rows and crop out the transcript IDs.
    tids[["a"]] <- count_data_A[[1]]
    nn <- names(count_data_A)
    count_data_A <- count_data_A[match(tx_filter$target_id, count_data_A[[1]]),  nn[seq.int(2, length(nn))],  with=FALSE]
    tids[["b"]] <- count_data_B[[1]]
    nn <- names(count_data_B)  # The number of columns may be different from A.
    count_data_B <- count_data_B[match(tx_filter$target_id, count_data_B[[1]]),  nn[seq.int(2, length(nn))],  with=FALSE]
  }
  
  return (list("ca"=count_data_A, "cb"=count_data_B, "ba"=boot_data_A, "bb"=boot_data_B, "ids"=tids))
}


#================================================================================
#' Apply scaling.
#' 
#' @param scaling Scaling factor(s).
#' @param threads Threads.
#' @param count_data_A Non-bootstrapped data for one condition.
#' @param count_data_B Non-bootstrapped data for other condition.
#' @param boot_data_A Boostrapped data for one condition.
#' @param bootsB Bootstrapped data for other condition.
#' @return List.
#'   
#' @import data.table
#' @importFrom parallel mclapply
#'
scale_data <- function(scaling, threads, steps, count_data_A, count_data_B, boot_data_A, boot_data_B) {
  # Scaling factors.
  lA <- dim(count_data_A)[2]
  lB <- dim(count_data_B)[2]
  sfA <- NULL
  sfB <- NULL
  if (length(scaling) == 1) {
    sfA <- rep(scaling, lA)
    sfB <- rep(scaling, lB)
  } else {
    sfA <- scaling[1:lA]
    sfB <- scaling[(1+lA):(lA+lB)]
  }
  # Scale counts.
  count_data_A <- as.data.table( mclapply(1:lA, function(x) { 
    return(count_data_A[[x]] * sfA[x])  # Can't apply scaling to whole table in one step, because each column/sample can have a different scaling factor.
  }, mc.cores=threads, mc.allow.recursive = FALSE) )
  count_data_B <- as.data.table( mclapply(1:lB, function(x) { 
    return(count_data_B[[x]] * sfB[x]) 
  }, mc.cores=threads, mc.allow.recursive = FALSE) )
  # Also scale the bootstraps, if they exist.
  if (steps > 1){
    boot_data_A <- lapply(1:lA , function (smpl){
      return(as.data.table( mclapply(boot_data_A[[smpl]], function(b) { 
        return(b * sfA[smpl]) # The bootstraps table belongs to a single sample, so here I can apply the factor to the whole table.
      }, mc.cores=threads, mc.allow.recursive = FALSE) )) })
    boot_data_B <- lapply(1:lB , function (smpl){
      return(as.data.table( mclapply(boot_data_B[[smpl]], function(b) { 
        return(b * sfB[smpl])
      }, mc.cores=threads, mc.allow.recursive = FALSE) )) })
  }
  
  return (list("ca"=count_data_A, "cb"=count_data_B, "ba"=boot_data_A, "bb"=boot_data_B))
}


#================================================================================
#' Create output structure.
#'
#' @param annot Pre-processed by \code{tidy_annot()}.
#' @param full Full-sized structure or core fields only. Either "full", "short" or "partial".
#' @param n How many scaling factors to reserve space for.
#' @return A list.
#'
#' @import data.table
#'
alloc_out <- function(annot, full, n=1){
  # Ensure data.table complies.
  # if (packageVersion("data.table") >= "1.9.8")
    setDTthreads(1)
  if (full == "full") {
    Parameters <- list("description"=NA_character_, "time"=date(),
                       "rats_version"=packageVersion("rats"), "R_version"=R.Version()[c("platform", "version.string")],
                       "var_name"=NA_character_, "cond_A"=NA_character_, "cond_B"=NA_character_,
                       "data_type"=NA_character_, "num_replic_A"=NA_integer_, "num_replic_B"=NA_integer_,
                       "num_genes"=NA_integer_, "num_transc"=NA_integer_,
                       "tests"=NA_character_, "use_sums"=NA, "correction"=NA_character_, 
                       "p_thresh"=NA_real_, "abund_thresh"=NA_real_, "dprop_thresh"=NA_real_, "abund_scaling"=numeric(length=n),
                       "quant_boot"=NA,"quant_reprod_thresh"=NA_real_,  "quant_bootnum"=NA_integer_,
                       "rep_boot"=NA, "rep_reprod_thresh"=NA_real_, "rep_bootnum"=NA_integer_, 
                       "seed"=NA_integer_, "reckless"=NA, "lean"=NA)
    Genes <- data.table("parent_id"=as.vector(unique(annot$parent_id)),
                        "elig"=FALSE, "sig"=NA, "elig_fx"=NA, "quant_reprod"=NA, "rep_reprod"=NA, "DTU"=FALSE, "transc_DTU"=NA,
                        "known_transc"=NA_integer_, "detect_transc"=NA_integer_, "elig_transc"=NA_integer_, "maxDprop"=NA_real_,
                        "pval"=NA_real_, "pval_corr"=NA_real_,
                        "quant_p_median"=NA_real_, "quant_p_min"=NA_real_, "quant_p_max"=NA_real_,
                        "quant_na_freq"=NA_real_, "quant_dtu_freq"=NA_real_,
                        "rep_p_median"=NA_real_, "rep_p_min"=NA_real_, "rep_p_max"=NA_real_,
                        "rep_na_freq"=NA_real_, "rep_dtu_freq"=NA_real_)
    with(Genes,
         setkey(Genes, parent_id) )
    Transcripts <- data.table("target_id"=annot$target_id, "parent_id"=annot$parent_id,
                              "elig_xp"=NA, "elig"=FALSE, "sig"=NA, "elig_fx"=NA, "quant_reprod"=NA, "rep_reprod"=NA, "DTU"=FALSE, "gene_DTU"=NA,
                              "abundA"=NA_real_, "abundB"=NA_real_, "stdevA"=NA_real_, "stdevB"=NA_real_, "log2FC"=NA_real_, 
                              "totalA"=NA_real_, "totalB"=NA_real_, "propA"=NA_real_, "propB"=NA_real_, "Dprop"=NA_real_,
                              "pval"=NA_real_, "pval_corr"=NA_real_,
                              "quant_p_median"=NA_real_, "quant_p_min"=NA_real_, "quant_p_max"=NA_real_,
                              "quant_Dprop_mean"=NA_real_, "quant_Dprop_stdev"=NA_real_, "quant_Dprop_min"=NA_real_, "quant_Dprop_max"=NA_real_,
                              "quant_na_freq"=NA_real_, "quant_dtu_freq"=NA_real_,
                              "rep_p_median"=NA_real_, "rep_p_min"=NA_real_, "rep_p_max"=NA_real_,
                              "rep_Dprop_mean"=NA_real_, "rep_Dprop_stdev"=NA_real_, "rep_Dprop_min"=NA_real_, "rep_Dprop_max"=NA_real_,
                              "rep_na_freq"=NA_real_, "rep_dtu_freq"=NA_real_)
    with(Transcripts,
         setkey(Transcripts, parent_id, target_id) )
    CountData <- list()
  } else if (full == "short") {
    Parameters <- list("num_replic_A"=NA_integer_, "num_replic_B"=NA_integer_)
    Genes <- data.table("parent_id"=levels(as.factor(annot$parent_id)),
                        "elig_transc"=NA_integer_, "elig"=FALSE, "elig_fx"=NA,
                        "pval"=NA_real_, "pval_corr"=NA_real_, "sig"=NA, "DTU"=FALSE)
    with(Genes,
         setkey(Genes, parent_id) )
    Transcripts <- data.table("target_id"=annot$target_id, "parent_id"=annot$parent_id,
                              "abundA"=NA_real_, "abundB"=NA_real_,
                              "log2FC"=NA_real_,  #log2FC currently not used in decision making and bootstrapping, but maybe in the future.
                              "totalA"=NA_real_, "totalB"=NA_real_, "elig_xp"=NA, "elig"=FALSE,
                              "propA"=NA_real_, "propB"=NA_real_, "Dprop"=NA_real_, "elig_fx"=NA,
                              "pval"=NA_real_, "pval_corr"=NA_real_, "sig"=NA, "DTU"=FALSE)
    with(Transcripts,
         setkey(Transcripts, parent_id, target_id) )
    CountData <- NULL
  } else {  # Partial.
    Parameters <- NULL
    Genes <- data.table("parent_id"=as.vector(unique(annot$parent_id)),
                        "p_median"=NA_real_, "p_min"=NA_real_, "p_max"=NA_real_,
                        "na_freq"=NA_real_, "dtu_freq"=NA_real_)
    Transcripts <- data.table("target_id"=annot$target_id, "parent_id"=annot$parent_id,
                              "p_median"=NA_real_, "p_min"=NA_real_, "p_max"=NA_real_,
                              "Dprop_mean"=NA_real_, "Dprop_stdev"=NA_real_, "Dprop_min"=NA_real_, "Dprop_max"=NA_real_,
                              "na_freq"=NA_real_, "dtu_freq"=NA_real_)
    CountData <- NULL
  }

  return(list("Parameters"= Parameters, "Genes"= Genes, "Transcripts"= Transcripts, "Abundances"= CountData))
}


#================================================================================
#' Set up and execute the tests.
#'
#' @param counts_A A data.table of counts for condition A. x: sample, y: transcript.
#' @param counts_B A data.table of counts for condition B. x: sample, y: transcript.
#' @param tx_filter A data.table with target_id and parent_id. Pre-processed with \code{tidy_annot()}.
#' @param test_transc Whether to do transcript-level test.
#' @param test_genes Whether to do gene-level test.
#' @param full Either "full" (for complete output structure) or "short" (for bootstrapping).
#' @param count_thresh Minimum average count across replicates.
#' @param p_thresh The p-value threshold.
#' @param dprop_thresh Minimum difference in proportions.
#' @param correction Multiple testing correction type.
#' @param threads Number of threads (POSIX systems only).
#' @param use_sums Use sums instead of means. (sums is the behaviour of RATs up to 0.6.5 inclusive, more sensitivie but more prone to false positives)
#' @return list
#'
#' @import data.table
#' @import utils
#' @import parallel
#'
calculate_DTU <- function(counts_A, counts_B, tx_filter, test_transc, test_genes, full, count_thresh, p_thresh, dprop_thresh, correction, threads= 1, use_sums=FALSE) {
  # Ensure data.table complies.
  setDTthreads(threads)

  #---------- PRE-ALLOCATE

  # Pre-allocate results object.
  resobj <- alloc_out(tx_filter, full, dim(counts_A)[2] + dim(counts_B)[2])
  resobj$Parameters["num_replic_A"] <- dim(counts_A)[2]
  resobj$Parameters["num_replic_B"] <- dim(counts_B)[2]

  with(resobj, {
    # Set key to gene ids.
    setkey(Transcripts, parent_id)
    setkey(Genes, parent_id)

    #---------- STATS

    # Isoforms abundance combined across replicates.
    if (use_sums){
      Transcripts[, abundA :=  rowSums(counts_A, na.rm=TRUE) ]
      Transcripts[, abundB :=  rowSums(counts_B, na.rm=TRUE) ]
    } else {
      Transcripts[, abundA :=  rowSums(counts_A, na.rm=TRUE) ]
      Transcripts[, abundB :=  rowSums(counts_B, na.rm=TRUE) ]
    }
    # Gene expression
    Transcripts[, totalA := sum(abundA, na.rm=TRUE), by=parent_id]
    Transcripts[, totalB := sum(abundB, na.rm=TRUE), by=parent_id]
    # Relative isoform abundance
    Transcripts[, propA := abundA/totalA]
    Transcripts[, propB := abundB/totalB]
    Transcripts[(is.nan(propA)), propA := NA_real_]  # Replace NaN with NA.
    Transcripts[(is.nan(propB)), propB := NA_real_]
    Transcripts[, Dprop := propB - propA]

    #---------- FILTER

    # Filter transcripts and genes to reduce number of tests:
    if (use_sums){
      ctA <- count_thresh * resobj$Parameters[["num_replic_A"]]  # Adjust count threshold for number of replicates.
      ctB <- count_thresh * resobj$Parameters[["num_replic_B"]]
    } else {
      ctA <- count_thresh
      ctB <- count_thresh
    }
    
    Transcripts[, elig_xp := (abundA >= ctA | abundB >= ctB)]
    Transcripts[, elig := (elig_xp &                          # isoform expressed above the noise threshold in EITHER condition
                             totalA >= ctA & totalB >= ctB &  # at least one eligibly expressed isoform exists in EACH condition (prevent gene-on/off from creating wild proportions from low counts)
                             totalA != 0 & totalB != 0 &      # gene expressed in BOTH conditions (failsafe, in case user sets noise threshold to 0)
                             (propA != 1 | propB != 1) )]     # not the only isoform expressed.
    
    Genes[, elig_transc := Transcripts[, as.integer(sum(elig, na.rm=TRUE)), by=parent_id][, V1] ]
    Genes[, elig := elig_transc >= 2]     # Need at least two expressed isoforms in order for DTU to even be possible.

    # Biologically significant.
    Transcripts[, elig_fx := abs(Dprop) >= dprop_thresh]
    Genes[, elig_fx := Transcripts[, any(elig_fx), by = parent_id][, V1] ]

    #---------- TESTS

    
    # Transcript-level test.
    if (test_transc) {
      if(any(Transcripts$elig)){  # If nothing is eligible, table subsetting by `elig` is nonsense and crashes. 
                                  # In that case we can just skip testing altegether, as all output fields are initialized with NA already.
      Transcripts[(elig), pval := unlist( mclapply( as.data.frame(t(Transcripts[(elig), .(abundA, abundB, totalA, totalB)])),
                                                      function(row) {
                                                        return( g.test.2(obsx= c(row[1], row[3]-row[1]), obsy= c(row[2], row[4]-row[2])) )
                                                      }, mc.cores= threads, mc.allow.recursive= FALSE, mc.preschedule= TRUE)
                                            ) ]
        Transcripts[(elig), pval_corr := p.adjust(pval, method=correction)]
        Transcripts[(elig), sig := pval_corr < p_thresh]
        Transcripts[(elig), DTU := sig & elig_fx]
      }
    }

    # Gene-level test.
    if (test_genes) {
      if(any(Genes$elig)){  # Skip testing if nothing is eligible, otherwise table subsetting by `elig` is nonsense and crashes. 
        Genes[(elig), pval := t( as.data.frame( mclapply(Genes[(elig), parent_id],
                                          function(parent) {
                                              # Extract all relevant data to avoid repeated look ups in the large table.
                                              subdt <- Transcripts[parent, .(abundA, abundB)]  # All isoforms, including not detected ones.
                                              return( g.test.2(obsx= subdt$abundA, obsy= subdt$abundB) )
                                          }, mc.cores= threads, mc.preschedule= TRUE, mc.allow.recursive= FALSE)
                  )) ]
        Genes[(elig), pval_corr := p.adjust(pval, method=correction)]
        Genes[(elig), sig := pval_corr < p_thresh]
        Genes[(elig), DTU := sig & elig_fx]
      }
    }
  })
  return(resobj)
}


#================================================================================
#' Perform the bootstraps.
#'
#' Abstracts the bootstrap process from the specific field names of each bootstrap
#' type and removes the associated code duplication.
#' 
#' Providing bootstrapped data will bootstrap the quantifications, whereas plain
#' data will bootstrap the replicates. Do not provide both.
#' 
#' @param lean Minimal insight reported.
#' @param tx_filter The annotation.
#' @param eltr Selection vector for eligible transcripts.
#' @param elge Selection vector for eligible genes.
#' @param boot_data_A Bootstrapped quantifications for first condition.
#' @param boot_data_B Bootstrapped quantifications for second condition.
#' @param count_data_A Plain quantifications for first condition.
#' @param count_data_B Plain quantifications for second condition.
#' @param niter Number of iterations (for quantification bootstraps only).
#' @param threads Number of threads.
#' @param verbose Progress messages.
#' @param test_transc Do transcript-level test.
#' @param test_genes Do gene-level test.
#' @param ... Other arguments to pass on to calculate_DTU().
#' 
#' @import utils
#' @import parallel
#' @import matrixStats
#' @import data.table
#
do_boot <- function(lean, tx_filter, eltr, elge,
                    boot_data_A=NULL, boot_data_B=NULL, count_data_A=NULL, count_data_B=NULL, 
                    niter=NULL, threads=1, verbose=TRUE, test_transc=TRUE, test_genes=TRUE, ...) 
{
  
  if (!is.null(count_data_A)) {
    pairs <- as.data.frame(t( expand.grid(1:dim(count_data_A)[2], 1:dim(count_data_B)[2]) ))
    niter <- length(pairs)
  }
  
  if (verbose)
    myprogress <- txtProgressBar(min = 0, max = niter, initial = 0, char = "=", width = NA, style = 3, file = stderr())
  
  # Structure to collect the results.
  tallies <- alloc_out(tx_filter, "partial")
  with(tallies, {
    Genes[, dtu_freq := 0]
    Genes[, na_freq := 0]
    Transcripts[, dtu_freq := 0]
    Transcripts[, na_freq := 0]
  })
  
  bootout <- lapply(1:niter, function(i) {  # Single-threaded. Forking happens within calculate_DTU(). This allows call_DTU() to track and display progress.
    if (verbose)
      setTxtProgressBar(myprogress, i)
    
    if(!is.null(boot_data_A)) {
      # Grab a random bootstrap from each replicate.
      counts_A <- as.data.table(lapply(boot_data_A, function(smpl) { smpl[[sample( names(smpl)[1:ncol(smpl)], 1)]] }))
      counts_B <- as.data.table(lapply(boot_data_B, function(smpl) { smpl[[sample( names(smpl)[1:ncol(smpl)], 1)]] }))
      # I have to use list syntax to get a vector back. Usual table syntax with "with=FALSE" returns a table and I fail to cast it.
    } else {
      # Grab the next replicate pair, one from each condition.
      # Scale them up for the number of respective samples, ie. "what if all my samples in this condition were very similar to the selected one".
      counts_A <- as.data.table( count_data_A[[ names(count_data_A)[pairs[[i]][1]] ]] ) * ncol(count_data_A)
      counts_B <- as.data.table( count_data_B[[ names(count_data_B)[pairs[[i]][2]] ]] ) * ncol(count_data_B)
    }
    
    # Test.
    iterout <- calculate_DTU(counts_A, counts_B, tx_filter, test_transc, test_genes, "short", ...)
    
    if (lean){
      # Calculate DTU rate in-place.
      with(tallies, {
        if (test_genes) {
          Genes[, dtu_freq := (dtu_freq * (i-1) + ifelse(iterout$Genes$DTU, 1, 0)) / i]  # incremental mean for i values based on the ith value and the mean of the other (i-1) values.
          Genes[, na_freq := (na_freq * (i-1) + ifelse(iterout$Genes$elig, 0, 1)) / i]
        }
        if (test_transc) {
          Transcripts[, dtu_freq := (dtu_freq * (i-1) + ifelse(iterout$Transcripts$DTU, 1, 0)) / i]
          Transcripts[, na_freq := (na_freq * (i-1) + ifelse(iterout$Transcripts$elig, 0, 1)) / i]
        }
      })
      return(NULL)  # The values of interest are updated directly in every iteration, nothing to return.
    } else {
      # Keep multiple columns of info per iteration so as to calculate stuff once bootstrapping is completed.
      return( list("tp" = iterout$Transcripts[, "pval_corr", with=FALSE],
                   "tdtu" = iterout$Transcripts[, "DTU", with=FALSE],
                   "ty" = iterout$Transcripts[, "Dprop", with=FALSE],
                   "gp" = iterout$Genes[, "pval_corr", with=FALSE],
                   "gdtu" = iterout$Genes[, "DTU", with=FALSE]) )
    }
  }) # End of iterations loop.
  
  if (verbose)
    message("Summarising bootstraps...")
  
  with(tallies, {
    if (test_transc) {
      if (!lean) {
        # Collate and summarise the like columns across iterations.
        td <- as.matrix(as.data.table(mclapply(bootout, function(iter) { iter[["tdtu"]] }, 
                                               mc.cores= threads, mc.preschedule = TRUE, mc.allow.recursive = FALSE)))
        Transcripts[(eltr), dtu_freq := rowCounts(td[eltr, ], value = TRUE, na.rm=TRUE) / niter]
        
        tp <- as.matrix(as.data.table(mclapply(bootout, function(iter) { iter[["tp"]] }, 
                                               mc.cores= threads, mc.preschedule = TRUE, mc.allow.recursive = FALSE)))
        Transcripts[(eltr), p_median := rowMedians(tp[eltr, ], na.rm = TRUE)]
        Transcripts[(eltr), p_min := rowMins(tp[eltr, ], na.rm = TRUE)]
        Transcripts[(eltr), p_max := rowMaxs(tp[eltr, ], na.rm = TRUE)]
        Transcripts[(eltr), na_freq := rowCounts(tp[eltr, ], value = NA, na.rm=FALSE) / niter]
        
        ty <- as.matrix(as.data.table(mclapply(bootout, function(iter) { iter[["ty"]] }, 
                                               mc.cores= threads, mc.preschedule = TRUE, mc.allow.recursive = FALSE)))
        Transcripts[(eltr), Dprop_mean := rowMeans(ty[eltr, ], na.rm = TRUE)]
        Transcripts[(eltr), Dprop_stdev := rowSds(ty[eltr, ], na.rm = TRUE)]
        Transcripts[(eltr), Dprop_min := rowMins(ty[eltr, ], na.rm = TRUE)]
        Transcripts[(eltr), Dprop_max := rowMaxs(ty[eltr, ], na.rm = TRUE)]
      }
    }
    
    if (test_genes) {
      if (!lean) {
        # Collate and summarise the like columns across iterations.
        gd <- as.matrix(as.data.table(mclapply(bootout, function(iter) { iter[["gdtu"]] }, 
                                               mc.cores= threads, mc.preschedule = TRUE, mc.allow.recursive = FALSE)))
        Genes[(elge), dtu_freq := rowCounts(gd[elge, ], value = TRUE, na.rm = TRUE) / niter]
        
        gp <- as.matrix(as.data.table(mclapply(bootout, function(iter) { iter[["gp"]] }, 
                                               mc.cores= threads, mc.preschedule = TRUE, mc.allow.recursive = FALSE)))
        Genes[(elge), p_median := rowMedians (gp[elge, ], na.rm = TRUE)]
        Genes[(elge), p_min := rowMins(gp[elge, ], na.rm = TRUE)]
        Genes[(elge), p_max := rowMaxs(gp[elge, ], na.rm = TRUE)]
        Genes[(elge), na_freq := rowCounts(gp[elge, ], value = NA, na.rm = FALSE) / niter]
      }
    }
  })
  
  return (tallies)
}


#================================================================================
#' Log-likelihood test of goodness of fit.
#'
#' For a set of observations against a set of probabilities.
#' 
#' General utility function.
#' 
#' @param obsx	a numeric vector of finite positive counts, with at least one non-zero value.
#' @param px	a vector of probabilities of the same length as xobs. ( sum(px) <= 1 )
#'
#' The order of values in the two vectors should be the same.
#' If any value of xp is zero and the corresponding xobs is not, g.test.1 will always reject the hypothesis.
#' No corrections are applied.
#' No input checks are applied.
#'
#' @importFrom stats pchisq
#' @export
#
g.test.1 <- function(obsx, px) {
  n = length(obsx)
  sx <- sum(obsx)
  ex <- sx * px  # expected values
  G <- 2 * sum( sapply (seq.int(1, n, 1),
                        function (i) { if (obsx[i] != 0) { return(obsx[i] * log(obsx[i]/ex[i])) } else { return(0) } }) )
  return( pchisq(G, df= n - 1, lower.tail= FALSE) )
}

#================================================================================
#' Log-likelihood test of independence.
#'
#' For two sets of observations.
#' 
#' General utility function.
#'
#' @param obsx	a numeric vector of positive counts, with at least one non-zero value.
#' @param obsy	a numeric vector of positive counts of same length as obsx, with at least one non-zero value.
#'
#' The order of values in the two vectors should be the same.
#' No corrections are applied.
#' No input checks are applied, as RATs needs to run this millions of times.
#'
#' @importFrom stats pchisq
#' @export
#
g.test.2 <- function(obsx, obsy) {
  n = length(obsx)
  # Row and column sums.
  sx <- sum(obsx)
  sy <- sum(obsy)
  sv <- obsx + obsy
  st <- sx + sy
  # Marginal probabilities.
  mpx <- sx / st
  mpy <- sy / st
  mpv <- sapply(sv, function(v) {v/st})
  # Expected values.
  ex <- mpx * mpv * st
  ey <- mpy * mpv * st
  # Statistic.
  G <- 2 * sum( sum( sapply (seq.int(1, n, 1),
                             function (i) { if (obsx[i] != 0) { return(obsx[i] * log(obsx[i]/ex[i])) } else { return(0) } }) ),
                sum( sapply (seq.int(1, n, 1),
                             function (i) { if (obsy[i] != 0) { return(obsy[i] * log(obsy[i]/ey[i])) } else { return(0) } }) )
              )
  return( pchisq(G, df= n - 1, lower.tail= FALSE) )
}


#================================================================================
#' Tidy up annotation.
#' 
#' With the general english meaning of the word "tidy", not the R Tidyverse meaning.
#' 
#' @param annot A data.frame matching transcript identifiers to gene identifiers. Any additional columns are allowed but ignored.
#' @param TARGET_COL The name of the column for the transcript identifiers in \code{annot}. (Default \code{"target_id"})
#' @param PARENT_COL The name of the column for the gene identifiers in \code{annot}. (Default \code{"parent_id"})
#' @return A data.table with genes under \code{parent_id} and transcripts under \code{target_id}, regardless of what the input annotation named them. Rows ordered first by gene and then by transcript.
#'
#' Extract the relevant columns as a \code{data.table} and order the rows by gene. This creates a fast
#' look-up structure that is consistent throughout RATs.
#'
#' @import data.table
#' @export
#
tidy_annot <- function(annot, TARGET_COL="target_id", PARENT_COL="parent_id") {
  tx_filter <- data.table(target_id= annot[[TARGET_COL]], parent_id= annot[[PARENT_COL]])
  with( tx_filter,
        setkey(tx_filter, parent_id, target_id) )

  return(tx_filter)
}


#================================================================================
#' Get largest value by absolute comparison.
#'
#' General utility function.
#'
#' Returns the original signed value with the largest absolute value in the vector. 
#' In case of equal absolutes between a positive and negative value, the positive 
#' value will always be the one that is returned.
#'
#' @param v A numeric vector.
#' @return A numeric value.
#' 
#' @export
maxabs <- function(v) {
	if (all(is.na(v)))
		return(NA_real_)
	x <- max(v, na.rm=TRUE)
	n <- min(v, na.rm=TRUE)
	if (abs(x) >= abs(n)) {
		return(x)
	} else {
		return(n)
	}
}


#================================================================================
#' Group samples by covariate value.
#'
#' *Legacy function* No longer in use, but kept for possible general utility.
#' 
#' Converts a table where each row is a sample and each column a variable into a
#' lookup structure that matches each value of each variable to a vector of corresponding samples.
#'
#' @param covariates A dataframe with different factor variables. Like the \code{sample_to_covariates} table of a \code{sleuth} object. Each row is a sample, each column is a covariate, each cell is a covariate value for the sample.
#' @return list of lists (per covariate) of vectors (per factor level).
#'
#' @export
#'
group_samples <- function(covariates) {
  samplesByVariable <- list()
  for(varname in names(covariates)) {
    categories <- unique(covariates[[varname]])
    samplesByVariable[[varname]] <- list()
    for (x in categories) {
      samplesByVariable[[varname]][[x]] <- which(covariates[[varname]] == x)
    }
  }
  return(samplesByVariable)
}


#EOF
