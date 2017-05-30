#================================================================================
#' Check input parameters.
#'
#' @param slo Sleuth object
#' @param annot Annotation dataframe.
#' @param name_A Condition name.
#' @param name_B Condition name.
#' @param varname Name of condition variable.
#' @param COUNTS_COL Name of counts column in bootstrap.
#' @param correction P-value correction method.
#' @param p_thresh Significance level.
#' @param TARGET_COL Name of transcript id column in annotation.
#' @param PARENT_COL Name of gene id column in annotation.
#' @param BS_TARGET_COL Name of transcript id column in bootstrap.
#' @param count_thresh Minimum transcript abundance per sample.
#' @param testmode Which test to run.
#' @param qboot Whether to bootstrap against quantifications.
#' @param qbootnum Number of bootstrap iterations.
#' @param qrep_thresh Confidence threshold.
#' @param dprop_thresh Minimum change in proportion.
#' @param count_data_A A dataframe of estimated counts.
#' @param count_data_B A dataframe of estimated counts.
#' @param boot_data_A A list of dataframes, one per sample, each with all the bootstrapped estimetes for the sample.
#' @param boot_data_B A list of dataframes, one per sample, each with all the bootstrapped estimetes for the sample.
#' @param rboot Whether to bootstrap against samples.
#' @param rrep_thresh Confidence threshold.
#' @param rrep_as_crit Whether to use rrep as a DTU criterion.
#' @param threads Number of threads.
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
#'
parameters_are_good <- function(slo, annot, name_A, name_B, varname, COUNTS_COL,
                            correction, p_thresh, TARGET_COL, PARENT_COL, BS_TARGET_COL,
                            count_thresh, testmode, qboot, qbootnum, dprop_thresh,
                            count_data_A, count_data_B, boot_data_A, boot_data_B, qrep_thresh,
                            threads, rboot, rrep_thresh, rrep_as_crit) {
  warnmsg <- list()

  # Input format.
  if(any(is.null(annot),
         all( any(is.null(slo), is.null(name_A), is.null(name_B), is.null(varname), is.null(BS_TARGET_COL), is.null(COUNTS_COL)),
              any(is.null(count_data_A), is.null(count_data_B)),
              any(is.null(boot_data_A), is.null(boot_data_B)) ) ))
    return(list("error"=TRUE, "message"="Insufficient parameters!"))

  # Annotation.
  if (!is.data.frame(annot))
    return(list("error"=TRUE, "message"="The provided annot is not a data.frame!"))
  if (any(!c(TARGET_COL, PARENT_COL) %in% names(annot)))
    return(list("error"=TRUE, "message"="The specified target and/or parent IDs field names do not exist in annot!"))
  if (length(annot$target_id) != length(unique(annot$target_id)))
    return(list("error"=TRUE, "message"="Some transcript identifiers are not unique!"))

  # Parameters.
  if (!correction %in% p.adjust.methods)
    return(list("error"=TRUE, "message"="Invalid p-value correction method name. Refer to stats::p.adjust.methods ."))
  if (!testmode %in% c("genes", "transc", "both"))
    return(list("error"=TRUE, "message"="Unrecognized value for testmode!"))
  if ((!is.numeric(p_thresh)) || p_thresh > 1 || p_thresh < 0)
    return(list("error"=TRUE, "message"="Invalid p-value threshold!"))
  if ((!is.numeric(dprop_thresh)) || dprop_thresh < 0 || dprop_thresh > 1)
    return(list("error"=TRUE, "message"="Invalid proportion difference threshold! Must be between 0 and 1."))
  if ((!is.numeric(count_thresh)) || count_thresh < 0)
    return(list("error"=TRUE, "message"="Invalid abundance threshold! Must be between 0 and 1."))
  if ((!is.numeric(qbootnum)) || qbootnum < 0)
    return(list("error"=TRUE, "message"="Invalid number of bootstraps! Must be a positive integer number."))
  if (!is.logical(qboot))
    return(list("error"=TRUE, "message"="Unrecognized value for qboot! Must be TRUE/FALSE."))
  if (!is.logical(rboot))
    return(list("error"=TRUE, "message"="Unrecognized value for rboot! Must be TRUE/FALSE."))
  if ((!is.numeric(qrep_thresh)) || qrep_thresh < 0 || qrep_thresh > 1)
    return(list("error"=TRUE, "message"="Invalid reproducibility threshold! Must be between 0 and 1."))
  if ((!is.numeric(rrep_thresh)) || rrep_thresh < 0 || rrep_thresh > 1)
    return(list("error"=TRUE, "message"="Invalid reproducibility threshold! Must be between 0 and 1."))
  if ((!is.numeric(threads)) || threads <= 0)
    return(list("error"=TRUE, "message"="Invalid number of threads! Must be positive integer."))
  if (threads > parallel::detectCores(logical= TRUE))
    return(list("error"=TRUE, "message"="Number of threads exceeds system's reported capacity."))
  if (!is.logical(rrep_as_crit))
    return(list("error"=TRUE, "message"="Unrecognized value for rrep_as_crit! Must be TRUE/FALSE."))

  # Sleuth
  if (!is.null(slo)) {
    if (any(! c("kal","sample_to_covariates") %in% names(slo), ! "bootstrap" %in% names(slo$kal[[1]]) ))
      return(list("error"=TRUE, "message"="The specified sleuth object is not valid!"))
    if (!COUNTS_COL %in% names(slo$kal[[1]]$bootstrap[[1]]))
      return(list("error"=TRUE, "message"="The specified counts field name does not exist!"))
    if (!varname %in% names(slo$sample_to_covariates))
      return(list("error"=TRUE, "message"="The specified covariate name does not exist!"))
    if (any(!c(name_A, name_B) %in% slo$sample_to_covariates[[varname]] ))
      return(list("error"=TRUE, "message"="One or both of the specified conditions do not exist!"))
    if (!BS_TARGET_COL %in% names(slo$kal[[1]]$bootstrap[[1]]))
      return(list("error"=TRUE, "message"="The specified target IDs field name does not exist in the bootstraps!"))
  }

  # Booted counts.
  if (!is.null(boot_data_A)) {  # If slo available ignore boot_data.
    if(is.null(slo)) {
      if (any(!is.list(boot_data_A), !is.list(boot_data_A), !is.data.table(boot_data_A[[1]]), !is.data.table(boot_data_B[[1]]) ))
        return(list("error"=TRUE, "message"="The bootstrap data are not lists of data.tables!"))
    } else {
      warnmsg["slo&boots"] <- "Received multiple input formats! Only the sleuth object will be used."
    }
  }

  # Counts.
  if (!is.null(count_data_A)) {  # If slo or boot_data available, ignore count_data.
    if(is.null(slo) && any(is.null(boot_data_A), is.null(boot_data_B))) {
      if (any(!is.data.table(count_data_A), !is.data.table(count_data_B)))
        return(list("error"=TRUE, "message"="The counts data are not data.tables!"))
      if (!all( count_data_A[[1]][order(count_data_A[[1]])] == count_data_B[[1]][order(count_data_B[[1]])] ))
        return(list("error"=TRUE, "message"="Inconsistent set of transcript IDs between conditions!"))
    } else {
      warnmsg["cnts&boots"] <- "Received multiple input formats! Only the bootstrapped data will be used."
    }
  }

  # Bootstrap.
  minboots <- NA_integer_
  samples_by_condition <- NULL
  numsamples <- NA_integer_
  if (!is.null(slo)) {
  	samples_by_condition <- group_samples(slo$sample_to_covariates)[[varname]]
  	numsamples <- length(samples_by_condition[[1]]) + length(samples_by_condition[[2]])
  } else if (!is.null(boot_data_A) && !is.null(boot_data_B)) {
  	numsamples <- length(boot_data_A) + length(boot_data_B)
  }
  maxmatrix <- 2^31 - 1
  if (qboot) {
    # Consistency,
    if (!is.null(boot_data_A) && !is.null(boot_data_B)) {
      # Direct data.
      tx <- boot_data_A[[1]][[1]][order(boot_data_A[[1]][[1]])]
      for (k in 2:length(boot_data_A)){
        if (!all( tx == boot_data_A[[k]][[1]][order(boot_data_A[[k]][[1]])] ))
          return(list("error"=TRUE, "message"="Inconsistent set of transcript IDs across samples!"))
      }
      for (k in 1:length(boot_data_B)){
        if (!all( tx == boot_data_B[[k]][[1]][order(boot_data_B[[k]][[1]])] ))
          return(list("error"=TRUE, "message"="Inconsistent set of transcript IDs across samples!"))
      }
    } else if (!is.null(slo)) {
      # Sleuth.
      tx <- slo$kal[[1]]$bootstrap[[1]][[BS_TARGET_COL]][ order(slo$kal[[1]]$bootstrap[[1]][[BS_TARGET_COL]]) ]
      for (k in 2:length(slo$kal)) {
        if (!all( tx == slo$kal[[k]]$bootstrap[[1]][[BS_TARGET_COL]][ order(slo$kal[[k]]$bootstrap[[1]][[BS_TARGET_COL]]) ] ))
          return(list("error"=TRUE, "message"="Inconsistent set of transcript IDs across samples!"))
      }
    } else {
      warnmsg["noboots"] <- "qboot is TRUE but no bootstrapped estimates were provided! Continuing without bootstrapping."
    }

    # Number of iterations.
    minboots <- infer_bootnum(slo, boot_data_A, boot_data_B)
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
  }
  if (rboot & any( length(samples_by_condition[[1]]) * length(samples_by_condition[[2]]) > maxmatrix/dim(annot)[1],
  				         length(boot_data_A)               * length(boot_data_B)               > maxmatrix/dim(annot)[1] ) )
    warnmsg["toomanyreplicates"] <- "The number of replicates is too high. Exhaustive 1vs1 would exceed maximum capacity of an R matrix."

  return(list("error"=FALSE, "message"="All good!", "maxboots"=minboots, "warn"=(length(warnmsg) > 0), "warnings"= warnmsg))
}


#================================================================================
#' Rule-of-thumb number of iterations.
#'
#' @param slo Sleuth object.
#' @param boot_data_A List of tables of bootstrapped counts.
#' @param boot_data_B List of tables of bootstrapped counts.
#' @return The least number of iterations seen in the input.
#'
infer_bootnum <- function(slo, boot_data_A, boot_data_B){
  minboot <- NA_integer_
  # Bootstrapped counts.
  if (!is.null(boot_data_A) && !is.null(boot_data_B)) {
    minboot <- length(boot_data_A[[1]]) - 1   # ID column
    for (k in 2:length(boot_data_A)){
      minboot <- min(minboot, length(boot_data_A[[k]]) - 1)
    }
    for (k in 1:length(boot_data_B)){
      minboot <- min(minboot, length(boot_data_B[[k]]) - 1)
    }
  # Sleuth.
  } else if (!is.null(slo)) {
    minboot <- length(slo$kal[[1]]$bootstrap)
    for (k in 2:length(slo$kal)) {
      minboot <- min(minboot, length(slo$kal[[k]]$bootstrap))
    }
  }
  return(minboot)
}


#================================================================================
#' Group samples by covariate value.
#' 
#' Sleuth records covariate values per sample. However, RATs needs the reverse: the 
#' samples that correspond to a given covariate value.
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


#================================================================================
#' Extract bootstrap counts into a less nested structure.
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
  # Ensure data.table complies.
  # if (packageVersion("data.table") >= "1.9.8")
    setDTthreads(threads)
  
  # Could just always use tidy_annot(), but that duplicates the annotation without reason. 
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
    dt[, target_id := tx]
    nn <- names(dt)
    ll <- length(nn)
    # Return reordered so that IDs are in first column.
    return(dt[, c(nn[ll], nn[seq.int(1, ll-1)]), with=FALSE])
  }, mc.cores= threads, mc.allow.recursive= FALSE, mc.preschedule= TRUE)
}


#================================================================================
#' Create output structure.
#'
#' @param annot Pre-processed by \code{tidy_annot()}.
#' @param full Full-sized structure or core fields only. Either "full" or "short".
#' @return A list.
#'
#' @import data.table
#'
alloc_out <- function(annot, full){
  # Ensure data.table complies.
  # if (packageVersion("data.table") >= "1.9.8")
    setDTthreads(1)
  if (full == "full") {
    Parameters <- list("description"=NA_character_, "time"=date(),
                       "rats_version"=packageVersion("rats"), "R_version"=R.Version()[c("platform", "version.string")],
                       "var_name"=NA_character_, "cond_A"=NA_character_, "cond_B"=NA_character_,
                       "data_type"=NA_character_, "num_replic_A"=NA_integer_, "num_replic_B"=NA_integer_,
                       "num_genes"=NA_integer_, "num_transc"=NA_integer_,
                       "tests"=NA_character_, "p_thresh"=NA_real_, "abund_thresh"=NA_real_, "dprop_thresh"=NA_real_,
                       "quant_reprod_thresh"=NA_real_, "quant_boot"=NA, "quant_bootnum"=NA_integer_,
                       "rep_reprod_thresh"=NA_real_, "rep_boot"=NA, "rep_bootnum"=NA_integer_, "rep_reprod_as_crit"=NA)
    Genes <- data.table("parent_id"=as.vector(unique(annot$parent_id)),
                        "elig"=NA, "sig"=NA, "elig_fx"=NA, "quant_reprod"=NA, "rep_reprod"=NA, "DTU"=NA, "transc_DTU"=NA,
                        "known_transc"=NA_integer_, "detect_transc"=NA_integer_, "elig_transc"=NA_integer_, "maxDprop"=NA_real_,
                        "pval"=NA_real_, "pval_corr"=NA_real_,
                        "quant_p_mean"=NA_real_, "quant_p_stdev"=NA_real_, "quant_p_min"=NA_real_, "quant_p_max"=NA_real_, 
                        "quant_na_freq"=NA_real_, "quant_dtu_freq"=NA_real_,
                        "rep_p_mean"=NA_real_, "rep_p_stdev"=NA_real_, "rep_p_min"=NA_real_, "rep_p_max"=NA_real_,
                        "rep_na_freq"=NA_real_, "rep_dtu_freq"=NA_real_)
    Transcripts <- data.table("target_id"=annot$target_id, "parent_id"=annot$parent_id,
                              "elig_xp"=NA, "elig"=NA, "sig"=NA, "elig_fx"=NA, "quant_reprod"=NA, "rep_reprod"=NA, "DTU"=NA, "gene_DTU"=NA,
                              "meanA"=NA_real_, "meanB"=NA_real_, "stdevA"=NA_real_, "stdevB"=NA_real_,
                              "sumA"=NA_real_, "sumB"=NA_real_, "log2FC"=NA_real_, "totalA"=NA_real_, "totalB"=NA_real_,
                              "propA"=NA_real_, "propB"=NA_real_, "Dprop"=NA_real_,
                              "pval"=NA_real_, "pval_corr"=NA_real_,
                              "quant_p_mean"=NA_real_, "quant_p_stdev"=NA_real_, "quant_p_min"=NA_real_,"quant_p_max"=NA_real_,
                              "quant_na_freq"=NA_real_, "quant_dtu_freq"=NA_real_,
                              "rep_p_mean"=NA_real_, "rep_p_stdev"=NA_real_, "rep_p_min"=NA_real_,"rep_p_max"=NA_real_,
                              "rep_na_freq"=NA_real_, "rep_dtu_freq"=NA_real_)
    CountData <- list()
  } else {
    Parameters <- list("num_replic_A"=NA_integer_, "num_replic_B"=NA_integer_)
    Genes <- data.table("parent_id"=levels(as.factor(annot$parent_id)), "DTU"=NA,
                        "elig_transc"=NA_integer_, "elig"=NA, "elig_fx"=NA,
                        "pval"=NA_real_, "pval_corr"=NA_real_, "sig"=NA)
    Transcripts <- data.table("target_id"=annot$target_id, "parent_id"=annot$parent_id, "DTU"=NA,
                              "sumA"=NA_real_, "sumB"=NA_real_, "log2FC"=NA_real_,  #log2FC currently not used in decision making and bootstrapping, but maybe in the future.
                              "totalA"=NA_real_, "totalB"=NA_real_,
                              "elig_xp"=NA, "elig"=NA,
                              "propA"=NA_real_, "propB"=NA_real_, "Dprop"=NA_real_, "elig_fx"=NA,
                              "pval"=NA_real_, "pval_corr"=NA_real_, "sig"=NA)
    CountData <- NULL
  }
  with(Genes,
       setkey(Genes, parent_id) )
  with(Transcripts,
       setkey(Transcripts, parent_id, target_id) )

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
#' @param count_thresh Minimum number of counts per replicate.
#' @param p_thresh The p-value threshold.
#' @param dprop_thresh Minimum difference in proportions.
#' @param correction Multiple testing correction type.
#' @param threads Number of threads (POSIX systems only).
#' @return list
#'
#' @import utils
#' @import parallel
#' @import data.table
#'
calculate_DTU <- function(counts_A, counts_B, tx_filter, test_transc, test_genes, full, count_thresh, p_thresh, dprop_thresh, correction, threads= 1) {
  # Ensure data.table complies.
  # if (packageVersion("data.table") >= "1.9.8")
    setDTthreads(threads)

  #---------- PRE-ALLOCATE

  # Pre-allocate results object.
  resobj <- alloc_out(tx_filter, full)
  resobj$Parameters["num_replic_A"] <- dim(counts_A)[2]
  resobj$Parameters["num_replic_B"] <- dim(counts_B)[2]

  with(resobj, {
    # Set key to gene ids.
    setkey(Transcripts, parent_id)
    setkey(Genes, parent_id)

    #---------- STATS

    # Statistics per transcript across all bootstraps per condition, for filtered targets only.
    Transcripts[, sumA :=  rowSums(counts_A, na.rm=TRUE) ]
    Transcripts[, sumB :=  rowSums(counts_B, na.rm=TRUE) ]
    # Sums and proportions, for filtered targets only.
    Transcripts[, totalA := sum(sumA, na.rm=TRUE), by=parent_id]
    Transcripts[, totalB := sum(sumB, na.rm=TRUE), by=parent_id]
    Transcripts[, propA := sumA/totalA]
    Transcripts[, propB := sumB/totalB]
    Transcripts[(is.nan(propA)), propA := NA_real_]  # Replace NaN with NA.
    Transcripts[(is.nan(propB)), propB := NA_real_]
    Transcripts[, Dprop := propB - propA]
    
    #---------- FILTER

    # Filter transcripts and genes to reduce number of tests:
    ctA <- count_thresh * resobj$Parameters[["num_replic_A"]]  # Adjust count threshold for number of replicates.
    ctB <- count_thresh * resobj$Parameters[["num_replic_B"]]
    Transcripts[, elig_xp := (sumA >= ctA | sumB >= ctB)]
    Transcripts[, elig := (elig_xp & totalA != 0 & totalB != 0 & (sumA != totalA | sumB != totalB))]  # If the entire gene is shut off in one condition, changes in proportion cannot be defined.
                                                                                                      # If sum and total are equal in both conditions the gene has only one expressed isoform, or one isoform altogether.
    # If sum and total are equal in both conditions, it has no detected siblings and thus cannot change in proportion.
    Genes[, elig_transc := Transcripts[, as.integer(sum(elig, na.rm=TRUE)), by=parent_id][, V1] ]
    Genes[, elig := elig_transc >= 2]

    # Biologically significant.
    Transcripts[, elig_fx := abs(Dprop) >= dprop_thresh]
    Genes[, elig_fx := Transcripts[, any(elig_fx), by = parent_id][, V1] ]

    #---------- TESTS

    # Transcript-level test.
    if (test_transc) {
      Transcripts[(elig), pval := unlist( mclapply( as.data.frame(t(Transcripts[(elig), .(sumA, sumB, totalA, totalB)])),
                                                    function(row) {
                                                      return( g.test.2(obsx= c(row[1], row[3]-row[1]), obsy= c(row[2], row[4]-row[2])) )
                                                    }, mc.cores= threads, mc.allow.recursive= FALSE, mc.preschedule= TRUE)
                                          ) ]
      Transcripts[(elig), pval_corr := p.adjust(pval, method=correction)]
      Transcripts[(elig), sig := pval_corr < p_thresh]
      Transcripts[(elig), DTU := sig & elig_fx]
    }

    # Gene-level test.
    if (test_genes) {
      Genes[(elig), pval := t( as.data.frame( mclapply(Genes[(elig), parent_id],
                                        function(parent) {
                                            # Extract all relevant data to avoid repeated look ups in the large table.
                                            subdt <- Transcripts[parent, .(sumA, sumB)]  # All isoforms, including not detected ones.
                                            return( g.test.2(obsx= subdt$sumA, obsy= subdt$sumB) )
                                        }, mc.cores= threads, mc.preschedule= TRUE, mc.allow.recursive= FALSE)
                )) ]
      Genes[(elig), pval_corr := p.adjust(pval, method=correction)]
      Genes[(elig), sig := pval_corr < p_thresh]
      Genes[(elig), DTU := sig & elig_fx]
    }
  })
  return(resobj)
}



#================================================================================
#' Log-likelihood test of goodness of fit.
#'
#' @param obsx	a numeric vector of finite positive counts, with at least one non-zero value.
#' @param px	a vector of probabilities of the same length as xobs. ( sum(px) <= 1 )
#'
#' The order of values in the two vectors should be the same.
#' If any value of xp is zero and the corresponding xobs is not, g.test.1 will always reject the hypothesis.
#' No corrections are applied. 
#' No input checks are applied, as RATs needs to run this millions of times.
#' 
#' @import stats
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
#' @param obsx	a numeric vector of positive counts, with at least one non-zero value.
#' @param obsy	a numeric vector of positive counts of same length as obsx, with at least one non-zero value.
#'
#' The order of values in the two vectors should be the same.
#' No corrections are applied. 
#' No input checks are applied, as RATs needs to run this millions of times. 
#'
#' @import stats
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
