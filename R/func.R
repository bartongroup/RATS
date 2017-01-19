#================================================================================
#' Check input parameters.
#' 
#' @param slo sleuth object 
#' @param annot annotation dataframe
#' @param name_A condition name
#' @param name_B condition name
#' @param varname name of condition variable
#' @param COUNTS_COL name of counts column in bootstrap
#' @param correction p-value correction method
#' @param p_thresh significance level
#' @param TARGET_COL name of transcript id column in annotation
#' @param PARENT_COL name of gene id column in annotation
#' @param BS_TARGET_COL name of transcript id column in bootstrap
#' @param count_thresh minimum frgments per transcript per sample
#' @param testmode which tests to run
#' @param boots which tests to bootstrap
#' @param bootnum number of bootstrap iterations
#' @param dprop_thresh minimum change in proportion
#' @param count_data_A A dataframe of estimated counts.
#' @param count_data_B A dataframe of estimated counts.
#' @param boot_data_A A list of dataframes, one per sample, each with all the bootstrapped estimetes for the sample.
#' @param boot_data_B A list of dataframes, one per sample, each with all the bootstrapped estimetes for the sample.
#' @param conf_thresh Confidence threshold.
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
parameters_are_good <- function(slo, annot, name_A, name_B, varname, COUNTS_COL,
                            correction, p_thresh, TARGET_COL, PARENT_COL, BS_TARGET_COL, 
                            count_thresh, testmode, boots, bootnum, dprop_thresh,
                            count_data_A, count_data_B, boot_data_A, boot_data_B, conf_thresh, threads) {
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
    return(list("error"=TRUE, "message"="Invalid read-count threshold! Must be between 0 and 1."))
  if ((!is.numeric(bootnum)) || bootnum < 0)
    return(list("error"=TRUE, "message"="Invalid number of bootstraps! Must be a positive integer number."))
  if (!boots %in% c("genes", "transc", "both", "none"))
    return(list("error"=TRUE, "message"="Unrecognized value for boots!"))
  if ((!is.numeric(conf_thresh)) || conf_thresh < 0 || conf_thresh > 1)
    return(list("error"=TRUE, "message"="Invalid confidence threshold! Must be between 0 and 1."))
  if ((!is.numeric(threads)) || threads <= 0)
    return(list("error"=TRUE, "message"="Invalid number of threads! Must be positive integer."))
  if (threads > parallel::detectCores(logical= TRUE))
    return(list("error"=TRUE, "message"="Number of threads exceeds system's reported capacity."))
  
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
  if (boots != "none") {
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
      return(list("error"=TRUE, "message"="No bootstrapped estimates were provided!"))
    }
    
    # Number of iterations.
    minboots <- infer_bootnum(slo, boot_data_A, boot_data_B)
    if (is.na(minboots) || minboots == 0)
      return(list("error"=TRUE, "message"="It appears some of your samples have no bootstraps!"))
    if (bootnum > minboots)
      warnmsg["toomanyboots"] <- "You are requesting more RATs bootstrap iterations than available in your quantification data! Are you sure about this?"
    if (100 > minboots)
      warnmsg["toofewboots"] <- "You are requesting fairly few bootstrap iterations, which reduces reproducibility of the calls. Are you sure about this?"
  } 
  
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
#' Group sample numbers by factor.
#'
#' @param covariates A dataframe with different factor variables.
#' @return list of lists (per covariate) of vectors (per factor level).
#'
#' Row number corresponds to sample number.
#'
group_samples <- function(covariates) {
  samplesByVariable <- list()
  for(varname in names(covariates)) {
    categories <- levels(as.factor(covariates[, varname]))
    samplesByVariable[[varname]] <- list()
    for (x in categories) {
      samplesByVariable[[varname]][[x]] <- which(covariates[, varname] == x)
    }
  }
  return(samplesByVariable)
}


#================================================================================
#' Extract bootstrap counts into a less nested structure.
#' 
#' @param slo A sleuth object.
#' @param tx A vector of transcript ids. The results will be ordered according to this vector.
#' @param samples A numeric vector of samples to extract counts for.
#' @param COUNTS_COL The name of the column with the counts.
#' @param BS_TARGET_COL The name of the column with the transcript IDs.
#' @return A list of data.tables, one per sample, containing all the bootstrap counts of the smaple. First column contains the transcript IDs.
#'
#' NA replaced with 0. 
#' 
#' Transcripts in \code{slo} that are missing from \code{tx} will be skipped completely.
#' Transcripts in \code{tx} that are missing from \code{slo} are automatically padded with NA, which we re-assign as 0.
#'
denest_sleuth_boots <- function(slo, tx, samples, COUNTS_COL, BS_TARGET_COL) {
  lapply(samples, function(smpl) {
    # Extract counts in the order of provided transcript vector, for safety and consistency.
    dt <- as.data.table( lapply(slo$kal[[smpl]]$bootstrap, function(boot) {
      roworder <- match(tx, boot[[BS_TARGET_COL]])
      boot[roworder, COUNTS_COL]
    }))
    # Replace any NAs with 0. Happens when annotation different from that used for DTE.
    dt[is.na(dt)] <- 0
    # Add transcript ID.
    dt[, "target_id" := tx]
    nn <- names(dt)
    ll <- length(nn)
    # Return reordered so that IDs are in first column.
    return(dt[, c(nn[ll], nn[seq.int(1, ll-1)]), with=FALSE])
  })
}


#================================================================================
#' Create output structure.
#' 
#' @param annot A dataframe with at least 2 columns: target_id and parent_id.
#' @param full Full-sized structure or core fields only. Either "full" or "short".
#' @return A list.
#' 
alloc_out <- function(annot, full){
  if (full == "full") {
    Parameters <- list("var_name"=NA_character_, "cond_A"=NA_character_, "cond_B"=NA_character_,
                       "num_replic_A"=NA_integer_, "num_replic_B"=NA_integer_,
                       "p_thresh"=NA_real_, "count_thresh"=NA_real_, "dprop_thresh"=NA_real_, "conf_thresh"=NA_real_,
                       "tests"=NA_character_, "bootstrap"=NA_character_, "bootnum"=NA_integer_,
                       "data_type"=NA_character_, "num_genes"=NA_integer_, "num_transc"=NA_integer_,
                       "description"=NA_character_,
                       "rats_version"=packageVersion("rats"), "R_version"=R.Version()[c("platform", "version.string")])
    Genes <- data.table("parent_id"=levels(as.factor(annot$parent_id)),
                        "DTU"=NA, "transc_DTU"=NA, 
                        "known_transc"=NA_integer_, "detect_transc"=NA_integer_, "elig_transc"=NA_integer_,
                        "elig"=NA, "elig_fx"=NA,
                        "pvalAB"=NA_real_, "pvalBA"=NA_real_, "pvalAB_corr"=NA_real_, "pvalBA_corr"=NA_real_, "sig"=NA, 
                        "boot_dtu_freq"=NA_real_, "conf"=NA, "boot_p_meanAB"=NA_real_, "boot_p_meanBA"=NA_real_, 
                        "boot_p_stdevAB"=NA_real_, "boot_p_stdevBA"=NA_real_, "boot_p_minAB"=NA_real_, "boot_p_minBA"=NA_real_, 
                        "boot_p_maxAB"=NA_real_, "boot_p_maxBA"=NA_real_, "boot_na"=NA_real_)
    Transcripts <- data.table("target_id"=annot$target_id, "parent_id"=annot$parent_id,
                              "DTU"=NA, "gene_DTU"=NA,
                              "meanA"=NA_real_, "meanB"=NA_real_,  # mean across replicates of means across bootstraps
                              "stdevA"=NA_real_, "stdevB"=NA_real_,  # standard deviation across replicates of means across bootstraps
                              "sumA"=NA_real_, "sumB"=NA_real_,  # sum across replicates of means across bootstraps
                              "totalA"=NA_real_, "totalB"=NA_real_,  # sum of all transcripts for that gene
                              "elig_xp"=NA, "elig"=NA,
                              "propA"=NA_real_, "propB"=NA_real_, "Dprop"=NA_real_, "elig_fx"=NA,
                              "pval"=NA_real_,  "pval_corr"=NA_real_, "sig"=NA, 
                              "boot_dtu_freq"=NA_real_, "conf"=NA, "boot_p_mean"=NA_real_, "boot_p_stdev"=NA_real_, 
                              "boot_p_min"=NA_real_,"boot_p_max"=NA_real_, "boot_na"=NA_real_)
    ReplicateData <- list()
  } else {
    Parameters <- list("num_replic_A"=NA_integer_, "num_replic_B"=NA_integer_)
    Genes <- data.table("parent_id"=levels(as.factor(annot$parent_id)), "DTU"=NA, 
                        "elig_transc"=NA_integer_, "elig"=NA, "elig_fx"=NA,
                        "pvalAB"=NA_real_, "pvalBA"=NA_real_,
                        "pvalAB_corr"=NA_real_, "pvalBA_corr"=NA_real_, "sig"=NA)
    Transcripts <- data.table("target_id"=annot$target_id, "parent_id"=annot$parent_id, "DTU"=NA, 
                              "sumA"=NA_real_, "sumB"=NA_real_,  # sum across replicates of means across bootstraps
                              "totalA"=NA_real_, "totalB"=NA_real_,  # sum of all transcripts for that gene
                              "elig_xp"=NA, "elig"=NA,
                              "propA"=NA_real_, "propB"=NA_real_, "Dprop"=NA_real_, "elig_fx"=NA,
                              "pval"=NA_real_, "pval_corr"=NA_real_, "sig"=NA)
    ReplicateData <- NULL
  }
  with(Genes,
       setkey(Genes, parent_id) )
  with(Transcripts, 
       setkey(Transcripts, parent_id, target_id) )
  
  return(list("Parameters"= Parameters, "Genes"= Genes, "Transcripts"= Transcripts, "ReplicateData"= ReplicateData))
}


#================================================================================
#' Set up and execute the tests.
#'
#' @param counts_A A data.table of counts for condition A. x: sample, y: transcript.
#' @param counts_B A data.table of counts for condition B. x: sample, y: transcript.
#' @param tx_filter A data.table with target_id and parent_id.
#' @param test_transc Whether to do transcript-level test.
#' @param test_genes Whether to do gene-level test.
#' @param full Either "full" (for complete output structure) or "short" (for bootstrapping).
#' @param count_thresh Minimum number of counts per replicate.
#' @param p_thresh The p-value threshold.
#' @param dprop_thresh Minimum difference in proportions.
#' @param correction Multiple testing correction type.
#' @return list
#' 
#' @import data.table
#' 
calculate_DTU <- function(counts_A, counts_B, tx_filter, test_transc, test_genes, full, count_thresh, p_thresh, dprop_thresh, correction) {
  
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
    Transcripts[, elig := (elig_xp & totalA != 0 & totalB != 0 & (sumA != totalA | sumB != totalB))]  # If the entire gene is shut off, changes in proportion cannot be defined.
    # If sum and total are equal in both conditions, it has no detected siblings and thus cannot change in proportion.
    Genes[, elig_transc := Transcripts[, as.integer(sum(elig, na.rm=TRUE)), by=parent_id][, V1] ]
    Genes[, elig := elig_transc >= 2]
    
    # Biologically significant.
    Transcripts[, elig_fx := abs(Dprop) >= dprop_thresh]
    Genes[, elig_fx := Transcripts[, any(elig_fx), by = parent_id][, V1] ]
    
    #---------- TESTS
    
    # Proportion test.
    if (test_transc) {
      Transcripts[(elig), pval := as.vector(apply(Transcripts[(elig), .(sumA, sumB, totalA, totalB)], MARGIN = 1, 
                                                  FUN = function(row) { prop.test(x = row[c("sumA", "sumB")], 
                                                                                  n = row[c("totalA", "totalB")], 
                                                                                  correct = TRUE)[["p.value"]] } )) ]
      Transcripts[(elig), pval_corr := p.adjust(pval, method=correction)]
      Transcripts[(elig), sig := pval_corr < p_thresh]
      Transcripts[(elig), DTU := sig & elig_fx]
    }
    
    # G test.
    if (test_genes) {
      Genes[(elig), c("pvalAB", "pvalBA") := 
              as.data.frame( t( as.data.frame( lapply(Genes[(elig), parent_id], function(parent) {
                # Extract all relevant data to avoid repeated look ups in the large table.
                subdt <- Transcripts[parent, .(sumA, sumB, propA, propB)]  # All isoforms, including not detected ones.
                pAB <- g.test(x = subdt[, sumA], p = subdt[, propB])
                pBA <- g.test(x = subdt[, sumB], p = subdt[, propA])
                return(c(pAB, pBA)) }) ))) ]
      Genes[(elig), pvalAB_corr := p.adjust(pvalAB, method=correction)]
      Genes[(elig), pvalBA_corr := p.adjust(pvalBA, method=correction)]
      Genes[(elig), sig := pvalAB_corr < p_thresh & pvalBA_corr < p_thresh]
      Genes[(elig), DTU := sig & elig_fx]
    }
  })
  return(resobj)
}


#================================================================================
#' Log-likelihood test of goodness of fit.
#'
#' @param x	a numeric vector of positive numbers, with at least one non-zero value.
#' @param p	a vector of probabilities of the same length of x.
#'
#' Sourced and adapted from from:
#' V3.3 Pete Hurd Sept 29 2001. phurd@ualberta.ca
#' http://www.psych.ualberta.ca/~phurd/cruft/g.test.r
#'
g.test <- function(x, p = rep(1/length(x), length(x)))
{
  n <- sum(x)
  E <- n * p
  names(E) <- names(x)
  g <- 0
  for (i in 1:length(x)){
    if (x[i] != 0) g <- g + x[i] * log(x[i]/E[i])
  }
  q <- 1
  STATISTIC <- G <- 2*g/q
  PARAMETER <- length(x) - 1
  PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
}

