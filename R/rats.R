#================================================================================
#' Calculate differential transcript usage.
#'
#' There are three options for input:
#' \itemize{
#'  \item{A sleuth object. This requires the following parameters: \code{slo}, \code{name_A}, \code{name_B} and optionally \code{varname}, \code{COUNT_COL}, \code{BS_TARGET_COL}.}
#'  \item{Bootstrapped count estimates. This requires the following parameters: \code{boot_data_A} and \code{boot_data_B}. \code{name_A} and \code{name_B} can optionally be used to name the conditions.}
#'  \item{Count estimates. This requires the following parameters: \code{count_data_A} and \code{count_data_B}. \code{name_A} and \code{name_B} can optionally be used to name the conditions.}
#' }
#' 
#' @param annot A data.frame matching the transcript identifiers to their corresponding gene identifiers. Any additional columns are allowed but ignored.
#' @param TARGET_COL The name of the transcript identifier column in the \code{annot} object. (Default \code{"target_id"})
#' @param PARENT_COL The name of the parent identifier column in the \code{annot} object. (Default \code{"parent_id"})
#' @param slo A Sleuth object.
#' @param name_A The name for one condition, as it appears in the \code{sample_to_covariates} table within the Sleuth object.
#' @param name_B The name for the other condition, as it appears in the \code{sample_to_covariates} table within the sleuth object.
#' @param varname The name of the covariate to which the two conditions belong, as it appears in the \code{sample_to_covariates} table within the sleuth object. (Default \code{"condition"}).
#' @param COUNTS_COL The name of the counts column to use for the DTU calculation (est_counts or tpm). (Default \code{"est_counts"})
#' @param BS_TARGET_COL The name of the transcript identifier column in the sleuth bootstrap tables. (Default \code{"target_id"})
#' @param count_data_A A data.table of estimated counts for condition A. One column per sample/replicate, one row per transcript. The first column should contain the transcript identifiers.
#' @param count_data_B A data.table of estimated counts for condition B. One column per sample/replicate, one row per transcript. The first column should contain the transcript identifiers.
#' @param boot_data_A A list of data.tables, one per sample/replicate of condition A. One bootstrap iteration's estimates per column, one transcript per row. The first column should contain the transcript identifiers.
#' @param boot_data_B A list of data.tables, one per sample/replicate of condition B. One bootstrap iteration's estimates per column, one transcript per row. The first column should contain the transcript identifiers.
#' @param p_thresh The p-value threshold, default 0.05.
#' @param count_thresh Minimum count of fragments per sample, in at least one of the conditions, for transcripts to be eligible for testing. (Default 10)
#' @param dprop_thresh Minimum change in proportion (effect size) of a transcript for it to be eligible to be significant. (Default 0.1)
#' @param conf_thresh Confidence threshold. The fraction of bootstrap iterations calling DTU required to have confidence in the final call. (Default 0.95) Ignored if no bootstraps.
#' @param correction The p-value correction to apply, as defined in \code{stats::p.adjust.methods}. (Default \code{"BH"})
#' @param testmode One of "genes", "transc", "both". (Default "both")
#' @param boots Bootstrap the p-values of either test. One of "genes", "transc", "both" or "none". (Default "both")
#' @param bootnum Number of bootstraps. (Default 100)
#' @param description Free-text description of the run. You can use this to add metadata to the results object. The results' description field can also be filled in after the run.
#' @param verbose Display progress updates. (Default \code{TRUE})
#' @return List of data tables, with gene-level and transcript-level information.
#'
#' @import utils
#' @import data.table
#' @import matrixStats
#' @export
call_DTU <- function(annot= NULL, TARGET_COL= "target_id", PARENT_COL= "parent_id",
                     slo= NULL, name_A= "Condition-A", name_B= "Condition-B", varname= "condition", COUNTS_COL= "est_counts", BS_TARGET_COL= "target_id",
                     count_data_A = NULL, count_data_B = NULL, boot_data_A = NULL, boot_data_B = NULL,
                     p_thresh= 0.05, count_thresh= 10, dprop_thresh= 0.1, conf_thresh= 0.95, correction= "BH", 
                     testmode= "both", boots= "both", bootnum= 100L, 
                     description=NA_character_, verbose= TRUE)
{
  #---------- PREP
  
  # Input checks.
  if (verbose) {
    message(paste0("Relative Abundance of Transcripts v.", packageVersion("rats")))
    message("Checking parameters...")
  }
  paramcheck <- parameters_good(slo, annot, name_A, name_B, varname, COUNTS_COL,
                                correction, p_thresh, TARGET_COL, PARENT_COL, BS_TARGET_COL, count_thresh, testmode, 
                                boots, bootnum, dprop_thresh, count_data_A, count_data_B, boot_data_A, boot_data_B, conf_thresh)
  if (paramcheck$error) 
    stop(paramcheck$message)
  
  bootnum <- as.integer(bootnum)
  test_transc <- any(testmode == c("transc", "both"))
  test_genes <- any(testmode == c("genes", "both"))
  boot_transc <- any(boots == c("transc", "both"))
  boot_genes <- any(boots == c("genes", "both"))
  
  # Determine data extraction steps.
  steps <- 1  # Assume estimated counts. Simplest case. Bypasses both extraction and averaging.
  if (!is.null(boot_data_A)) {
    steps <- 2  # Bootstrapped estimates. Bypasses extraction, only needs averaging.
  } else if (!is.null(slo)) {
    steps <- 3  # Sleuth object. Most steps. Requires both extraction and averaging.
  }
  
  #----------- LOOK-UP
  
  # Look-up from target_id to parent_id.
  if (verbose)
    message("Creating look-up structures...")
  tx_filter <- data.table(target_id = annot[[TARGET_COL]], parent_id = annot[[PARENT_COL]])
  with( tx_filter,
        setkey(tx_filter, parent_id, target_id) )  # Fixates the order of genes and transcripts to be used throughout the rest of this package.
  # Reverse look-up from replicates to covariates.
  samples_by_condition <- group_samples(slo$sample_to_covariates)[[varname]]
  
  #---------- EXTRACT DATA
  
  if (steps == 3) {
    if (verbose)
      message("Restructuring and aggregating bootstraps...")
    # Re-order rows and collate booted counts in a dataframe per sample. Put dataframes in a list per condition.
    # Target_id is included but NOT used as key so as to ensure that the order keeps matching tx_filter.  
    boot_data_A <- denest_sleuth_boots(slo, tx_filter$target_id, samples_by_condition[[name_A]], COUNTS_COL, BS_TARGET_COL )
    boot_data_B <- denest_sleuth_boots(slo, tx_filter$target_id, samples_by_condition[[name_B]], COUNTS_COL, BS_TARGET_COL )
  } else if (steps == 2) {
    # Just re-order rows.
    boot_data_A <- lapply(boot_data_A, function(x) { x[match(tx_filter$target_id, x[[1]]), ] })
    boot_data_B <- lapply(boot_data_B, function(x) { x[match(tx_filter$target_id, x[[1]]), ] })
  }
  
  if (steps > 1) {  # ID columns are removed so I don't have to constantly subset.
    # Calculate mean count across bootstraps.
    count_data_A <- as.data.table(lapply(boot_data_A, function(b) { n <- names(b); rowMeans(b[, n[2:length(n)], with=FALSE]) }))
    count_data_B <- as.data.table(lapply(boot_data_B, function(b) { n <- names(b); rowMeans(b[, n[2:length(n)], with=FALSE]) }))
    # (data.tables don't have access to column ranges by index, so I have to go roundabout via their names).
    # (Also, the first column is the target_id, so I leave it out).
  } else {
    # Just re-order rows and crop out the transcript IDs.
    nn <- names(count_data_A)
    count_data_A <- count_data_A[match(tx_filter$target_id, count_data_A[[1]]),  nn[seq.int(2, length(nn))],  with=FALSE]
    nn <- names(count_data_B)  # The number of columns may be different from A.
    count_data_B <- count_data_B[match(tx_filter$target_id, count_data_B[[1]]),  nn[seq.int(2, length(nn))],  with=FALSE]
  }
  
  
  
  #---------- TEST
  
  # Do the core work.
  if (verbose)
    message("Calculating significances...")
  suppressWarnings(
    resobj <- calculate_DTU(count_data_A, count_data_B, tx_filter, test_transc, test_genes, "full", count_thresh, p_thresh, dprop_thresh, correction) )
  
  #-------- ADD INFO
  
  if (verbose)
    message("Filling in info and descriptive statistics...")
  
  with(resobj, {
    # Fill in results detais.
    Genes[, known_transc :=  Transcripts[, length(target_id), by=parent_id][, V1] ]  # V1 is the automatic column name for the lengths in the subsetted data.table
    Genes[, detect_transc :=  Transcripts[, .(parent_id, ifelse(sumA + sumB > 0, 1, 0))][, as.integer(sum(V2)), by = parent_id][, V1] ]  # Sum returns double..
    Genes[(is.na(detect_transc)), detect_transc := 0]
    Transcripts[, meanA :=  rowMeans(count_data_A) ]
    Transcripts[, meanB :=  rowMeans(count_data_B) ]
    Transcripts[, stdevA :=  rowSds(as.matrix(count_data_A)) ]
    Transcripts[, stdevB :=  rowSds(as.matrix(count_data_B)) ]
  })
  
  # Fill in run info. (if done within the with() block, changes are local-scoped and don't take effect)
  resobj$Parameters["var_name"] <- varname
  resobj$Parameters["cond_A"] <- name_A
  resobj$Parameters["cond_B"] <- name_B
  resobj$Parameters["p_thresh"] <- p_thresh 
  resobj$Parameters["count_thresh"] <- count_thresh
  resobj$Parameters["dprop_thresh"] <- dprop_thresh
  resobj$Parameters["conf_thresh"] <- conf_thresh
  resobj$Parameters["tests"] <- testmode
  resobj$Parameters["bootstrap"] <- boots
  resobj$Parameters["bootnum"] <- bootnum
  if (steps==3) { 
    resobj$Parameters["data_type"] <- "sleuth" 
  } else if (steps==2) {
    resobj$Parameters["data_type"] <- "bootstrapped abundance estimates"
  } else  if (steps==1) {
    resobj$Parameters["data_type"] <- "plain abundance estimates"
  }
  resobj$Parameters["num_genes"] <- length(levels(annot$parent_id))
  resobj$Parameters["num_transc"] <- length(annot$target_id)
  resobj$Parameters["description"] <- description
  
  #---------- BOOTSTRAP
  
  if (any(boot_transc, boot_genes)) {
    if (verbose)
      message("Bootstrapping...")
    
    # Bootstrapping can take long, so showing progress is nice.
    if (verbose)
      myprogress <- utils::txtProgressBar(min = 0, max = bootnum, initial = 0, char = "=", width = NA, style = 3, file = "")
    
    #----- Iterations
    
    bootres <- lapply(1:bootnum, function(b) {
                  # Update progress.
                  if (verbose)
                    setTxtProgressBar(myprogress, b)
      
                  # Grab a bootstrap from each replicate. 
                  counts_A <- as.data.table(lapply(boot_data_A, function(smpl) { smpl[[sample( names(smpl)[2:(dim(smpl)[2])], 1)]] }))
                  counts_B <- as.data.table(lapply(boot_data_B, function(smpl) { smpl[[sample( names(smpl)[2:(dim(smpl)[2])], 1)]] }))
                  # I have to use list syntax to get a vector back. Usual table syntax with "with=FALSE" returns a table and I fail to cast it.
                  # Also, the first column is the target_id, so I leave it out.
                  
                  # Do the work.
                  # Ignore warning. Chi-square test generates warnings for counts <5. This is expected behaviour. Transcripts changing between off and on are often culprits.
                  suppressWarnings(
                    bout <- calculate_DTU(counts_A, counts_B, tx_filter, test_transc, test_genes, "short", count_thresh, p_thresh, dprop_thresh, correction))
                  
                  with(bout, {
                    return(list("pp" = Transcripts[, pval_corr],
                                "pdtu" = Transcripts[, DTU],
                                "gpab" = Genes[, pvalAB_corr],
                                "gpba" = Genes[, pvalBA_corr],
                                "gdtu" = Genes[, DTU] )) }) })
    
    #----- Stats
    
    if (verbose)
      message("Summarising bootstraps...")
    
    with(resobj, {
      if (boot_transc) {
        # !!! POSSIBLE source of ERRORS if bootstraps * transcripts exceed R's maximum matrix size. (due to number of either) !!!
        pd <- as.matrix(as.data.table(lapply(bootres, function(b) { b[["pdtu"]] })))
        Transcripts[(elig), boot_dtu_freq := rowCounts(pd[Transcripts[, elig], ], value = TRUE, na.rm=TRUE) / bootnum]
        pp <- as.matrix(as.data.table(lapply(bootres, function(b) { b[["pp"]] })))
        Transcripts[(elig), boot_p_mean := rowMeans(pp[Transcripts[, elig], ], na.rm = TRUE)]
        Transcripts[(elig), boot_p_stdev := rowSds(pp[Transcripts[, elig], ], na.rm = TRUE)]
        Transcripts[(elig), boot_p_min := rowMins(pp[Transcripts[, elig], ], na.rm = TRUE)]
        Transcripts[(elig), boot_p_max := rowMaxs(pp[Transcripts[, elig], ], na.rm = TRUE)]
        Transcripts[(elig), boot_na := rowCounts(pp[Transcripts[, elig], ], value = NA, na.rm=FALSE) / bootnum]
        Transcripts[(elig & DTU), conf := (boot_dtu_freq >= conf_thresh)]
        Transcripts[(elig & !DTU), conf := (boot_dtu_freq <= 1-conf_thresh)]
        # Adjust DTU calls.
        Transcripts[(elig), DTU := (DTU & conf)]
      }
      if (boot_genes) {
        # !!! POSSIBLE source of ERRORS if bootstraps * genes exceed R's maximum matrix size. (due to number of bootstraps) !!!
        gabres <- as.matrix(as.data.table(lapply(bootres, function(b) { b[["gpab"]] })))
        gbares <- as.matrix(as.data.table(lapply(bootres, function(b) { b[["gpba"]] })))
        gdres <- as.matrix(as.data.table(lapply(bootres, function(b) { b[["gdtu"]] })))
        Genes[(elig), boot_dtu_freq := rowCounts(gdres[Genes[, elig], ], value = TRUE, na.rm = TRUE) / bootnum]
        Genes[(elig), boot_p_meanAB := rowMeans(gabres[Genes[, elig], ], na.rm = TRUE)]
        Genes[(elig), boot_p_meanBA := rowMeans(gbares[Genes[, elig], ], na.rm = TRUE)]
        Genes[(elig), boot_p_stdevAB := rowSds(gabres[Genes[, elig], ], na.rm = TRUE)]
        Genes[(elig), boot_p_stdevBA := rowSds(gbares[Genes[, elig], ], na.rm = TRUE)]
        Genes[(elig), boot_p_minAB := rowMins(gabres[Genes[, elig], ], na.rm = TRUE)]
        Genes[(elig), boot_p_minBA := rowMins(gbares[Genes[, elig], ], na.rm = TRUE)]
        Genes[(elig), boot_p_maxAB := rowMaxs(gabres[Genes[, elig], ], na.rm = TRUE)]
        Genes[(elig), boot_p_maxBA := rowMaxs(gbares[Genes[, elig], ], na.rm = TRUE)]
        Genes[(elig), boot_na := rowCounts(gabres[Genes[, elig], ], value = NA, na.rm = FALSE) / bootnum]  # It doesn't matter if AB or BA, affected identically by gene eligibility.
        Genes[(elig & DTU), conf := (boot_dtu_freq >= conf_thresh)]
        Genes[(elig & !DTU), conf := (boot_dtu_freq <= 1-conf_thresh)]
        # Adjust DTU calls.
        Genes[(elig), DTU := (DTU & conf)]
      }
    })
  }
  
  #---------- DONE
  
  if (verbose)
    message("Tidying up...")
  
  # Store the replicate means adter re-adding the IDs.
  with(count_data_A, {
    count_data_A[,  target_id := tx_filter$target_id]
    count_data_A[,  parent_id := tx_filter$parent_id]
    setkey(count_data_A, parent_id)
  })
  with(count_data_B, {
    count_data_B[,  target_id := tx_filter$target_id]
    count_data_B[,  parent_id := tx_filter$parent_id]
    setkey(count_data_B, parent_id)
  })
  resobj$ReplicateData <- list("condA"= count_data_A, "condB"= count_data_B)
  
  with(resobj, {
    # Cross-display the DTU calls.
    if (test_transc)
      Genes[, transc_DTU := Transcripts[, any(DTU), by = parent_id][, V1] ]
    if (test_genes)
      Transcripts[, gene_DTU := merge(Genes[, .(parent_id, DTU)], Transcripts[, .(parent_id)])[, DTU] ]
    
    # Drop the bootstrap columns, if unused.
    if (!boot_transc)
      Transcripts[, c("boot_dtu_freq", "boot_p_mean", "boot_p_stdev", "boot_p_min", "boot_p_max", "boot_na") := NULL]
    if(!boot_genes)
      Genes[, c("boot_dtu_freq", "boot_p_meanAB", "boot_p_meanBA", "boot_p_stdevAB", "boot_p_stdevBA", "boot_p_minAB",
              "boot_p_minBA", "boot_p_maxAB", "boot_p_maxBA", "boot_na") := NULL]
  })
  
  if(verbose) {
    message("All done!")
    print(dtu_summary(resobj))
  }
  
  return(resobj)
}









#================================================================================
#================================================================================





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
#' 
#' @return List with a logical value and a message.
#'
parameters_good <- function(slo, annot, name_A, name_B, varname, COUNTS_COL,
                            correction, p_thresh, TARGET_COL, PARENT_COL, BS_TARGET_COL, 
                            count_thresh, testmode, boots, bootnum, dprop_thresh,
                            count_data_A, count_data_B, boot_data_A, boot_data_B, conf_thresh) {
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
  if ((!is.numeric(bootnum)) || bootnum < 1)
    return(list("error"=TRUE, "message"="Invalid number of bootstraps! Must be a positive number."))
  if (!boots %in% c("genes", "transc", "both", "none"))
    return(list("error"=TRUE, "message"="Unrecognized value for boots!"))
  if ((!is.numeric(conf_thresh)) || conf_thresh < 0 || conf_thresh > 1)
    return(list("error"=TRUE, "message"="Invalid confidence threshold! Must be between 0 and 1."))
  
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
  if (!is.null(boot_data_A)  && is.null(slo))  # If slo available ignore boot_data.
    if (any(!is.list(boot_data_A), !is.list(boot_data_A), !is.data.table(boot_data_A[[1]]), !is.data.table(boot_data_B[[1]]) ))
        return(list("error"=TRUE, "message"="The bootstrap data are not lists of data.tables!"))
  
  # Counts.
  if (!is.null(count_data_A)  && is.null(slo) && any(is.null(boot_data_A), is.null(boot_data_B))){  # If slo or boot_data available, ignore count_data.
    if (any(!is.data.table(count_data_A), !is.data.table(count_data_B)))
      return(list("error"=TRUE, "message"="The counts data are not data.tables!"))
    if (!all( count_data_A[[1]][order(count_data_A[[1]])] == count_data_B[[1]][order(count_data_B[[1]])] ))
      return(list("error"=TRUE, "message"="Inconsistent set of transcript IDs between conditions!"))
  }
  
  # Bootstrap.
  if (boots != "none") {
    # Direct data.
    if (!is.null(boot_data_A) && !is.null(boot_data_B)) {
      tx <- boot_data_A[[1]][[1]][order(boot_data_A[[1]][[1]])]
      for (k in 2:length(boot_data_A)){
        if (!all( tx == boot_data_A[[k]][[1]][order(boot_data_A[[k]][[1]])] ))
          return(list("error"=TRUE, "message"="Inconsistent set of transcript IDs across samples!"))
      }
      for (k in 2:length(boot_data_B)){
        if (!all( tx == boot_data_B[[k]][[1]][order(boot_data_B[[k]][[1]])] ))
          return(list("error"=TRUE, "message"="Inconsistent set of transcript IDs across samples!"))
      }
    # Sleuth.  
    } else if (!is.null(slo)) {
      tx <- slo$kal[[1]]$bootstrap[[1]][[BS_TARGET_COL]][ order(slo$kal[[1]]$bootstrap[[1]][[BS_TARGET_COL]]) ]
      for (k in 2:length(slo$kal)) {
        if (!all( tx == slo$kal[[k]]$bootstrap[[1]][[BS_TARGET_COL]][ order(slo$kal[[k]]$bootstrap[[1]][[BS_TARGET_COL]]) ] ))
          return(list("error"=TRUE, "message"="Inconsistent set of transcript IDs across samples!"))
      }
    } else {
      return(list("error"=TRUE, "message"="No bootstrapped estimates were provided!"))
    }
  } 
  
  return(list("error"=FALSE, "message"="All good!"))
}


#================================================================================
#' Group sample numbers by factor.
#'
#' @param covariates a dataframe with different factor variables.
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
#' NA replaced with 0. Means counts per transcripts are calculated and included as
#' a column per sample.
#' 
#' @param slo A sleuth object.
#' @param tx A vector of transcript ids. The results will be ordered according to this vector.
#' @param samples A numeric vector of samples to extract counts for.
#' @param COUNTS_COL The name of the column with the counts.
#' @param BS_TARGET_COL The name of the column with the transcript IDs.
#' @return A list of data.tables, one per sample, containing all the bootstrap counts of the smaple. First column contains the transcript IDs.
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
                       "rats_version"=packageVersion("rats"), "R_version"=R.Version())
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
    Transcripts[, sumA :=  rowSums(counts_A) ]
    Transcripts[, sumB :=  rowSums(counts_B) ]
    # Sums and proportions, for filtered targets only.
    Transcripts[, totalA := sum(sumA), by=parent_id]
    Transcripts[, totalB := sum(sumB), by=parent_id]
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
    Genes[, elig_transc := Transcripts[, as.integer(sum(elig)), by=parent_id][, V1] ]
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
                       subdt <- Transcripts[parent, .(sumA, sumB, propA, propB)]
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

