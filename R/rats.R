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
#' @param testmode One of \itemize{\item{"genes"}, \item{"transc"}, \item{"both" (default)}}.
#' @param boots Bootstrap the p-values of either test. One of \itemize{\item{"genes"}, \item{"transc"}, \item{"both" (default)} \item{"none"}}.
#' @param bootnum Number of bootstraps. (if 0, bootnum will be infered from the data)
#' @param description Free-text description of the run. You can use this to add metadata to the results object. The results' description field can also be filled in after the run.
#' @param verbose Display progress updates and warnings. (Default \code{TRUE})
#' @param threads Number of threads to use (only on POSIX architectures). (Default 1)
#' @param dbg Prematurely terminate execution at the specified stage. Used to speed up tests by avoiding unnecessary downstream processing. (Default 0: do not interrupt)
#' @return List of data tables, with gene-level and transcript-level information.
#'
#' @import utils
#' @import parallel
#' @import data.table
#' @import matrixStats
#' @export
call_DTU <- function(annot= NULL, TARGET_COL= "target_id", PARENT_COL= "parent_id",
                     slo= NULL, name_A= "Condition-A", name_B= "Condition-B", varname= "condition", COUNTS_COL= "est_counts", BS_TARGET_COL= "target_id",
                     count_data_A = NULL, count_data_B = NULL, boot_data_A = NULL, boot_data_B = NULL,
                     p_thresh= 0.05, count_thresh= 10, dprop_thresh= 0.1, conf_thresh= 0.95, correction= "BH", 
                     testmode= "both", boots= "both", bootnum= 0L, 
                     description= NA_character_, verbose= TRUE, threads= 1L, dbg= 0)
{
  #---------- PREP
  
  
  if (verbose) {
    message(paste0("Relative Abundance of Transcripts v.", packageVersion("rats")))
    message("Checking parameters...")
  }
  # Input checks.
  paramcheck <- parameters_are_good(slo, annot, name_A, name_B, varname, COUNTS_COL,
                                correction, p_thresh, TARGET_COL, PARENT_COL, BS_TARGET_COL, count_thresh, testmode, 
                                boots, bootnum, dprop_thresh, count_data_A, count_data_B, boot_data_A, boot_data_B, conf_thresh, threads)
  if (paramcheck$error) 
    stop(paramcheck$message)
  if (verbose)
    if(paramcheck$warn)
      for (w in paramcheck$warnings) {
        warning(w)  # So it displays as warning at the end of the run.
        message(w)  # So it displays at runtime.
      }
  
  if (bootnum == 0 && boots != "none")   # Use smart default.
    bootnum = paramcheck$maxboots
  bootnum <- as.integer(bootnum)  # Can't be decimal.
  threads <- as.integer(threads)  # Can't br decimal.
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
  
  if (dbg == 1)
    return(NULL)
  
  #----------- LOOK-UP
  
  # Look-up from target_id to parent_id.
  if (verbose)
    message("Creating look-up structures...")
  tx_filter <- data.table(target_id = annot[[TARGET_COL]], parent_id = annot[[PARENT_COL]])
  with( tx_filter,
        setkey(tx_filter, parent_id, target_id) )  # Fixates the order of genes and transcripts to be used throughout the rest of this package.
  # Reverse look-up from replicates to covariates.
  samples_by_condition <- group_samples(slo$sample_to_covariates)[[varname]]
  
  if (dbg == 2)
    return(NULL)
  
  #---------- EXTRACT DATA
  
  if (steps == 3) {   # From Sleuth
    if (verbose)
      message("Restructuring and aggregating bootstraps...")
    # Re-order rows and collate booted counts in a dataframe per sample. Put dataframes in a list per condition.
    # Target_id is included but NOT used as key so as to ensure that the order keeps matching tx_filter.  
    boot_data_A <- denest_sleuth_boots(slo, tx_filter$target_id, samples_by_condition[[name_A]], COUNTS_COL, BS_TARGET_COL )
    boot_data_B <- denest_sleuth_boots(slo, tx_filter$target_id, samples_by_condition[[name_B]], COUNTS_COL, BS_TARGET_COL )
  } else if (steps == 2) {    # From generic bootstrapped data  
    # Just re-order rows.
    boot_data_A <- lapply(boot_data_A, function(x) { x[match(tx_filter$target_id, x[[1]]), ] })
    boot_data_B <- lapply(boot_data_B, function(x) { x[match(tx_filter$target_id, x[[1]]), ] })
  }
  
  if (dbg == 3)
    return(list("bootA"=count_data_A, "bootB"=count_data_B))
  
  # ID columns are removed so I don't have to constantly subset.
  if (steps > 1) {    # Either Sleuth or generic bootstraps.
    # Calculate mean count across bootstraps.
    count_data_A <- as.data.table(lapply(boot_data_A, function(b) { n <- names(b); rowMeans(b[, n[2:length(n)], with=FALSE]) }))
    count_data_B <- as.data.table(lapply(boot_data_B, function(b) { n <- names(b); rowMeans(b[, n[2:length(n)], with=FALSE]) }))
    # (data.tables don't have access to column ranges by index, so I have to go roundabout via their names).
    # (Also, the first column is the target_id, so I leave it out).
  } else {    # From generic unbootstrapped data.
    # Just re-order rows and crop out the transcript IDs.
    nn <- names(count_data_A)
    count_data_A <- count_data_A[match(tx_filter$target_id, count_data_A[[1]]),  nn[seq.int(2, length(nn))],  with=FALSE]
    nn <- names(count_data_B)  # The number of columns may be different from A.
    count_data_B <- count_data_B[match(tx_filter$target_id, count_data_B[[1]]),  nn[seq.int(2, length(nn))],  with=FALSE]
  }
  
  if (dbg == 4)
    return(list("countA"=count_data_A, "countB"=count_data_B))
  
  #---------- TEST
  
  # Do the core work.
  if (verbose)
    message("Calculating significances...")
  suppressWarnings(
    resobj <- calculate_DTU(count_data_A, count_data_B, tx_filter, test_transc, test_genes, "full", count_thresh, p_thresh, dprop_thresh, correction, threads) )
  
  if (dbg == 5)
    return(resobj)
  
  #-------- ADD INFO
  
  if (verbose)
    message("Filling in info and descriptive statistics...")
  
  with(resobj, {
    # Fill in results detais.
    Genes[, known_transc :=  Transcripts[, length(target_id), by=parent_id][, V1] ]  # V1 is the automatic column name for the lengths in the subsetted data.table
    Genes[, detect_transc :=  Transcripts[, .(parent_id, ifelse(sumA + sumB > 0, 1, 0))][, as.integer(sum(V2)), by = parent_id][, V1] ]  # Sum returns type double.
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
  
  if (dbg == 6)
    return(resobj)
  
  #---------- BOOTSTRAP
  
  if (any(boot_transc, boot_genes)) {
    if (verbose)
      message("Bootstrapping...")
    
    # Bootstrapping can take long, so showing progress is nice.
    # if (verbose)
    #   myprogress <- utils::txtProgressBar(min = 0, max = bootnum, initial = 0, char = "=", width = NA, style = 3, file = "")
    
    #----- Iterations
    
    bootres <- mclapply(1:bootnum, function(b) {
                  # Update progress.
                  # if (verbose)
                    # setTxtProgressBar(myprogress, b)
      
                  # Grab a bootstrap from each replicate. 
                  counts_A <- as.data.table(lapply(boot_data_A, function(smpl) { smpl[[sample( names(smpl)[2:(dim(smpl)[2])], 1)]] }))
                  counts_B <- as.data.table(lapply(boot_data_B, function(smpl) { smpl[[sample( names(smpl)[2:(dim(smpl)[2])], 1)]] }))
                  # I have to use list syntax to get a vector back. Usual table syntax with "with=FALSE" returns a table and I fail to cast it.
                  # Also, the first column is the target_id, so I leave it out.
                  
                  # Do the work.
                  # Ignore warning. Chi-square test generates warnings for counts <5. This is expected behaviour. Transcripts changing between off and on are often culprits.
                  suppressWarnings(
                    bout <- calculate_DTU(counts_A, counts_B, tx_filter, test_transc, test_genes, "short", count_thresh, p_thresh, dprop_thresh, correction, threads=1) )
                  
                  with(bout, {
                    return(list("pp" = Transcripts[, pval_corr],
                                "pdtu" = Transcripts[, DTU],
                                "gpab" = Genes[, pvalAB_corr],
                                "gpba" = Genes[, pvalBA_corr],
                                "gdtu" = Genes[, DTU] )) }) 
                },
                mc.cores= threads, mc.allow.recursive= FALSE, mc.preschedule= TRUE)
    if (verbose)  # Forcing a new line after the progress bar.
      message("")
    
    if (dbg == 7)
      return(resobj)
    
    #----- Stats
    
    if (verbose)
      message("Summarising bootstraps...")
    
    with(resobj, {
      if (boot_transc) {
        # !!! POSSIBLE source of ERRORS if bootstraps * transcripts exceed R's maximum matrix size. (due to number of either) !!!
        pd <- as.matrix(as.data.table(mclapply(bootres, function(b) { b[["pdtu"]] }, mc.cores= threads)))
        Transcripts[(elig), boot_dtu_freq := rowCounts(pd[Transcripts[, elig], ], value = TRUE, na.rm=TRUE) / bootnum]
        pp <- as.matrix(as.data.table(mclapply(bootres, function(b) { b[["pp"]] }, mc.cores= threads)))
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
        gabres <- as.matrix(as.data.table(mclapply(bootres, function(b) { b[["gpab"]] }, mc.cores= threads)))
        gbares <- as.matrix(as.data.table(mclapply(bootres, function(b) { b[["gpba"]] }, mc.cores= threads)))
        gdres <- as.matrix(as.data.table(mclapply(bootres, function(b) { b[["gdtu"]] }, mc.cores= threads)))
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
  
  if (dbg == 8)
    return(resobj)
  
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
      Genes[, transc_DTU := Transcripts[, any(DTU, na.rm=TRUE), by = parent_id][, V1] ]
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


