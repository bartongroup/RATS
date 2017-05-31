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
#' @param annot A data.table matching transcript identifiers to gene identifiers. Any additional columns are allowed but ignored.
#' @param TARGET_COL The name of the column for the transcript identifiers in \code{annot}. (Default \code{"target_id"})
#' @param PARENT_COL The name of the column for the gene identifiers in \code{annot}. (Default \code{"parent_id"})
#' @param slo A Sleuth object.
#' @param name_A The name for one condition, as it appears in the \code{sample_to_covariates} table within the Sleuth object.  (Default "Condition-A")
#' @param name_B The name for the other condition, as it appears in the \code{sample_to_covariates} table within the sleuth object.  (Default "Condition-B")
#' @param varname The name of the covariate to which the two conditions belong, as it appears in the \code{sample_to_covariates} table within the sleuth object. (Default \code{"condition"}).
#' @param COUNTS_COL For Sleuth objects only. The name of the counts column to use for the DTU calculation (est_counts or tpm). (Default \code{"est_counts"})
#' @param BS_TARGET_COL For Sleuth objects only. The name of the transcript identifiers column in the bootstrap tables. (Default \code{"target_id"})
#' @param count_data_A A data.table of estimated counts for condition A. One column per sample/replicate, one row per transcript. The first column should contain the transcript identifiers.
#' @param count_data_B A data.table of estimated counts for condition B. One column per sample/replicate, one row per transcript. The first column should contain the transcript identifiers.
#' @param boot_data_A A list of data.tables, one per sample/replicate of condition A. One bootstrap iteration's estimates per column, one transcript per row. The first column should contain the transcript identifiers.
#' @param boot_data_B A list of data.tables, one per sample/replicate of condition B. One bootstrap iteration's estimates per column, one transcript per row. The first column should contain the transcript identifiers.
#' @param p_thresh The p-value threshold. (Default 0.05)
#' @param abund_thresh Noise threshold. Minimum mean abundance for transcripts to be eligible for testing. (Default 5)
#' @param dprop_thresh Effect size threshold. Minimum change in proportion of a transcript for it to be considered meaningful. (Default 0.20)
#' @param correction The p-value correction to apply, as defined in \code{stats::p.adjust.methods}. (Default \code{"BH"})
#' @param scaling A scaling factor to be applied to the abundances, *prior* to any thresholding and testing. Useful for scaling TPM abundances to the actual number of transcripts in the sequencing sample, if you use TPMs and if you can infer that information from your sample preparation. Improper use of the scaling factor will artificially inflate/deflate the significances obtained. (Default 1)
#' @param testmode One of \itemize{\item{"genes"}, \item{"transc"}, \item{"both" (default)}}.
#' @param qboot Bootstrap the DTU robustness against bootstrapped quantifications data. (Default \code{TRUE}) Ignored if input is \code{count_data}.
#' @param qbootnum Number of iterations for \code{qboot}. (Default 0) If 0, RATs will try to infer a value from the data.
#' @param qrep_thresh Reproducibility threshold for quantification bootsrapping. (Default 0.95)
#' @param rboot Bootstrap the DTU robustness against the replicates. Does ALL 1 vs 1 combinations. (Default \code{FALSE})
#' @param rrep_thresh Reproducibility threshold for replicate bootsrapping. (Default 0.85)
#' @param rrep_as_crit Whether DTU calls should include replicate reproducibility as a criterion. (Default FALSE)
#' @param description Free-text description of the run. You can use this to add metadata to the results object.
#' @param verbose Display progress updates and warnings. (Default \code{TRUE})
#' @param threads Number of threads to use. (Default 1) Multi-threading will be ignored on non-POSIX systems.
#' @param dbg Debugging mode. Interrupt execution at the specified flag-point. Used to speed up code-tests by avoiding irrelevant downstream processing. (Default 0: do not interrupt)
#' @return List of mixed types. Contains a list of runtime settings, a table of gene-level results, a table of transcript-level results, and a list of two tables with the transcript abundaces.
#'
#' @import utils
#' @import parallel
#' @import data.table
#' @import matrixStats
#' @export
call_DTU <- function(annot= NULL, TARGET_COL= "target_id", PARENT_COL= "parent_id",
                     slo= NULL, name_A= "Condition-A", name_B= "Condition-B", varname= "condition", COUNTS_COL= "est_counts", BS_TARGET_COL= "target_id",
                     count_data_A = NULL, count_data_B = NULL, boot_data_A = NULL, boot_data_B = NULL,
                     p_thresh= 0.05, abund_thresh= 1, dprop_thresh= 0.2, correction= "BH", scaling= 1,
                     testmode= "both", qboot= TRUE, qbootnum= 0L, qrep_thresh= 0.95, rboot=FALSE, rrep_thresh= 0.85, rrep_as_crit= FALSE,
                     description= NA_character_, verbose= TRUE, threads= 1L, dbg= 0)
{
  #---------- PREP


  if (verbose) {
    message(paste0("Relative Abundance of Transcripts v.", packageVersion("rats")))
    message("Checking parameters...")
  }
  # Input checks.
  paramcheck <- parameters_are_good(slo, annot, name_A, name_B, varname, COUNTS_COL,
                                correction, p_thresh, TARGET_COL, PARENT_COL, BS_TARGET_COL, abund_thresh, testmode,
                                qboot, qbootnum, dprop_thresh, count_data_A, count_data_B, boot_data_A, boot_data_B,
                                qrep_thresh, threads, rboot, rrep_thresh, rrep_as_crit, scaling)
  if (paramcheck$error)
    stop(paramcheck$message)
  if (verbose)
    if(paramcheck$warn)
      for (w in paramcheck$warnings) {
        warning(w)  # So it displays as warning at the end of the run.
        message(w)  # So it displays at runtime.
      }

  threads <- as.integer(threads)  # Can't be decimal.
  # Ensure data.table complies.
  # if (packageVersion("data.table") >= "1.9.8")
    setDTthreads(threads)

  if (qbootnum == 0 && qboot)   # Use smart default.
    qbootnum = paramcheck$maxboots
  qbootnum <- as.integer(qbootnum)  # Can't be decimal.
  test_transc <- any(testmode == c("transc", "both"))
  test_genes <- any(testmode == c("genes", "both"))

  # Determine data extraction steps.
  steps <- 1  # Assume estimated counts. Simplest case. Bypasses both extraction and averaging.
  if (!is.null(boot_data_A)) {
    steps <- 2  # Bootstrapped estimates. Bypasses extraction, only needs averaging.
  } else if (!is.null(slo)) {
    steps <- 3  # Sleuth object. Most steps. Requires both extraction and averaging.
  }
  if (steps == 1 || is.na(qbootnum))
    qboot <- FALSE  # No quantification bootstraps data was provided.

  if (dbg == "prep")
    return(steps)


  #----------- LOOK-UP


  # Look-up from target_id to parent_id.
  if (verbose)
    message("Creating look-up structures...")
  tx_filter <- tidy_annot(annot, TARGET_COL, PARENT_COL)

  if (dbg == "indx")
    return(tx_filter)


  #---------- EXTRACT DATA

  
  if (verbose)
    message("Preparing data...")
  
  if (steps == 3) {   # From Sleuth
    # Reverse look-up from replicates to covariates.
  	samples_by_condition <- group_samples(slo$sample_to_covariates)[[varname]]
  	# Re-order rows and collate booted counts in a dataframe per sample. Put dataframes in a list per condition.
    # Target_id is included but NOT used as key so as to ensure that the order keeps matching tx_filter.
    boot_data_A <- denest_sleuth_boots(slo, tx_filter, samples_by_condition[[name_A]], COUNTS_COL, BS_TARGET_COL, "target_id", "parent_id", threads )
    boot_data_B <- denest_sleuth_boots(slo, tx_filter, samples_by_condition[[name_B]], COUNTS_COL, BS_TARGET_COL, "target_id", "parent_id", threads )
  } else if (steps == 2) {    # From generic bootstrapped data
    # Just re-order rows.
    boot_data_A <- lapply(boot_data_A, function(x) { x[match(tx_filter$target_id, x[[1]]), ] })
    boot_data_B <- lapply(boot_data_B, function(x) { x[match(tx_filter$target_id, x[[1]]), ] })
  }
  
  if (dbg == "bootin")
    return(list("bootA"=count_data_A, "bootB"=count_data_B))

  if (steps > 1) {    # Either Sleuth or generic bootstraps.
    # Remove ID columns so I don't have to always subset for math operations.
    for (smpl in c(boot_data_A, boot_data_B)) {
        smpl[, 1 := NULL]
    }
    
    # Scale abundances before aggregating.
    if (scaling != 1) {
      if (verbose)
        message("Applying scaling factor...")
      boot_data_A <- mclapply(boot_data_A, function(dt) { return(dt * scaling) }, mc.cores = threads, mc.preschedule = TRUE, mc.allow.recursive = FALSE)
      boot_data_B <- mclapply(boot_data_B, function(dt) { return(dt * scaling) }, mc.cores = threads, mc.preschedule = TRUE, mc.allow.recursive = FALSE)
    }
    
    # Calculate mean count across bootstraps.
    count_data_A <- as.data.table(mclapply(boot_data_A, function(b) { return(rowMeans(b)) },
                                           mc.cores = threads, mc.preschedule = TRUE, mc.allow.recursive = FALSE))
    count_data_B <- as.data.table(mclapply(boot_data_B, function(b) { return(rowMeans(b)) },
                                           mc.cores = threads, mc.preschedule = TRUE, mc.allow.recursive = FALSE))
  } else {    # From generic unbootstrapped data.
    # Just re-order rows and crop out the transcript IDs.
    nn <- names(count_data_A)
    count_data_A <- count_data_A[match(tx_filter$target_id, count_data_A[[1]]),  nn[seq.int(2, length(nn))],  with=FALSE]
    nn <- names(count_data_B)  # The number of columns may be different from A.
    count_data_B <- count_data_B[match(tx_filter$target_id, count_data_B[[1]]),  nn[seq.int(2, length(nn))],  with=FALSE]
    
    # Scale abundances.
    if (scaling != 1) {
      if (verbose)
        message("Applying scaling factor...")
      count_data_A <- count_data_A * scaling
      count_data_B <- count_data_B * scaling
    }
  }

  if (dbg == "data")
    return(list("countA"=count_data_A, "countB"=count_data_B))


  #---------- TEST


  # Do the core work.
  if (verbose)
    message("Calculating significances...")
  resobj <- calculate_DTU(count_data_A, count_data_B, tx_filter, test_transc, test_genes, "full", abund_thresh, p_thresh, dprop_thresh, correction, threads)

  if (dbg == "test")
    return(resobj)

  #-------- ADD INFO

  if (verbose)
    message("Filling in info and descriptive statistics...")

  with(resobj, {
    # Fill in results detais.
    Genes[, known_transc :=  Transcripts[, length(target_id), by=parent_id][, V1] ]  # V1 is the automatic column name for the lengths in the subsetted data.table
    Genes[, detect_transc :=  Transcripts[, .(parent_id, ifelse(sumA + sumB > 0, 1, 0))][, as.integer(sum(V2)), by = parent_id][, V1] ]  # Sum returns type double.
    Genes[(is.na(detect_transc)), detect_transc := 0]
    Transcripts[, meanA := rowMeans(count_data_A) ]
    Transcripts[, meanB := rowMeans(count_data_B) ]
    Transcripts[, stdevA := rowSds(as.matrix(count_data_A)) ]
    Transcripts[, stdevB := rowSds(as.matrix(count_data_B)) ]
    Transcripts[, log2FC := log2(sumB / sumA)]
  })

  # Fill in run info. (if done within the with() block, changes are local-scoped and don't take effect)
  resobj$Parameters["var_name"] <- varname
  resobj$Parameters["cond_A"] <- name_A
  resobj$Parameters["cond_B"] <- name_B
  resobj$Parameters["p_thresh"] <- p_thresh
  resobj$Parameters["abund_thresh"] <- abund_thresh
  resobj$Parameters["dprop_thresh"] <- dprop_thresh
  resobj$Parameters["abund_scaling"] <- scaling
  resobj$Parameters["tests"] <- testmode
  resobj$Parameters["rep_boot"] <- rboot
  resobj$Parameters["quant_boot"] <- qboot
  if (steps==3) {
    resobj$Parameters["data_type"] <- "sleuth"
  } else if (steps==2) {
    resobj$Parameters["data_type"] <- "bootstrapped abundance estimates"
  } else  if (steps==1) {
    resobj$Parameters["data_type"] <- "plain abundance estimates"
  }
  resobj$Parameters["num_genes"] <- length(levels(annot[[PARENT_COL]]))
  resobj$Parameters["num_transc"] <- length(annot[[TARGET_COL]])
  resobj$Parameters["description"] <- description

  if (dbg == "info")
    return(resobj)


  #---------- INTER-REPLICATE VARIABILITY BOOTSTRAP


  if (rboot) {
    if (verbose)
      message("Bootstrapping replicates...")

    pairs <- as.data.frame(t( expand.grid(1:dim(count_data_A)[2], 1:dim(count_data_B)[2]) ))
    numpairs <- length(pairs)
    resobj$Parameters["rep_bootnum"] <- numpairs
    resobj$Parameters["rep_reprod_thresh"] <- rrep_thresh
    resobj$Parameters["rep_reprod_as_crit"] <- rrep_as_crit

    if (verbose)
      myprogress <- utils::txtProgressBar(min = 0, max = numpairs, initial = 0, char = "=", width = NA, style = 3, file = "")

    repres <- lapply(1:numpairs, function(p) {  # Single-threaded. Forking happens within calculate_DTU().
                  # Update progress.
                  if (verbose)
                    setTxtProgressBar(myprogress, p)

                  # Grab a replicate from each condition.
                  counts_A <- as.data.table( count_data_A[[ names(count_data_A)[pairs[[p]][1]] ]] )
                  counts_B <- as.data.table( count_data_A[[ names(count_data_B)[pairs[[p]][2]] ]] )

                  # Do the work.
                  pout <- calculate_DTU(counts_A, counts_B, tx_filter, test_transc, test_genes, "short", abund_thresh, p_thresh, dprop_thresh, correction, threads)

                  with(pout, {
                    return(list("pp" = Transcripts[, pval_corr],
                                "pdtu" = Transcripts[, DTU],
                                "gp" = Genes[, pval_corr],
                                "gdtu" = Genes[, DTU] )) })
                })

    if (verbose)  # Forcing a new line after the progress bar.
      message("")

    if (dbg == "rboot")
      return(repres)

    #----- Stats

    if (verbose)
      message("Summarising bootstraps...")

    with(resobj, {
      if (test_transc) {
        eltr <- Transcripts[, elig]
        pd <- as.matrix(as.data.table(mclapply(repres, function(p) { p[["pdtu"]] }, mc.cores= threads, mc.preschedule = TRUE, mc.allow.recursive = FALSE)))
        Transcripts[(elig), rep_dtu_freq := rowCounts(pd[eltr, ], value = TRUE, na.rm=TRUE) / numpairs]
        pp <- as.matrix(as.data.table(mclapply(repres, function(p) { p[["pp"]] }, mc.cores= threads, mc.preschedule = TRUE, mc.allow.recursive = FALSE)))
        Transcripts[(elig), rep_p_mean := rowMeans(pp[eltr, ], na.rm = TRUE)]
        Transcripts[(elig), rep_p_stdev := rowSds(pp[eltr, ], na.rm = TRUE)]
        Transcripts[(elig), rep_p_min := rowMins(pp[eltr, ], na.rm = TRUE)]
        Transcripts[(elig), rep_p_max := rowMaxs(pp[eltr, ], na.rm = TRUE)]
        Transcripts[(elig), rep_na_freq := rowCounts(pp[eltr, ], value = NA, na.rm=FALSE) / numpairs]
        Transcripts[(elig & DTU), rep_reprod := (rep_dtu_freq >= rrep_thresh)]
        Transcripts[(elig & !DTU), rep_reprod := (rep_dtu_freq <= 1-rrep_thresh)]
      }

      if (test_genes) {
        elge <- Genes[, elig]
        gres <- as.matrix(as.data.table(mclapply(repres, function(p) { p[["gp"]] }, mc.cores= threads, mc.preschedule = TRUE, mc.allow.recursive = FALSE)))
        gdres <- as.matrix(as.data.table(mclapply(repres, function(p) { p[["gdtu"]] }, mc.cores= threads, mc.preschedule = TRUE, mc.allow.recursive = FALSE)))
        Genes[(elig), rep_dtu_freq := rowCounts(gdres[elge, ], value = TRUE, na.rm = TRUE) / numpairs]
        Genes[(elig), rep_p_mean := rowMeans(gres[elge, ], na.rm = TRUE)]
        Genes[(elig), rep_p_stdev := rowSds(gres[elge, ], na.rm = TRUE)]
        Genes[(elig), rep_p_min := rowMins(gres[elge, ], na.rm = TRUE)]
        Genes[(elig), rep_p_max := rowMaxs(gres[elge, ], na.rm = TRUE)]
        Genes[(elig), rep_na_freq := rowCounts(gres[elge, ], value = NA, na.rm = FALSE) / numpairs]  # It doesn't matter if AB or BA, affected identically by gene eligibility.
        Genes[(elig & DTU), rep_reprod := (rep_dtu_freq >= rrep_thresh)]
        Genes[(elig & !DTU), rep_reprod := (rep_dtu_freq <= 1-rrep_thresh)]
      }
    })

    if (dbg == "rbootsum")
      return(resobj)
  }


  #---------- QUANTIFICATION BOOTSTRAP


  if (qboot) {
    if (verbose) {
      message("Bootstrapping quantifications...")
      # Bootstrapping can take long, so showing progress is nice.
      myprogress <- utils::txtProgressBar(min = 0, max = qbootnum, initial = 0, char = "=", width = NA, style = 3, file = "")
    }

    resobj$Parameters["quant_reprod_thresh"] <- qrep_thresh
    resobj$Parameters["quant_bootnum"] <- qbootnum

    #----- Iterations

    bootres <- lapply(1:qbootnum, function(b) {  # Single-threaded. Forking happens within calculate_DTU().
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
                  bout <- calculate_DTU(counts_A, counts_B, tx_filter, test_transc, test_genes, "short", abund_thresh, p_thresh, dprop_thresh, correction, threads)

                  with(bout, {
                    return(list("pp" = Transcripts[, pval_corr],
                                "pdtu" = Transcripts[, DTU],
                                "gp" = Genes[, pval_corr],
                                "gdtu" = Genes[, DTU] )) })
              })
    if (verbose)  # Forcing a new line after the progress bar.
      message("")

    if (dbg == "qboot")
      return(bootres)

    #----- Stats

    if (verbose)
      message("Summarising bootstraps...")

    with(resobj, {
      if (test_transc) {
        # !!! POSSIBLE source of ERRORS if bootstraps * transcripts exceed R's maximum matrix size. (due to number of either) !!!
        pd <- as.matrix(as.data.table(mclapply(bootres, function(b) { b[["pdtu"]] }, mc.cores= threads, mc.preschedule = TRUE, mc.allow.recursive = FALSE)))
        Transcripts[(elig), quant_dtu_freq := rowCounts(pd[Transcripts[, elig], ], value = TRUE, na.rm=TRUE) / qbootnum]
        pp <- as.matrix(as.data.table(mclapply(bootres, function(b) { b[["pp"]] }, mc.cores= threads, mc.preschedule = TRUE, mc.allow.recursive = FALSE)))
        Transcripts[(elig), quant_p_mean := rowMeans(pp[Transcripts[, elig], ], na.rm = TRUE)]
        Transcripts[(elig), quant_p_stdev := rowSds(pp[Transcripts[, elig], ], na.rm = TRUE)]
        Transcripts[(elig), quant_p_min := rowMins(pp[Transcripts[, elig], ], na.rm = TRUE)]
        Transcripts[(elig), quant_p_max := rowMaxs(pp[Transcripts[, elig], ], na.rm = TRUE)]
        Transcripts[(elig), quant_na_freq := rowCounts(pp[Transcripts[, elig], ], value = NA, na.rm=FALSE) / qbootnum]
        Transcripts[(elig & DTU), quant_reprod := (quant_dtu_freq >= qrep_thresh)]
        Transcripts[(elig & !DTU), quant_reprod := (quant_dtu_freq <= 1-qrep_thresh)]
      }
      if (test_genes) {
        # !!! POSSIBLE source of ERRORS if bootstraps * genes exceed R's maximum matrix size. (due to number of bootstraps) !!!
        gres <- as.matrix(as.data.table(mclapply(bootres, function(b) { b[["gp"]] }, mc.cores= threads, mc.preschedule = TRUE, mc.allow.recursive = FALSE)))
        gdres <- as.matrix(as.data.table(mclapply(bootres, function(b) { b[["gdtu"]] }, mc.cores= threads, mc.preschedule = TRUE, mc.allow.recursive = FALSE)))
        Genes[(elig), quant_dtu_freq := rowCounts(gdres[Genes[, elig], ], value = TRUE, na.rm = TRUE) / qbootnum]
        Genes[(elig), quant_p_mean := rowMeans(gres[Genes[, elig], ], na.rm = TRUE)]
        Genes[(elig), quant_p_stdev := rowSds(gres[Genes[, elig], ], na.rm = TRUE)]
        Genes[(elig), quant_p_min := rowMins(gres[Genes[, elig], ], na.rm = TRUE)]
        Genes[(elig), quant_p_max := rowMaxs(gres[Genes[, elig], ], na.rm = TRUE)]
        Genes[(elig), quant_na_freq := rowCounts(gres[Genes[, elig], ], value = NA, na.rm = FALSE) / qbootnum]  # It doesn't matter if AB or BA, affected identically by gene eligibility.
        Genes[(elig & DTU), quant_reprod := (quant_dtu_freq >= qrep_thresh)]
        Genes[(elig & !DTU), quant_reprod := (quant_dtu_freq <= 1-qrep_thresh)]
      }
    })
  }

  if (dbg == "qbootsum")
    return(resobj)


  #---------- DONE


  if (verbose)
    message("Tidying up...")

  # Reject low-reproducibility DTU calls.
  with(resobj, {
    if (qboot) {
      Transcripts[(elig), DTU := (DTU & quant_reprod)]
      Genes[(elig), DTU := (DTU & quant_reprod)]
    }
    if (rrep_as_crit & rboot) {
      Transcripts[(elig), DTU := (DTU & rep_reprod)]
      Genes[(elig), DTU := (DTU & rep_reprod)]
    }
  })

  # Store the replicate means after re-adding the IDs.
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
  resobj$Abundances <- list("condA"= count_data_A, "condB"= count_data_B)

  with(resobj, {
    # Cross-display the DTU calls.
    Genes[, maxDprop := Transcripts[, maxabs(Dprop), by=parent_id][, V1] ]
    if (test_transc)
      Genes[, transc_DTU := Transcripts[, any(DTU), by = parent_id][, V1] ]
    if (test_genes)
      Transcripts[, gene_DTU := merge(Genes[, .(parent_id, DTU)], Transcripts[, .(parent_id)])[, DTU] ]

    # Drop the bootstrap columns, if unused.
    if (!qboot || !test_transc)
        Transcripts[, c("quant_dtu_freq", "quant_p_mean", "quant_p_stdev", "quant_p_min", "quant_p_max", "quant_na_freq", "quant_reprod") := NULL]
    if(!qboot || !test_genes)
      Genes[, c("quant_dtu_freq", "quant_p_mean", "quant_p_stdev", "quant_p_min", "quant_p_max", "quant_na_freq", "quant_reprod") := NULL]
    if (!rboot || !test_transc)
      Transcripts[, c("rep_dtu_freq", "rep_p_mean", "rep_p_stdev", "rep_p_min", "rep_p_max", "rep_na_freq", "rep_reprod") := NULL]
    if(!rboot || !test_genes)
      Genes[, c("rep_dtu_freq", "rep_p_mean", "rep_p_stdev", "rep_p_min", "rep_p_max", "rep_na_freq", "rep_reprod") := NULL]
  })

  if(verbose) {
    message("All done!")
    message("Summary of DTU results:")
    dtusum <- dtu_summary(resobj)
    print(dtusum)
    message("Isoform-switching subset of DTU:")
    switchsum <- dtu_switch_summary(resobj)
    print(switchsum)
    
  }

  return(resobj)
}


