#================================================================================
#' Calculate differential transcript usage.
#'
#' There are two modes for input:
#' \itemize{
#  \item{A sleuth object. This requires the following parameters: \code{slo}, \code{name_A}, \code{name_B} and optionally \code{varname}, \code{COUNT_COL}, \code{BS_TARGET_COL}.}
#'  \item{Bootstrapped count estimates. This requires the following parameters: \code{boot_data_A} and \code{boot_data_B}.}
#'  \item{Count estimates. This requires the following parameters: \code{count_data_A} and \code{count_data_B}.}
#' }
#'
#' @param annot A data.table matching transcript identifiers to gene identifiers. Any additional columns are allowed but ignored.
#' @param count_data_A A data.table of estimated counts for condition A. One column per sample/replicate, one row per transcript. The first column should contain the transcript identifiers.
#' @param count_data_B A data.table of estimated counts for condition B. One column per sample/replicate, one row per transcript. The first column should contain the transcript identifiers.
#' @param boot_data_A A list of data.tables, one per sample/replicate of condition A. One bootstrap iteration's estimates per column, one transcript per row. The first column should contain the transcript identifiers.
#' @param boot_data_B A list of data.tables, one per sample/replicate of condition B. One bootstrap iteration's estimates per column, one transcript per row. The first column should contain the transcript identifiers.
#' @param TARGET_COL The name of the column for the transcript identifiers in \code{annot}. (Default \code{"target_id"})
#' @param PARENT_COL The name of the column for the gene identifiers in \code{annot}. (Default \code{"parent_id"})
#' @param name_A The name for one condition. (Default "Condition-A")
#' @param name_B The name for the other condition. (Default "Condition-B")
#' @param varname The name of the covariate to which the two conditions belong. (Default \code{"condition"}).
#' @param p_thresh The p-value threshold. (Default 0.05)
#' @param dprop_thresh Effect size threshold. Minimum change in proportion of a transcript for it to be considered meaningful. (Default 0.20)
#' @param abund_thresh Noise threshold. Minimum mean (across replicates) abundance for transcripts (and genes) to be eligible for testing. (Default 5)
#' @param correction The p-value correction to apply, as defined in  \code{\link[stats]{p.adjust.methods}}. (Default \code{"BH"})
#' @param scaling A scaling factor or vector of scaling factors, to be applied to the abundances *prior* to any thresholding and testing. Useful for scaling TPMs (transcripts per 1 million reads) to the actual library sizes of the samples. If a vector is supplied, the order should correspond to the samples in group A followed by the samples in group B. WARNING: Improper use of the scaling factor will artificially inflate/deflate the significances obtained.
#' @param testmode One of \itemize{\item{"genes"}, \item{"transc"}, \item{"both" (default)}}.
#' @param lean Reduce memory footprint by not tracking mean/median/max/min/stdev for Dprop and pval across bootstrap iterations. The respective columns will be absent from the output structure. (Default TRUE)
#' @param qboot Bootstrap the DTU robustness against bootstrapped quantifications data. (Default \code{TRUE}) Ignored if input is \code{count_data}.
#' @param qbootnum Number of iterations for \code{qboot}. (Default 0) If 0, RATs will try to infer a value from the data.
#' @param qrep_thresh Reproducibility threshold for quantification bootsrapping. (Default 0.95)
#' @param rboot Bootstrap the DTU robustness against the replicates. Does *ALL* 1 vs 1 combinations. (Default \code{TRUE})
#' @param rrep_thresh Reproducibility threshold for replicate bootsrapping. (Default 0.85) With few replicates per condition, the reproducibility takes heavily quantized values. For 3x3, there are 9 possible 1v1 comparisons, and a consistency of 8/9 = 0.88.
#' @param description Free-text description of the run. You can use this to add metadata to the results object.
#' @param verbose Display progress updates and warnings. (Default \code{TRUE})
#' @param threads Number of threads to use. (Default 1) Multi-threading will be ignored on non-POSIX systems.
#' @param seed A numeric integer used to initialise the random number engine. Use this only if reproducible bootstrap selections are required. (Default NA)
#' @param reckless RATs normally aborts if any inconsistency is detected among the transcript IDs found in the annotation and the quantifications. Enabling reckless mode will downgrade this error to a warning and allow RATs to continue the run. Not recommended unless you know why the inconsistency exists and how it will affect the results. (Default FALSE)
#' @param dbg Debugging mode only. Interrupt execution at the specified flag-point. Used to speed up code-tests by avoiding irrelevant downstream processing. (Default 0: do not interrupt)
#' @return List of mixed types. Contains a list of runtime settings, a table of gene-level results, a table of transcript-level results, and a list of two tables with the transcript abundaces.
#'
#' @import utils
#' @import parallel
#' @import data.table
#' @import matrixStats
#' @export
call_DTU <- function(annot= NULL, TARGET_COL= "target_id", PARENT_COL= "parent_id",
                     count_data_A= NULL, count_data_B= NULL, boot_data_A= NULL, boot_data_B= NULL,
                     name_A= "Condition-A", name_B= "Condition-B", varname= "condition",
                     p_thresh= 0.05, abund_thresh= 5, dprop_thresh= 0.2, correction= "BH", scaling= 1.0,
                     testmode= "both", lean= TRUE, qboot= TRUE, qbootnum= 0L, qrep_thresh= 0.95, rboot=TRUE, rrep_thresh= 0.85,
                     description= NA_character_, verbose= TRUE, threads= 1L, seed= NA_integer_, reckless= FALSE, dbg= "0") {
  
  
  #---------- PARAMETERS

    
  {
    # Input checks.
    paramcheck <- parameters_are_good(annot, count_data_A, count_data_B, boot_data_A, boot_data_B,
                                      TARGET_COL, PARENT_COL,
                                      correction, testmode, scaling, threads, seed,
                                      p_thresh, abund_thresh, dprop_thresh, 
                                      qboot, qbootnum, qrep_thresh, rboot, rrep_thresh, reckless, lean, verbose)
    if (paramcheck$error)
      stop(paramcheck$message)
    
    if (verbose) {
      message(paste0("Relative Abundance of Transcripts v.", packageVersion("rats")))
      message("Preparing...")
      
      if(paramcheck$warn) {
        for (w in paramcheck$warnings) {
          warning(w)  # So it displays as warning at the end of the run.
          message(w)  # So it displays at runtime.
        }
      }
    }
    if (!is.na(seed))
      set.seed(as.integer(seed))
    threads <- as.integer(threads)          # Can't be decimal.
    setDTthreads(threads)
    if (qbootnum == 0 && qboot)             # Use smart default.
      qbootnum <- paramcheck$maxboots
    qbootnum <- as.integer(qbootnum)        # Can't be decimal.
    test_transc <- any(testmode == c("transc", "both"))
    test_genes <- any(testmode == c("genes", "both"))
  
    # Determine type of data input to handle accordingly.
    steps <- 1  # Assume plain estimated counts. Simplest case.
    if (!is.null(boot_data_A))
      steps <- 2  # Bootstrapped estimates. 
    # Make sure no mixed signals.
    if (steps == 1 || is.na(qbootnum) || qbootnum==0)
      qboot <- FALSE
  
    if (dbg == "prep")
      return(list(steps, qbootnum, test_transc, test_genes))
  }

  
  #----------- ANNOTATION


  {
    if (verbose)
      message("Creating look-up structure...")
    
    # Look-up from target_id to parent_id.
    tx_filter <- tidy_annot(annot, TARGET_COL, PARENT_COL)
  
    if (dbg == "indx")
      return(tx_filter)
  }

  
  #---------- DATA


  {
    if (verbose)
      message("Preparing data...")
    
    # Preprocess bootstrapped data.
    if (steps == 2) {
      # Re-order rows to match the annotation.
      boot_data_A <- lapply(boot_data_A, function(x) { x[match(tx_filter$target_id, x[[1]]), ] })
      boot_data_B <- lapply(boot_data_B, function(x) { x[match(tx_filter$target_id, x[[1]]), ] })
      # Remove ID columns so I don't have to always subset for math operations.
      for (smpl in c(boot_data_A, boot_data_B))
          smpl[, 1 := NULL]
      if (dbg == "bootin")
        return(list("bootA"=boot_data_A, "bootB"=boot_data_B))
      # Calculate mean count across bootstraps.
      count_data_A <- as.data.table(mclapply(boot_data_A, function(b) { return(rowMeans(b)) },
                                             mc.cores = threads, mc.preschedule = TRUE, mc.allow.recursive = FALSE))
      count_data_B <- as.data.table(mclapply(boot_data_B, function(b) { return(rowMeans(b)) },
                                             mc.cores = threads, mc.preschedule = TRUE, mc.allow.recursive = FALSE))
    # Preprocess unbootstrapped data.
    } else {
      # Just re-order rows and crop out the transcript IDs.
      nn <- names(count_data_A)
      count_data_A <- count_data_A[match(tx_filter$target_id, count_data_A[[1]]),  nn[seq.int(2, length(nn))],  with=FALSE]
      nn <- names(count_data_B)  # The number of columns may be different from A.
      count_data_B <- count_data_B[match(tx_filter$target_id, count_data_B[[1]]),  nn[seq.int(2, length(nn))],  with=FALSE]
    }
  
    if (dbg == "countin")
      return(list("countA"=count_data_A, "countB"=count_data_B))
    
      
    if (any(scaling != 1)) {
      if (verbose)
        message("Applying scaling factor(s)...")
      
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
    }
    
    if (dbg == "scale")
      return(list("countA"=count_data_A, "countB"=count_data_B, "bootA"=boot_data_A, "bootB"=boot_data_B))
  }

  
  #---------- TEST


  {
    if (verbose)
      message("Calculating significances...")
      
    # Plain test on the abundances.
    resobj <- calculate_DTU(count_data_A, count_data_B, tx_filter, test_transc, test_genes, "full", abund_thresh, p_thresh, dprop_thresh, correction, threads)
  
    if (dbg == "test")
      return(resobj)
  }
  
  
  #-------- STATS & META-DATA

  
  {
    if (verbose)
      message("Filling in info...")
  
    with(resobj, {
      # Fill in data stats.
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
    resobj$Parameters[["var_name"]] <- varname
    resobj$Parameters[["cond_A"]] <- name_A
    resobj$Parameters[["cond_B"]] <- name_B
    resobj$Parameters[["p_thresh"]] <- p_thresh
    resobj$Parameters[["abund_thresh"]] <- abund_thresh
    resobj$Parameters[["dprop_thresh"]] <- dprop_thresh
    resobj$Parameters[["abund_scaling"]] <- c(scaling)
    resobj$Parameters[["tests"]] <- testmode
    resobj$Parameters[["rep_boot"]] <- rboot
    resobj$Parameters[["quant_boot"]] <- qboot
    if (steps==2) {
      resobj$Parameters[["data_type"]] <- "bootstrapped abundance estimates"
    } else  if (steps==1) {
      resobj$Parameters[["data_type"]] <- "plain abundance estimates"
    }
    resobj$Parameters[["num_genes"]] <- length(unique(annot[[PARENT_COL]]))
    resobj$Parameters[["num_transc"]] <- length(annot[[TARGET_COL]])
    resobj$Parameters[["description"]] <- description
    resobj$Parameters[["seed"]] <- seed
    resobj$Parameters[["correction"]] <- correction
    resobj$Parameters[["reckless"]] <- reckless
    resobj$Parameters[["lean"]] <- lean
    
    
    if (dbg == "info")
      return(resobj)
  }
  

  #---------- INTER-REPLICATE VARIABILITY BOOTSTRAP
  
  
  if (rboot) {
    if (verbose)
      message("Bootstrapping against replicates...")
    
    resobj$Parameters[["rep_bootnum"]] <- ncol(count_data_A) * ncol(count_data_B)
    resobj$Parameters[["rep_reprod_thresh"]] <- rrep_thresh
    # Avoid repetitive cluttered syntax.
    eltr <- resobj$Transcripts$elig
    elge <- resobj$Genes$elig
    # Bootstrap.
    tallies <- do_boot(lean= lean, tx_filter= tx_filter, eltr= eltr, elge= elge,
                       count_data_A= count_data_A, count_data_B= count_data_B,  
                       niter= resobj$Parameters[["rep_bootnum"]], threads= threads, verbose=TRUE, test_transc= test_transc, test_genes= test_genes, 
                       count_thresh= abund_thresh, p_thresh= p_thresh, dprop_thresh= dprop_thresh, correction= correction)
    
    if (dbg == "rboot")
      return(tallies)
    
    # Copy results into the output object.
    with(resobj, {
      Transcripts[(eltr), rep_dtu_freq := tallies$Transcripts$dtu_freq[eltr]]
      Transcripts[(eltr), rep_na_freq := tallies$Transcripts$na_freq[eltr]]
      Transcripts[(eltr & DTU), rep_reprod := (rep_dtu_freq >= rrep_thresh)]
      Transcripts[(eltr & !DTU), quant_reprod := (rep_dtu_freq < 1 - rrep_thresh)]
      Genes[(elge), rep_dtu_freq := tallies$Genes$dtu_freq[elge]]
      Genes[(elge), rep_na_freq := tallies$Genes$na_freq[elge]]
      Genes[(elge & DTU), rep_reprod := (rep_dtu_freq >= rrep_thresh)]
      Genes[(elge & !DTU), rep_reprod := (rep_dtu_freq < 1 - rrep_thresh)]
      if (!lean) {
        Transcripts[(eltr), rep_p_median := tallies$Transcripts$p_median[eltr]]
        Transcripts[(eltr), rep_p_min := tallies$Transcripts$p_min[eltr]]
        Transcripts[(eltr), rep_p_max := tallies$Transcripts$p_max[eltr]]
        Transcripts[(eltr), rep_Dprop_mean := tallies$Transcripts$Dprop_mean[eltr]]
        Transcripts[(eltr), rep_Dprop_stdev := tallies$Transcripts$Dprop_stdev[eltr]]
        Transcripts[(eltr), rep_Dprop_min := tallies$Transcripts$Dprop_min[eltr]]
        Transcripts[(eltr), rep_Dprop_max := tallies$Transcripts$Dprop_max[eltr]]
        Genes[(elge), rep_p_median := tallies$Genes$p_median[elge]]
        Genes[(elge), rep_p_min := tallies$Genes$p_min[elge]]
        Genes[(elge), rep_p_max := tallies$Genes$p_max[elge]]
      }
    })
    
    if (verbose)  # Forcing a new line after the progress bar.
      message("")
    
    if (dbg == "rbootsum")
      return(resobj)
  }
  
  
  #---------- QUANTIFICATION BOOTSTRAP
  
  
  if (qboot) {
    if (verbose)
      message("Bootstrapping against quantifications...")
    
    resobj$Parameters[["quant_bootnum"]] <- qbootnum
    resobj$Parameters[["quant_reprod_thresh"]] <- qrep_thresh
    # Avoid repetitive cluttered syntax.
    eltr <- resobj$Transcripts$elig
    elge <- resobj$Genes$elig
    # Bootstrap.
    tallies <- do_boot(lean- lean, tx_filter= tx_filter, eltr= eltr, elge= elge,
                       boot_data_A= boot_data_A, boot_data_B= boot_data_B,  
                       niter= qbootnum, threads= threads, verbose= verbose, test_transc= test_transc, test_genes= test_genes, 
                       count_thresh= abund_thresh, p_thresh= p_thresh, dprop_thresh= dprop_thresh, correction= correction)
    
    if (dbg == "qboot")
      return(tallies)
    
    # Copy results into the output object.
    with(resobj, {
      Transcripts[(eltr), quant_dtu_freq := tallies$Transcripts$dtu_freq[eltr]]
      Transcripts[(eltr), quant_na_freq := tallies$Transcripts$na_freq[eltr]]
      Transcripts[(eltr & DTU), quant_reprod := (quant_dtu_freq >= qrep_thresh)]
      Transcripts[(eltr & !DTU), quant_reprod := (quant_dtu_freq < 1 - qrep_thresh)]
      Genes[(elge), quant_dtu_freq := tallies$Genes$dtu_freq[elge]]
      Genes[(elge), quant_na_freq := tallies$Genes$na_freq[elge]]
      Genes[(elge & DTU), quant_reprod := (quant_dtu_freq >= qrep_thresh)]
      Genes[(elge & !DTU), quant_reprod := (quant_dtu_freq < 1 - qrep_thresh)]
      if (!lean) {
          Transcripts[(eltr), quant_p_median := tallies$Transcripts$p_median[eltr]]
          Transcripts[(eltr), quant_p_min := tallies$Transcripts$p_min[eltr]]
          Transcripts[(eltr), quant_p_max := tallies$Transcripts$p_max[eltr]]
          Transcripts[(eltr), quant_Dprop_mean := tallies$Transcripts$Dprop_mean[eltr]]
          Transcripts[(eltr), quant_Dprop_stdev := tallies$Transcripts$Dprop_stdev[eltr]]
          Transcripts[(eltr), quant_Dprop_min := tallies$Transcripts$Dprop_min[eltr]]
          Transcripts[(eltr), quant_Dprop_max := tallies$Transcripts$Dprop_max[eltr]]
          Genes[(elge), quant_p_median := tallies$Genes$p_median[elge]]
          Genes[(elge), quant_p_min := tallies$Genes$p_min[elge]]
          Genes[(elge), quant_p_max := tallies$Genes$p_max[elge]]
      }
    })
    
    if (verbose)  # Forcing a new line after the progress bar.
      message("")
    
    if (dbg == "qbootsum")
      return(resobj)
  }
  

  #---------- DONE


  {
    if (verbose)
      message("Tidying up...")
  
    # Reject low-reproducibility DTU calls.
    with(resobj, {
      if(qboot){
        Transcripts[(elig), DTU := (DTU & quant_reprod)]
        Genes[(elig), DTU := (DTU & quant_reprod)]
      }
      if (rboot) {
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
    resobj$Abundances <- list("condA" = count_data_A, "condB" = count_data_B)
  
    # Tidy up fields.
    with(resobj, {
      # Cross-display the DTU calls.
      Genes[, maxDprop := Transcripts[, maxabs(Dprop), by = parent_id][, V1] ]
      Genes[, transc_DTU := Transcripts[, any(DTU), by = parent_id][, V1] ]
      Transcripts[, gene_DTU := merge(Genes[, .(parent_id, DTU)], Transcripts[, .(parent_id)])[, DTU] ]
  
      # Eradicate NAs from flag columns, as they mess up sub-setting of the tables.
      Genes[is.na(elig), elig := c(FALSE)]
      Genes[is.na(sig), sig := c(FALSE)]
      Genes[is.na(elig_fx), elig_fx := c(FALSE)]
      Genes[is.na(quant_reprod), quant_reprod := c(FALSE)]
      Genes[is.na(rep_reprod), rep_reprod := c(FALSE)]
      Genes[is.na(DTU), DTU := c(FALSE)]
      Genes[is.na(transc_DTU), transc_DTU := c(FALSE)]
      Transcripts[is.na(elig_xp), elig_xp := c(FALSE)]
      Transcripts[is.na(elig), elig := c(FALSE)]
      Transcripts[is.na(sig), sig := c(FALSE)]
      Transcripts[is.na(elig_fx), elig_fx := c(FALSE)]
      Transcripts[is.na(quant_reprod), quant_reprod := c(FALSE)]
      Transcripts[is.na(rep_reprod), rep_reprod := c(FALSE)]
      Transcripts[is.na(DTU), DTU := c(FALSE)]
      Transcripts[is.na(gene_DTU), gene_DTU := c(FALSE)]
      
      # Eradicate NAs from count columns, as they mess up plotting. These would be caused by inconsistencies between annotation and quantifications, and would only exist in "reckless" mode.
      for (i in 1:Parameters$num_replic_A){
        idx <- is.na(Abundances[[1]][[paste0('V',i)]])
        Abundances[[1]][idx, paste0('V',i) := 0]
      }
      for (i in 1:Parameters$num_replic_B){
        idx <- is.na(Abundances[[2]][[paste0('V',i)]])
        Abundances[[2]][idx, paste0('V',i) := 0]
      }
      Transcripts[is.na(meanA), meanA := c(0)]
      Transcripts[is.na(meanB), meanB := c(0)]
      Transcripts[is.na(propA), propA := c(0)]
      Transcripts[is.na(propB), propB := c(0)]
      
      # Drop unused bootstrap columns.
      if (!qboot || !test_transc)
        Transcripts[, c("quant_dtu_freq", "quant_na_freq", "quant_reprod") := NULL]
      if(!qboot || !test_genes)
        Genes[, c("quant_dtu_freq", "quant_na_freq", "quant_reprod") := NULL]
      if (!rboot || !test_transc)
        Transcripts[, c("rep_dtu_freq", "rep_na_freq", "rep_reprod") := NULL]
      if(!rboot || !test_genes)
        Genes[, c("rep_dtu_freq", "rep_na_freq", "rep_reprod") := NULL]
      
      # Drop unused bootstrap statistics columns.
      if (lean || !qboot || !test_transc)
        Transcripts[, c("quant_p_median", "quant_p_min", "quant_p_max", "quant_Dprop_mean", "quant_Dprop_stdev", "quant_Dprop_min", "quant_Dprop_max") := NULL]
      if(lean || !qboot || !test_genes)
        Genes[, c( "quant_p_median", "quant_p_min", "quant_p_max") := NULL]
      if (lean || !rboot || !test_transc)
        Transcripts[, c("rep_p_median", "rep_p_min", "rep_p_max", "rep_Dprop_mean", "rep_Dprop_stdev", "rep_Dprop_min", "rep_Dprop_max") := NULL]
      if(lean || !rboot || !test_genes)
        Genes[, c("rep_p_median", "rep_p_min", "rep_p_max") := NULL]
      
    })
  
    if(verbose) {
      message("All done!")
      cat(noquote("# Summary of DTU results:\n"))
      print( dtu_summary(resobj) )
      cat(noquote("\n# Isoform-switching subset of DTU:\n"))
      print( dtu_switch_summary(resobj) )
    }
  }

  return(resobj)
}

#EOF
