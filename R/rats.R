#================================================================================
#' Calculate differential transcript usage.
#'
#' @param slo A sleuth object.
#' @param annot A dataframe matching the transcript identifiers to their corresponding gene identifiers.
#' @param name_A The name for one condition, as it appears in the \code{sample_to_covariates} table within the sleuth object.
#' @param name_B The name for the other condition, as it appears in the \code{sample_to_covariates} table within the sleuth object.
#' @param varname The name of the covariate to which the two conditions belong, as it appears in the \code{sample_to_covariates} table within the sleuth object. Default \code{"condition"}.
#' @param verbose Display progress updates, default \code{FALSE}.
#' @param p_thresh The p-value threshold, default 0.05.
#' @param count_thresh Minimum count of fragments per sample, in at least one of the conditions, for transcripts to be eligible for testing (default 5).
#' @param dprop_thresh Minimum change in proportion (effect size) of a transcript for it to be eligible to be significant. (default 0.1).
#' @param testmode One of "G-test", "prop-test", "both" (default both).
#' @param correction The p-value correction to apply, as defined in \code{stats::p.adjust.methods}, default \code{"BH"}.
#' @param boots Bootstrap the p-values of either test. One of "G-test", "prop-test", "both", "none". (default "none").
#' @param bootnum Number of bootstraps (default 10000).
#' @param COUNTS_COL The name of the counts column to use for the DTU calculation (est_counts or tpm), default \code{"est_counts"}.
#' @param TARGET_COL The name of the transcript identifier column in the annot object, default \code{"target_id"}
#' @param PARENT_COL The name of the parent identifier column in the annot object, default \code{"parent_id"}.
#' @param BS_TARGET_COL The name of the transcript identifier column in the sleuth bootstrap tables, default \code{"target_id"}.
#' @return List of data tables, with gene-level and transcript-level information.
#'
#' @import data.table
#' @import matrixStats
#' @export
call_DTU <- function(slo, annot, name_A, name_B, varname= "condition", 
                          p_thresh= 0.05, count_thresh= 5, dprop_thresh= 0.1, testmode= "both", correction= "BH", 
                          verbose= FALSE, boots= "none", bootnum= 100L,
                          COUNTS_COL= "est_counts", TARGET_COL= "target_id", PARENT_COL= "parent_id", BS_TARGET_COL= "target_id")
{
  #---------- PREP
  
  # Input checks.
  paramcheck <- parameters_good(slo, annot, name_A, name_B, varname, COUNTS_COL,
                                correction, p_thresh, TARGET_COL, PARENT_COL, BS_TARGET_COL, verbose, count_thresh, testmode, 
                                boots, bootnum, dprop_thresh)
  if (paramcheck$error) stop(paramcheck$message)
  bootnum <- as.integer(bootnum)
  
  boot_transc <- any(boots == c("prop-test", "both"))
  
  boot_genes <- any(boots == c("G-test", "g-test", "both"))
  
  # Set up progress bar
  progress <- init_progress(verbose)
  progress <- update_progress(progress)
  
  #----------- LOOK-UP
  
  progress <- update_progress(progress)
  # Look-up from target_id to parent_id.
  tx_filter <- data.table(target_id = annot[[TARGET_COL]], parent_id = annot[[PARENT_COL]])
  with( tx_filter,
        setkey(tx_filter, parent_id, target_id) )  # Fixates the order of genes and transcripts to be used throughout the rest of this package.
  # Reverse look-up from replicates to covariates.
  samples_by_condition <- group_samples(slo$sample_to_covariates)[[varname]]
  
  #---------- EXTRACT DATA
  
  progress <- update_progress(progress)
  # De-nest, average and index the counts from the bootstraps. 
  # For each condition, a list holds a dataframe for each replicate. 
  # Each dataframe holds the counts from ALL the bootstraps. Target_id is included but NOT used as key so as to ensure the order keeps matching tx_filter.
  data_A <- denest_boots(slo, tx_filter$target_id, samples_by_condition[[name_A]], COUNTS_COL, BS_TARGET_COL )
  data_B <- denest_boots(slo, tx_filter$target_id, samples_by_condition[[name_B]], COUNTS_COL, BS_TARGET_COL )
  bootmeans_A <- as.data.table(lapply(data_A, function(b) { n <- names(b); rowMeans(b[, n[1:length(n)-1], with=FALSE]) }))  # Tables don't have access to column ranges by index so I have to fish our the names.
  bootmeans_B <- as.data.table(lapply(data_B, function(b) { n <- names(b); rowMeans(b[, n[1:length(n)-1], with=FALSE]) }))
  
  #---------- TEST
  
  progress <- update_progress(progress)
  # Do the core work.
  suppressWarnings(
    resobj <- calculate_DTU(bootmeans_A, bootmeans_B, tx_filter, testmode, "full", count_thresh, p_thresh, dprop_thresh, correction) )
  
  #-------- ADD INFO
  
  progress <- update_progress(progress)
    # Fill in run info. (if done withing the with() block, changes don't take effect)
  resobj$Parameters["var_name"] <- varname
  resobj$Parameters["cond_A"] <- name_A
  resobj$Parameters["cond_B"] <- name_B
  resobj$Parameters["p_thresh"] <- p_thresh 
  resobj$Parameters["count_thresh"] <- count_thresh
  resobj$Parameters["dprop_thresh"] <- dprop_thresh
  resobj$Parameters["tests"] <- testmode
  resobj$Parameters["bootstrap"] <- boots
  resobj$Parameters["bootnum"] <- bootnum
  
  with(resobj, {
    # Fill in results detais.
    Genes[, known_transc :=  Transcripts[, length(target_id), by=parent_id][, V1] ]  # V1 is the automatic column name for the lengths in the subsetted data.table
    detected = Transcripts[, abs((propA + propB))] > 0
    Genes[, detect_transc :=  Transcripts[, .(parent_id, ifelse(abs(propA + propB) > 0, 1, 0))][, as.integer(sum(V2)), by = parent_id][, V1] ]  # Sum returns doube..
    Genes[(is.na(detect_transc)), detect_transc := 0]
    Transcripts[, meanA :=  rowMeans(bootmeans_A) ]
    Transcripts[, meanB :=  rowMeans(bootmeans_B) ]
    Transcripts[, stdevA :=  rowSds(as.matrix(bootmeans_A)) ]
    Transcripts[, stdevB :=  rowSds(as.matrix(bootmeans_B)) ]
    if (any(testmode == c("prop-test", "both")))
      Genes[, transc_DTU := Transcripts[, any(DTU), by = parent_id][, V1] ]
    if (any(testmode == c("G-test", "g-test", "both")))
      Transcripts[, gene_DTU := merge(Genes[, .(parent_id, DTU)], Transcripts[, .(parent_id)])[, DTU] ]
  })
  
  #---------- BOOTSTRAP
  
  progress <- update_progress(progress)
  if (any(boot_transc, boot_genes)) {
    
    #----- Iterations
    
    bootres <- lapply(1:bootnum, function(b) {
                  # Grab a bootstrap from each replicate. 
                  counts_A <- as.data.table(lapply(data_A, function(smpl) { smpl[[sample( names(smpl)[1:(dim(smpl)[2]-1)], 1)]] }))  # Have to use list syntax to get a vector back. 
                                                                                                                                     # Usual table syntax with "with=FALSE" returns a table and I fail to cast it.
                                                                                                                                     # Also, the last column is the target_id, so I leave it out.
                  counts_B <- as.data.table(lapply(data_B, function(smpl) { smpl[[sample( names(smpl)[1:(dim(smpl)[2]-1)], 1)]] }))
                  # Do the work.
                  # Ignore warning. Chi-square test generates warnings for counts <5. This is expected behaviour. Transcripts changing between off and on are often culprits.
                  suppressWarnings(
                    bout <- calculate_DTU(counts_A, counts_B, tx_filter, testmode, "short", count_thresh, p_thresh, dprop_thresh, correction))
                  
                  with(bout, {
                    return(list("pp" = Transcripts[, pval_corr],
                                "pdtu" = Transcripts[, DTU],
                                "gpab" = Genes[, pvalAB_corr],
                                "gpba" = Genes[, pvalBA_corr],
                                "gdtu" = Genes[, DTU] )) }) })
    
    #----- Stats
    
    with(resobj, {
      if (boot_transc) {
        # !!! POSSIBLE source of ERRORS if bootstraps * transcripts exceed R's maximum matrix size. (due to number of either) !!!
        pd <- as.matrix(as.data.table(lapply(bootres, function(b) { b[["pdtu"]] })))
        Transcripts[(elig), boot_freq := rowCounts(pd[Transcripts[, elig], ], value = TRUE) / bootnum]
        pp <- as.matrix(as.data.table(lapply(bootres, function(b) { b[["pp"]] })))
        Transcripts[(elig), boot_mean := rowMeans(pp[Transcripts[, elig], ], na.rm = TRUE)]
        Transcripts[(elig), boot_stdev := rowSds(pp[Transcripts[, elig], ], na.rm = TRUE)]
        Transcripts[(elig), boot_min := rowMins(pp[Transcripts[, elig], ], na.rm = TRUE)]
        Transcripts[(elig), boot_max := rowMaxs(pp[Transcripts[, elig], ], na.rm = TRUE)]
        Transcripts[(elig), boot_na := rowCounts(pp[Transcripts[, elig], ], value = NA) / bootnum]
      }
      if (boot_genes) {
        # !!! POSSIBLE source of ERRORS if bootstraps * genes exceed R's maximum matrix size. (due to number of bootstraps) !!!
        gabres <- as.matrix(as.data.table(lapply(bootres, function(b) { b[["gpab"]] })))
        gbares <- as.matrix(as.data.table(lapply(bootres, function(b) { b[["gpba"]] })))
        gdres <- as.matrix(as.data.table(lapply(bootres, function(b) { b[["gdtu"]] })))
        Genes[(elig), boot_freq := rowCounts(gdres[Genes[, elig], ], value = TRUE) / bootnum]
        Genes[(elig), boot_meanAB := rowMeans(gabres[Genes[, elig], ], na.rm = TRUE)]
        Genes[(elig), boot_meanBA := rowMeans(gbares[Genes[, elig], ], na.rm = TRUE)]
        Genes[(elig), boot_stdevAB := rowSds(gabres[Genes[, elig], ], na.rm = TRUE)]
        Genes[(elig), boot_stdevBA := rowSds(gbares[Genes[, elig], ], na.rm = TRUE)]
        Genes[(elig), boot_minAB := rowMins(gabres[Genes[, elig], ], na.rm = TRUE)]
        Genes[(elig), boot_minBA := rowMins(gbares[Genes[, elig], ], na.rm = TRUE)]
        Genes[(elig), boot_maxAB := rowMaxs(gabres[Genes[, elig], ], na.rm = TRUE)]
        Genes[(elig), boot_maxBA := rowMaxs(gbares[Genes[, elig], ], na.rm = TRUE)]
        Genes[(elig), boot_na := rowCounts(gabres[Genes[, elig], ], value = NA) / bootnum]  # It doesn't matter if AB or BA, affected identically by gene eligibility.
      }
    })
  }
  
  #---------- DONE
 
  with(resobj, {   
    # Drop columns that exist for convenient calculations but are not useful to end users.
    Transcripts[, c("totalA", "totalB") := NULL]  # Some plots need it but it's ugly in the report tables. It takes a blink to recalculate.
    # Drop the columns.
    if (!boot_transc)
      Transcripts[, c("boot_freq", "boot_mean", "boot_stdev", "boot_min", "boot_max", "boot_na") := NULL]
    if(!boot_genes)
      Genes[, c("boot_freq", "boot_meanAB", "boot_meanBA", "boot_stdevAB", "boot_stdevBA", "boot_minAB",
              "boot_minBA", "boot_maxAB", "boot_maxBA", "boot_na") := NULL]
  })
  
  progress <- update_progress(progress)

  if(verbose)
    print(dtu_summary(resobj))
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
#' @param verbose print progress report
#' @param count_thresh minimum frgments per transcript per sample
#' @param testmode which tests to run
#' @param boots which tests to bootstrap
#' @param bootnum number of bootstrap iterations
#' @param dprop_thresh minimum change in proportion
#' 
#' @return List with a logical value and a message.
#'
parameters_good <- function(slo, annot, name_A, name_B, varname, COUNTS_COL,
                            correction, p_thresh, TARGET_COL, PARENT_COL, BS_TARGET_COL, verbose, 
                            count_thresh, testmode, boots, bootnum, dprop_thresh) 
{
  if (!is.data.frame(annot))
    return(list("error"=TRUE, "message"="The provided annot is not a data.frame!"))
  
  if (any(! c(TARGET_COL, PARENT_COL) %in% names(annot)))
    return(list("error"=TRUE, "message"="The specified target and/or parent IDs field-names do not exist in annot!"))
  if (! BS_TARGET_COL %in% names(slo$kal[[1]]$bootstrap[[1]]))
    return(list("error"=TRUE, "message"="The specified target IDs field-name does not exist in the bootstraps!"))
  
  if (! COUNTS_COL %in% names(slo$kal[[1]]$bootstrap[[1]]))
    return(list("error"=TRUE, "message"="The specified counts field-name does not exist!"))
  
  if (! correction %in% p.adjust.methods)
    return(list("error"=TRUE, "message"="Invalid p-value correction method name. Refer to stats::p.adjust.methods!"))
  
  if ((!is.numeric(p_thresh)) || p_thresh > 1 || p_thresh < 0)
    return(list("error"=TRUE, "message"="Invalid p-value threshold!"))
  
  if (! varname %in% names(slo$sample_to_covariates))
    return(list("error"=TRUE, "message"="The specified covariate name does not exist!"))
  
  if (any(! c(name_A, name_B) %in% slo$sample_to_covariates[[varname]] ))
    return(list("error"=TRUE, "message"="One or both of the specified conditions do not exist!"))
  
  if (!is.logical(verbose))
    return(list("error"=TRUE, "message"="verbose must be a logical value!"))
  
  if (! testmode %in% c("G-test", "g-test", "prop-test", "both"))
    return(list("error"=TRUE, "message"="Unrecognized value for testmode!"))
  
  if (! boots %in% c("G-test", "g-test","prop-test", "both", "none"))
    return(list("error"=TRUE, "message"="Unrecognized value for boots!"))
  
  if ((!is.numeric(bootnum)) || bootnum < 1)
    return(list("error"=TRUE, "message"="Invalid number of bootstraps! Must be a positive number."))
  
  tx <- slo$kal[[1]]$bootstrap[[1]][[BS_TARGET_COL]][ order(slo$kal[[1]]$bootstrap[[1]][[BS_TARGET_COL]]) ]
  for (k in 2:length(slo$kal)) {
    if (!all( tx == slo$kal[[k]]$bootstrap[[1]][[BS_TARGET_COL]][ order(slo$kal[[k]]$bootstrap[[1]][[BS_TARGET_COL]]) ] ))
      return(list("error"=TRUE, "message"="Inconsistent set of transcript IDs! Please try again, using the same annotation for all your samples!"))
  }
  
  if ((!is.numeric(dprop_thresh)) || dprop_thresh < 0 || dprop_thresh > 1)
    return(list("error"=TRUE, "message"="Invalid proportion difference threshold! Must be between 0 and 1."))
  
  if ((!is.numeric(count_thresh)) || count_thresh < 0)
    return(list("error"=TRUE, "message"="Invalid read-count threshold! Must be between 0 and 1."))
  
  return(list("error"=FALSE, "message"="All good!"))
}


#================================================================================
#' Initialise progress updates
#'
#' @param on Flag indicating whether updates are on (TRUE) or not (FALSE)
#' @return The progress update object
#'
init_progress <- function(on)
{
  progress_steps <- data.frame(c(1, 2, 3, 10, 11, 15, 30, 31, 100),
                               c("Checked parameters...",
                                 "Creating look-up structures...",
                                 "Extracting counts from bootstraps...",
                                 "Allocating output structure...",
                                 "Calculating counts statistics...",
                                 "Calculating p-values...",
                                 "Filling in more info...",
                                 "Bootstrapping p-values (if applicable)...",
                                 "All done!"),
                               stringsAsFactors = FALSE)
  progress <- TxtProgressUpdate(steps=progress_steps, on=on)
  return(progress)
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
#' @return a list of data.tables, one per sample, containing all the bootstrap counts.
#'
#' Transcripts in \code{slo} that are missing from \code{tx} will be skipped completely.
#' Transcripts in \code{tx} that are missing from \code{slo} are automatically padded with NA, which we re-assign as 0.
#'
denest_boots <- function(slo, tx, samples, COUNTS_COL, BS_TARGET_COL) {
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
  })
}


#================================================================================
#' Create output structure.
#' 
#' @param annot a dataframe with at least 2 columns: target_id and parent_id.
#' @param full Full-sized structure or core fields only. One of "full", "short".
#' @return a list-based structure of class dtu.
#' 
alloc_out <- function(annot, full){
  if (full == "full") {
    Parameters <- list("var_name"=NA_character_, "cond_A"=NA_character_, "cond_B"=NA_character_,
                       "num_replic_A"=NA_integer_, "num_replic_B"=NA_integer_,
                       "p_thresh"=NA_real_, "count_thresh"=NA_real_, "dprop_thresh"=NA_real_,
                       "tests"=NA_character_, "bootstrap"=NA_character_, "bootnum"=NA_integer_)
    Genes <- data.table("parent_id"=levels(as.factor(annot$parent_id)),
                        "DTU"=NA, "transc_DTU"=NA, 
                        "known_transc"=NA_integer_, "detect_transc"=NA_integer_, "elig_transc"=NA_integer_,
                        "elig"=NA, "elig_fx"=NA,
                        "pvalAB"=NA_real_, "pvalBA"=NA_real_, "pvalAB_corr"=NA_real_, "pvalBA_corr"=NA_real_, "sig"=NA, 
                        "boot_freq"=NA_real_, "boot_meanAB"=NA_real_, "boot_meanBA"=NA_real_, 
                        "boot_stdevAB"=NA_real_, "boot_stdevBA"=NA_real_, "boot_minAB"=NA_real_, "boot_minBA"=NA_real_, 
                        "boot_maxAB"=NA_real_, "boot_maxBA"=NA_real_, "boot_na"=NA_real_)
    Transcripts <- data.table("target_id"=annot$target_id, "parent_id"=annot$parent_id,
                              "DTU"=NA, "gene_DTU"=NA,
                              "meanA"=NA_real_, "meanB"=NA_real_,    # mean across replicates of means across bootstraps
                              "stdevA"=NA_real_, "stdevB"=NA_real_,  # standard deviation across replicates of means across bootstraps
                              "sumA"=NA_real_, "sumB"=NA_real_, "elig_xp"=NA, "elig"=NA,    # sum across replicates of means across bootstraps
                              "propA"=NA_real_, "propB"=NA_real_, "Dprop"=NA_real_, "elig_fx"=NA,
                              "pval"=NA_real_,  "pval_corr"=NA_real_, "sig"=NA, 
                              "boot_freq"=NA_real_, "boot_mean"=NA_real_, "boot_stdev"=NA_real_, 
                              "boot_min"=NA_real_,"boot_max"=NA_real_, "boot_na"=NA_real_,
                              "totalA"=NA_real_, "totalB"=NA_real_)  # sum of all transcripts for that gene
  } else {
    Parameters <- list("num_replic_A"=NA_integer_, "num_replic_B"=NA_integer_)
    Genes <- data.table("parent_id"=levels(as.factor(annot$parent_id)), "DTU"=NA, 
                        "elig_transc"=NA_integer_, "elig"=NA, "elig_fx"=NA,
                        "pvalAB"=NA_real_, "pvalBA"=NA_real_,
                        "pvalAB_corr"=NA_real_, "pvalBA_corr"=NA_real_, "sig"=NA)
    Transcripts <- data.table("target_id"=annot$target_id, "parent_id"=annot$parent_id, "DTU"=NA, 
                              "sumA"=NA_real_, "sumB"=NA_real_, "elig_xp"=NA, "elig"=NA,      # sum across replicates of means across bootstraps
                              "propA"=NA_real_, "propB"=NA_real_, "Dprop"=NA_real_, "elig_fx"=NA,
                              "pval"=NA_real_,  "pval_corr"=NA_real_, "sig"=NA,
                              "totalA"=NA_real_, "totalB"=NA_real_)  # sum of all transcripts for that gene
  }
  with(Genes,       setkey(Genes, parent_id))
  with(Transcripts, setkey(Transcripts, parent_id, target_id))
  
  return(list("Parameters"=Parameters, "Genes"=Genes, "Transcripts"=Transcripts))
}


#================================================================================
#' Set up and execute the tests.
#'
#' @param counts_A A data.table of counts for condition A. x: sample, y: transcript.
#' @param counts_B A data.table of counts for condition B. x: sample, y: transcript.
#' @param tx_filter A data.table with target_id and parent_id.
#' @param testmode One of "G-test", "prop-test", "both".
#' @param full Either "full" (for complete output structure) or "short" (for bootstrapping).
#' @param count_thresh Minimum number of counts per replicate.
#' @param p_thresh The p-value threshold.
#' @param dprop_thresh Minimum difference in proportions.
#' @param correction Multiple testing correction type.
#' @return list
#' 
#' @import data.table
#' 
calculate_DTU <- function(counts_A, counts_B, tx_filter, testmode, full, count_thresh, p_thresh, dprop_thresh, correction) {
  
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
    Transcripts[(is.na(elig_xp)), elig_xp := FALSE]
    Transcripts[, elig := (elig_xp & totalA != 0 & totalB != 0 & (sumA != totalA | sumB != totalB))]  # If the entire gene is shut off, changes in proportion cannot be defined.
                                                                                                      # If sum and total are equal in both conditions, it has no detected siblings and thus cannot change in proportion.
    
    Genes[, elig_transc := Transcripts[, .(parent_id, ifelse(elig, 1, 0))][, as.integer(sum(V2)), by = parent_id][, V1] ]  # Sum of 1's is integer. Otherwise sum() changes the column to double, defeating the pre-allocation.
    Genes[, elig := elig_transc >= 2]
    
    # Biologically significant.
    Transcripts[, elig_fx := abs(Dprop) >= dprop_thresh]
    Genes[, elig_fx := Transcripts[, .(parent_id, any(elig_fx)), by = parent_id][, V2] ]

    #---------- TESTS
  
    # Proportion test.
    if ( any(testmode == c("prop-test", "both"))) {
      Transcripts[(elig), pval := as.vector(apply(Transcripts[(elig), .(sumA, sumB, totalA, totalB)], MARGIN = 1, 
                                                     FUN = function(row) { prop.test(x = row[c("sumA", "sumB")], 
                                                                                     n = row[c("totalA", "totalB")], 
                                                                                     correct = TRUE)[["p.value"]] } )) ]
      Transcripts[(elig), pval_corr := p.adjust(pval, method=correction)]
      Transcripts[(elig), sig := pval_corr < p_thresh]
      Transcripts[(elig), DTU := sig & elig_fx]
    }
    
    # G test.
    if ( any(testmode == c("G-test", "g-test", "both"))) {
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

