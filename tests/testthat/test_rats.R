context("DTU internal structures")

#==============================================================================
test_that("The reporting structures are created correctly", {
  sim <- sim_sleuth_data()
  full <- alloc_out(sim$annot, "full")
  short <- alloc_out(sim$annot, "short")
  
  expect_type(full, "list")
  expect_type(short, "list")
  expect_equal(length(full), 3)
  expect_equal(length(full), length(short))
  expect_named(full, c("Parameters", "Genes", "Transcripts"))
  expect_true(all(names(full) == names(short)))
  
  expect_type(full$Parameters, "list")
  expect_true(typeof(full$Parameters) == typeof(short$Parameters))
  expect_length(full$Parameters, 11)
  expect_length(short$Parameters, 2)
  expect_named(full$Parameters, c("var_name", "cond_A", "cond_B", "num_replic_A", "num_replic_B", "p_thresh", 
                                  "count_thresh", "dprop_thresh", "tests", "bootstrap", "bootnum"))
  expect_true(all(names(short$Parameters) %in% names(full$Parameters)))
  
  expect_true(is.data.frame(full$Genes))
  expect_true(is.data.frame(short$Genes))
  expect_equal(dim(full$Genes)[2], 23)
  expect_equal(dim(short$Genes)[2], 10)
  expect_named(full$Genes, c("parent_id", "DTU", "transc_DTU", "known_transc", "detect_transc", "elig_transc",  
                             "elig", "elig_fx", "pvalAB", "pvalBA", "pvalAB_corr", "pvalBA_corr", "sig", "boot_freq", 
                             "boot_meanAB", "boot_meanBA", "boot_stdevAB", "boot_stdevBA", "boot_minAB", "boot_minBA", 
                             "boot_maxAB", "boot_maxBA", "boot_na"))
  expect_true(all(names(short$Genes) %in% names(full$Genes)))
  
  expect_true(is.data.frame(full$Transcripts))
  expect_true(is.data.frame(short$Transcripts))
  expect_equal(dim(full$Transcripts)[2], 26)
  expect_equal(dim(short$Transcripts)[2], 15)
  expect_named(full$Transcripts, c("target_id", "parent_id", "DTU", "gene_DTU", "meanA", "meanB", "stdevA", "stdevB",
                                   "sumA", "sumB", "elig_xp", "propA", "propB", "Dprop", "elig_fx", "pval", "pval_corr",
                                   "sig", "boot_freq", "boot_mean", "boot_stdev", "boot_min", "boot_max", "boot_na",
                                   "totalA", "totalB"))
  expect_true(all(names(short$Transcripts) %in% names(full$Transcripts)))
})


context("DTU internal data munging")

#==============================================================================
test_that("Samples are grouped correctly", {
  sim <- sim_sleuth_data()
  r <- group_samples(sim$slo$sample_to_covariates)
  
  # number of covariates
  expect_equal(length(r), length(sim$slo$sample_to_covariates))
  # names of covariates
  expect_named(r, names(sim$slo$sample_to_covariates))
  # number of values of each covariate
  expect_equal(length(r[[1]]), length(levels(as.factor(sim$slo$sample_to_covariates[[1]]))))
  expect_equal(length(r[[2]]), length(levels(as.factor(sim$slo$sample_to_covariates[[2]]))))
  # total number of samples
  expect_equal(sum(sapply(r[[1]],length)), length(sim$slo$sample_to_covariates[[1]]))
  expect_equal(sum(sapply(r[[2]],length)), length(sim$slo$sample_to_covariates[[2]]))
})


#==============================================================================
test_that("Bootstrapped counts are extracted correctly", {
  samples <- c(1,3)
  bst <- "id"
  cnt <- "counts"
  sim <- sim_sleuth_data(COUNTS_COL = cnt, BS_TARGET_COL = bst)
  lr <- denest_boots(sim$slo, sim$annot[[1]], samples, cnt, bst)
  
  for (i in 1:length(lr)) {
    # The transcripts supposed to be there are there.
    expect_true(all(sim$isx %in% as.character(lr[[i]]$target_id)))
    # No NA.
    expect_false(any(is.na(lr[[i]])))
    
    # Number of bootstraps per sample.
    expect_equal(length(lr[[i]]) - 1, length(sim$slo$kal[[samples[i]]]$bootstrap))  # the last column in lr[] is ID
    
    # All target counts pulled from the correct bootstraps and the correct transcripts.
    counts_ok <- sapply(1:(length(lr[[i]]) - 1), function(j) {
      fltr1 <- match(sim$isx, sim$slo$kal[[samples[i]]]$bootstrap[[j]][[bst]])  # Where in the boot are the expected IDs.
      fltr2 <- match(sim$isx, lr[[i]][["target_id"]])                           # Where in the extracted counts are the expected IDs.
      all(sim$slo$kal[[samples[i]]]$bootstrap[[j]][[cnt]][fltr1] == lr[[i]][[j]][fltr2])  # Both vectors' elements are ordered by the same IDs.
    })
    expect_true(all(counts_ok))
    
    # IDs in annotation, but not in bootstraps, should be 0.
    missing_from_boots_ok <- sapply(1:(length(lr[[i]]) - 1), function(j) {
      nib <- setdiff(sim$annot$target_id, sim$isx)
      nib <- match(nib, lr[[i]][["target_id"]])
      all(lr[[i]][[j]][nib] == 0)
    })
    expect_true(all(missing_from_boots_ok))
      
    # IDs in bootstrap, but not in annotation, should be absent.
    missing_from_annot_ok <- sapply(1:(length(lr[[i]]) - 1), function(j) {
      nia <- setdiff(sim$slo$kal[[i]]$bootstrap[[j]], sim$isx)
      any(nia %in% lr[[i]][["target_id"]])
    })
    expect_false(any(missing_from_annot_ok))
  }
})


context("DTU Input checks.")

#==============================================================================
test_that("The input checks work", {
  name_A <- "one"
  name_B <- "two"
  wrong_name <- "RUBBISH_COLUMN_NAME"
  
  # No false alarms with valid parameters.
  sim <- sim_sleuth_data(varname= "waffles", COUNTS_COL= "counts", TARGET_COL= "target" , PARENT_COL= "parent", 
                         BS_TARGET_COL= "id", cnames= c("AAAA","BBBB"))
  expect_silent(call_DTU(sim$slo, sim$annot, "AAAA", "BBBB", varname= "waffles", p_thresh= 0.01, count_thresh= 10,
                              testmode= "prop-test", correction= "bonferroni", verbose= FALSE, boots= "g-test",
                              bootnum= 2, COUNTS_COL= "counts", TARGET_COL= "target", 
                              PARENT_COL= "parent", BS_TARGET_COL= "id"))
  # No false alarms with defaults.
  sim <- sim_sleuth_data(cnames=c(name_A, name_B))
  expect_silent(call_DTU(sim$slo, sim$annot, name_A, name_B))
  
  # Annottaion is not a dataframe.
  expect_error(call_DTU(sim$slo, c("not", "a", "dataframe"), name_A, name_B), "annot is not a data.frame.")
  # Annotation field names.
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, TARGET_COL= wrong_name),
               "target and/or parent IDs field-names do not exist in annot", fixed= TRUE)
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, PARENT_COL= wrong_name),
               "target and/or parent IDs field-names do not exist in annot", fixed= TRUE)
  
  # Bootstrap field names.
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, BS_TARGET_COL= wrong_name),
               "target IDs field-name does not exist in the bootstraps", fixed= TRUE)
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, COUNTS_COL= wrong_name),
               "counts field-name does not exist", fixed= TRUE)
  
  # Correction method.
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, correction= wrong_name),
               "Invalid p-value correction method name", fixed= TRUE)
  
  # Covariate name.
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, varname= wrong_name),
               "covariate name does not exist", fixed= TRUE)
  
  # Condition names.
  expect_error(call_DTU(sim$slo, sim$annot, wrong_name, name_B),
               "conditions do not exist", fixed= TRUE)
  expect_error(call_DTU(sim$slo, sim$annot, name_A, wrong_name),
               "conditions do not exist", fixed= TRUE)
  
  # Verbose is bool.
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, verbose="yes"),
               "verbose must be a logical", fixed= TRUE)
  
  # Probability threshold.
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, p_thresh = 666),
               "Invalid p-value threshold", fixed= TRUE)
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, p_thresh = -0.05),
               "Invalid p-value threshold", fixed= TRUE)
  
  # Read counts threshold.
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, count_thresh = -5),
               "Invalid read-count threshold", fixed= TRUE)
  
  # Tests.
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, testmode="GCSE"),
               "Unrecognized value for testmode", fixed= TRUE)
  expect_silent(call_DTU(sim$slo, sim$annot, name_A, name_B, testmode="g-test"))
  expect_silent(call_DTU(sim$slo, sim$annot, name_A, name_B, testmode="prop-test"))
  
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, boots="GCSE"),
               "Unrecognized value for boots", fixed= TRUE)
  expect_silent(call_DTU(sim$slo, sim$annot, name_A, name_B, boots="g-test", bootnum = 2))
  expect_silent(call_DTU(sim$slo, sim$annot, name_A, name_B, boots="prop-test", bootnum = 2))
  
  # Number of bootstraps.
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, bootnum = -5),
               "Invalid number of bootstraps", fixed= TRUE)
  
  # Proportion change threshold.
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, dprop_thresh = -2),
               "Invalid proportion difference threshold", fixed= TRUE)
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, dprop_thresh = 2),
               "Invalid proportion difference threshold", fixed= TRUE)
  
  # Inconsistent annotation.
  sim <- sim_sleuth_data(errannot_inconsistent= TRUE, cnames= c(name_A, name_B))
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B),
               "Inconsistent set of transcript IDs", fixed= TRUE)
})


context("DTU Output")

#==============================================================================
test_that("The output structure is correct", {
  sim <- sim_sleuth_data(cnames=c("ONE","TWO"))
  full <- call_DTU(sim$slo, sim$annot, "ONE", "TWO", boots="both", bootnum=2)
  
  expect_type(full, "list")
  expect_equal(length(full), 3)
  expect_named(full, c("Parameters", "Genes", "Transcripts"))
  
  expect_type(full$Parameters, "list")
  expect_length(full$Parameters, 11)
  expect_named(full$Parameters, c("var_name", "cond_A", "cond_B", "num_replic_A", "num_replic_B", "p_thresh", 
                                  "count_thresh", "dprop_thresh", "tests", "bootstrap", "bootnum"))
  
  expect_true(is.data.frame(full$Genes))
  expect_equal(dim(full$Genes)[2], 23)
  expect_named(full$Genes, c("parent_id", "DTU", "transc_DTU", "known_transc", "detect_transc", "elig_transc",  
                             "elig", "elig_fx", "pvalAB", "pvalBA", "pvalAB_corr", "pvalBA_corr", "sig", "boot_freq", 
                             "boot_meanAB", "boot_meanBA", "boot_stdevAB", "boot_stdevBA", "boot_minAB", "boot_minBA", 
                             "boot_maxAB", "boot_maxBA", "boot_na"))
  expect_true(is.numeric(full$Genes[["known_transc"]]))
  expect_true(is.numeric(full$Genes[["detect_transc"]]))
  expect_true(is.numeric(full$Genes[["pvalAB"]]))
  expect_true(is.numeric(full$Genes[["pvalBA"]]))
  expect_true(is.numeric(full$Genes[["pvalAB_corr"]]))
  expect_true(is.numeric(full$Genes[["pvalBA_corr"]]))
  expect_true(is.numeric(full$Genes[["boot_freq"]]))
  expect_true(is.numeric(full$Genes[["boot_meanAB"]]))
  expect_true(is.numeric(full$Genes[["boot_meanBA"]]))
  expect_true(is.numeric(full$Genes[["boot_stdevAB"]]))
  expect_true(is.numeric(full$Genes[["boot_stdevBA"]]))
  expect_true(is.numeric(full$Genes[["boot_minAB"]]))
  expect_true(is.numeric(full$Genes[["boot_minBA"]]))
  expect_true(is.numeric(full$Genes[["boot_maxAB"]]))
  expect_true(is.numeric(full$Genes[["boot_maxBA"]]))
  expect_true(is.numeric(full$Genes[["boot_na"]]))
  expect_true(is.logical(full$Genes[["elig"]]))
  expect_true(is.logical(full$Genes[["elig_fx"]]))
  expect_true(is.logical(full$Genes[["sig"]]))
  expect_true(is.logical(full$Genes[["DTU"]]))
  expect_true(is.logical(full$Genes[["transc_DTU"]]))
  
  expect_true(is.data.frame(full$Transcripts))
  expect_equal(dim(full$Transcripts)[2], 24)
  expect_named(full$Transcripts, c("target_id", "parent_id", "DTU", "gene_DTU", "meanA", "meanB", "stdevA", "stdevB",
                                   "sumA", "sumB", "elig_xp", "propA", "propB", "Dprop", "elig_fx", "pval", "pval_corr",
                                   "sig", "boot_freq", "boot_mean", "boot_stdev", "boot_min", "boot_max", "boot_na"))
  expect_true(is.logical(full$Transcripts[["elig_xp"]]))
  expect_true(is.logical(full$Transcripts[["elig_fx"]]))
  expect_true(is.logical(full$Transcripts[["sig"]]))
  expect_true(is.logical(full$Transcripts[["DTU"]]))
  expect_true(is.logical(full$Transcripts[["gene_DTU"]]))
  expect_true(is.numeric(full$Transcripts[["propA"]]))
  expect_true(is.numeric(full$Transcripts[["propB"]]))
  expect_true(is.numeric(full$Transcripts[["Dprop"]]))
  expect_true(is.numeric(full$Transcripts[["sumA"]]))
  expect_true(is.numeric(full$Transcripts[["sumB"]]))
  expect_true(is.numeric(full$Transcripts[["meanA"]]))
  expect_true(is.numeric(full$Transcripts[["meanB"]]))
  expect_true(is.numeric(full$Transcripts[["stdevA"]]))
  expect_true(is.numeric(full$Transcripts[["stdevB"]]))
  expect_true(is.numeric(full$Transcripts[["pval"]]))
  expect_true(is.numeric(full$Transcripts[["pval_corr"]]))
  expect_true(is.numeric(full$Transcripts[["boot_freq"]]))
  expect_true(is.numeric(full$Transcripts[["boot_mean"]]))
  expect_true(is.numeric(full$Transcripts[["boot_stdev"]]))
  expect_true(is.numeric(full$Transcripts[["boot_min"]]))
  expect_true(is.numeric(full$Transcripts[["boot_max"]]))
  expect_true(is.numeric(full$Transcripts[["boot_na"]]))
})
