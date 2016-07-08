context("DTU internal structures")

#==============================================================================
test_that("The reporting structures are not created correctly", {
  full <- alloc_out(mini_anno, "full")
  short <- alloc_out(mini_anno, "short")
  
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
                                  "count_thresh", "tests", "bootstrap", "bootnum", "threads"))
  expect_true(all(names(short$Parameters) %in% names(full$Parameters)))
  
  expect_true(is.data.frame(full$Genes))
  expect_true(is.data.frame(short$Genes))
  expect_equal(dim(full$Genes)[2], 23)
  expect_equal(dim(short$Genes)[2], 7)
  expect_named(full$Genes, c("parent_id", "known_transc", "detect_transc", "eligible", "Pt_DTU", "Gt_DTU", 
                             "Gt_dtuAB", "Gt_dtuBA", "Gt_pvalAB", "Gt_pvalBA", "Gt_pvalAB_corr", "Gt_pvalBA_corr", 
                             "Gt_boot_dtuAB", "Gt_boot_dtuBA", "Gt_boot_meanAB", "Gt_boot_meanBA", "Gt_boot_stdevAB", 
                             "Gt_boot_stdevBA", "Gt_boot_minAB", "Gt_boot_minBA", "Gt_boot_maxAB", "Gt_boot_maxBA", 
                             "Gt_boot_na"))
  expect_true(all(names(short$Genes) %in% names(full$Genes)))
  
  expect_true(is.data.frame(full$Transcripts))
  expect_true(is.data.frame(short$Transcripts))
  expect_equal(dim(full$Transcripts)[2], 24)
  expect_equal(dim(short$Transcripts)[2], 12)
  expect_named(full$Transcripts, c("target_id", "parent_id", "propA", "propB", "Dprop", "eligible", "Gt_DTU", "Pt_DTU", 
                                   "Pt_pval", "Pt_pval_corr", "Pt_boot_dtu", "Pt_boot_mean", "Pt_boot_stdev", 
                                   "Pt_boot_min", "Pt_boot_max", "Pt_boot_na",
                                   "sumA", "sumB", "meanA", "meanB", "stdevA", "stdevB", "totalA", "totalB"))
  expect_true(all(names(short$Transcripts) %in% names(full$Transcripts)))
})


context("DTU internal data munging")

#==============================================================================
test_that("Samples are not grouped correctly", {
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
test_that("Bootstrapped counts are not extracted correctly", {
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
test_that("The input checks don't work", {
  name_A <- "one"
  name_B <- "two"
  wrong_name <- "RUBBISH_COLUMN_NAME"
  
  # No false alarms.
  sim <- sim_sleuth_data(varname="waffles", COUNTS_COL="counts", TARGET_COL="target" , PARENT_COL="parent", BS_TARGET_COL="id", cnames=c("AAAA","BBBB"))
  expect_silent(calculate_DTU(sim$slo, sim$annot, "AAAA", "BBBB", varname = "waffles", p_thresh = 0.01, count_thresh = 10,
                              testmode = "prop-test", correction = "bonferroni", verbose = FALSE, boots = "g-test",
                              bootnum = 2, threads = 1, COUNTS_COL = "counts", TARGET_COL = "target", 
                              PARENT_COL = "parent", BS_TARGET_COL = "id"))
  sim <- sim_sleuth_data(cnames=c(name_A, name_B))
  expect_silent(calculate_DTU(sim$slo, sim$annot, name_A, name_B))
  
  # Annottaion is not a dataframe.
  expect_error(calculate_DTU(sim$slo, c("not", "a", "dataframe"), name_A, name_B), "annot is not a data.frame.")
  # Annotation field names.
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, TARGET_COL= wrong_name),
               "target and/or parent IDs field-names do not exist in annot", fixed= TRUE)
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, PARENT_COL= wrong_name),
               "target and/or parent IDs field-names do not exist in annot", fixed= TRUE)
  
  # Bootstrap field names.
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, BS_TARGET_COL= wrong_name),
               "target IDs field-name does not exist in the bootstraps", fixed= TRUE)
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, COUNTS_COL= wrong_name),
               "counts field-name does not exist", fixed= TRUE)
  
  # Correction method.
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, correction= wrong_name),
               "Invalid p-value correction method name", fixed= TRUE)
  
  # Covariate name.
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, varname= wrong_name),
               "covariate name does not exist", fixed= TRUE)
  
  # Condition names.
  expect_error(calculate_DTU(sim$slo, sim$annot, wrong_name, name_B),
               "conditions do not exist", fixed= TRUE)
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, wrong_name),
               "conditions do not exist", fixed= TRUE)
  
  # Verbose is bool.
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, verbose="yes"),
               "verbose must be a logical", fixed= TRUE)
  
  # Probability threshold.
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, p_thresh = 666),
               "Invalid p-value threshold", fixed= TRUE)
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, p_thresh = -0.05),
               "Invalid p-value threshold", fixed= TRUE)
  
  # Read counts threshold.
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, count_thresh = -5),
               "Invalid read-count threshold", fixed= TRUE)
  
  # Tests.
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, testmode="GCSE"),
               "Unrecognized value for testmode", fixed= TRUE)
  expect_silent(calculate_DTU(sim$slo, sim$annot, name_A, name_B, testmode="g-test"))
  expect_silent(calculate_DTU(sim$slo, sim$annot, name_A, name_B, testmode="prop-test"))
  
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, boots="GCSE"),
               "Unrecognized value for boots", fixed= TRUE)
  expect_silent(calculate_DTU(sim$slo, sim$annot, name_A, name_B, boots="g-test", bootnum = 2))
  expect_silent(calculate_DTU(sim$slo, sim$annot, name_A, name_B, boots="prop-test", bootnum = 2))
  
  # Number of bootstraps.
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, bootnum = -5),
               "Invalid number of bootstraps", fixed= TRUE)
  
  # Inconsistent annotation.
  sim <- sim_sleuth_data(errannot_inconsistent = TRUE)
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B),
               "Inconsistent set of transcript IDs", fixed= TRUE)
})


context("DTU Output")

#==============================================================================
test_that("The output structure is not correct", {
  sim <- sim_sleuth_data(cnames=c("ONE","TWO"))
  full <- calculate_DTU(sim$slo, sim$annot, "ONE", "TWO", boots="both", bootnum=2)
  
  expect_type(full, "list")
  expect_equal(length(full), 3)
  expect_named(full, c("Parameters", "Genes", "Transcripts"))
  
  expect_type(full$Parameters, "list")
  expect_length(full$Parameters, 11)
  expect_named(full$Parameters, c("var_name", "cond_A", "cond_B", "num_replic_A", "num_replic_B", "p_thresh", 
                                  "count_thresh", "tests", "bootstrap", "bootnum", "threads"))
  
  expect_true(is.data.frame(full$Genes))
  expect_equal(dim(full$Genes)[2], 23)
  expect_named(full$Genes, c("parent_id", "known_transc", "detect_transc", "eligible", "Pt_DTU", "Gt_DTU", 
                             "Gt_dtuAB", "Gt_dtuBA", "Gt_pvalAB", "Gt_pvalBA", "Gt_pvalAB_corr", "Gt_pvalBA_corr", 
                             "Gt_boot_dtuAB", "Gt_boot_dtuBA", "Gt_boot_meanAB", "Gt_boot_meanBA", "Gt_boot_stdevAB", 
                             "Gt_boot_stdevBA", "Gt_boot_minAB", "Gt_boot_minBA", "Gt_boot_maxAB", "Gt_boot_maxBA", 
                             "Gt_boot_na"))
  expect_true(is.numeric(full$Genes[["known_transc"]]))
  expect_true(is.numeric(full$Genes[["detect_transc"]]))
  expect_true(is.numeric(full$Genes[["Gt_pvalAB"]]))
  expect_true(is.numeric(full$Genes[["Gt_pvalBA"]]))
  expect_true(is.numeric(full$Genes[["Gt_pvalAB_corr"]]))
  expect_true(is.numeric(full$Genes[["Gt_pvalBA_corr"]]))
  expect_true(is.numeric(full$Genes[["Gt_boot_dtuAB"]]))
  expect_true(is.numeric(full$Genes[["Gt_boot_dtuBA"]]))
  expect_true(is.numeric(full$Genes[["Gt_boot_meanAB"]]))
  expect_true(is.numeric(full$Genes[["Gt_boot_meanBA"]]))
  expect_true(is.numeric(full$Genes[["Gt_boot_stdevAB"]]))
  expect_true(is.numeric(full$Genes[["Gt_boot_stdevBA"]]))
  expect_true(is.numeric(full$Genes[["Gt_boot_minAB"]]))
  expect_true(is.numeric(full$Genes[["Gt_boot_minBA"]]))
  expect_true(is.numeric(full$Genes[["Gt_boot_maxAB"]]))
  expect_true(is.numeric(full$Genes[["Gt_boot_maxBA"]]))
  expect_true(is.numeric(full$Genes[["Gt_boot_na"]]))
  expect_true(is.logical(full$Genes[["eligible"]]))
  expect_true(is.logical(full$Genes[["Pt_DTU"]]))
  expect_true(is.logical(full$Genes[["Gt_DTU"]]))
  expect_true(is.logical(full$Genes[["Gt_dtuAB"]]))
  expect_true(is.logical(full$Genes[["Gt_dtuBA"]]))
  
  expect_true(is.data.frame(full$Transcripts))
  expect_equal(dim(full$Transcripts)[2], 24)
  expect_named(full$Transcripts, c("target_id", "parent_id", "propA", "propB", "Dprop", "eligible", "Gt_DTU", "Pt_DTU", 
                                   "Pt_pval", "Pt_pval_corr", "Pt_boot_dtu", "Pt_boot_mean", "Pt_boot_stdev", 
                                   "Pt_boot_min", "Pt_boot_max", "Pt_boot_na",
                                   "sumA", "sumB", "meanA", "meanB", "stdevA", "stdevB", "totalA", "totalB"))
  expect_true(is.logical(full$Transcripts[["eligible"]]))
  expect_true(is.logical(full$Transcripts[["Gt_DTU"]]))
  expect_true(is.logical(full$Transcripts[["Pt_DTU"]]))
  expect_true(is.numeric(full$Transcripts[["propA"]]))
  expect_true(is.numeric(full$Transcripts[["propB"]]))
  expect_true(is.numeric(full$Transcripts[["Dprop"]]))
  expect_true(is.numeric(full$Transcripts[["sumA"]]))
  expect_true(is.numeric(full$Transcripts[["sumB"]]))
  expect_true(is.numeric(full$Transcripts[["meanA"]]))
  expect_true(is.numeric(full$Transcripts[["meanB"]]))
  expect_true(is.numeric(full$Transcripts[["stdevA"]]))
  expect_true(is.numeric(full$Transcripts[["stdevB"]]))
  expect_true(is.numeric(full$Transcripts[["totalA"]]))
  expect_true(is.numeric(full$Transcripts[["totalB"]]))
  expect_true(is.numeric(full$Transcripts[["Pt_pval"]]))
  expect_true(is.numeric(full$Transcripts[["Pt_pval_corr"]]))
  expect_true(is.numeric(full$Transcripts[["Pt_boot_dtu"]]))
  expect_true(is.numeric(full$Transcripts[["Pt_boot_mean"]]))
  expect_true(is.numeric(full$Transcripts[["Pt_boot_stdev"]]))
  expect_true(is.numeric(full$Transcripts[["Pt_boot_min"]]))
  expect_true(is.numeric(full$Transcripts[["Pt_boot_max"]]))
  expect_true(is.numeric(full$Transcripts[["Pt_boot_na"]]))
})
