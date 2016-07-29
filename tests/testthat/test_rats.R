#==============================================================================
#==============================================================================
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
  expect_equal(dim(full$Transcripts)[2], 27)
  expect_equal(dim(short$Transcripts)[2], 16)
  expect_named(full$Transcripts, c("target_id", "parent_id", "DTU", "gene_DTU", "meanA", "meanB", "stdevA", "stdevB",
                                   "sumA", "sumB", "totalA", "totalB", "elig_xp", "elig", "propA", "propB", "Dprop", 
                                   "elig_fx", "pval", "pval_corr", "sig", "boot_freq", "boot_mean", "boot_stdev", 
                                   "boot_min", "boot_max", "boot_na"))
  expect_true(all(names(short$Transcripts) %in% names(full$Transcripts)))
})


#==============================================================================
#==============================================================================
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

#==============================================================================
test_that("Filters work correctly", {
  
  # !!! This test is tightly dependent on the data used for the test, in order to
  # !!! ensure correct response to specific scenarios.
  
  sim <- sim_sleuth_data(cnames=c("ONE","TWO"))
  mydtu <- call_DTU(sim$slo, sim$annot, "ONE", "TWO", verbose = FALSE)
  
  expect_equivalent(as.list(mydtu$Genes["1A1N", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(1, 1, 0, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Genes["1B1C", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(2, 1, 0, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Genes["1D1C", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(2, 1, 0, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Genes["ALLA", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(1, 1, 0, FALSE, NA))
  expect_equivalent(as.list(mydtu$Genes["ALLB", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(2, 2, 0, FALSE, NA))
  expect_equivalent(as.list(mydtu$Genes["CC", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(2, 2, 2, TRUE, TRUE))
  expect_equivalent(as.list(mydtu$Genes["LC", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(2, 2, 1, FALSE, TRUE))
  expect_equivalent(as.list(mydtu$Genes["MIX6", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(6, 5, 5, TRUE, TRUE))
  expect_equivalent(as.list(mydtu$Genes["NIB", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(1, 0, 0, FALSE, NA))
  expect_equivalent(as.list(mydtu$Genes["NN", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(2, 2, 2, TRUE, FALSE))
  
  setkey(mydtu$Transcripts, target_id)
  expect_equivalent(as.list(mydtu$Transcripts["1A1N-2", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Transcripts["1B1C.1", list(elig_xp, elig, elig_fx)]),
                    list(FALSE, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Transcripts["1B1C.2", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Transcripts["1D1C:one", list(elig_xp, elig, elig_fx)]),
                    list(FALSE, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Transcripts["1D1C:two", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Transcripts["ALLA1", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, FALSE, NA))
  expect_equivalent(as.list(mydtu$Transcripts["ALLB1", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, FALSE, NA))
  expect_equivalent(as.list(mydtu$Transcripts["ALLB2", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, FALSE, NA))
  expect_equivalent(as.list(mydtu$Transcripts["CC_a", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, TRUE))
  expect_equivalent(as.list(mydtu$Transcripts["CC_b", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, TRUE))
  expect_equivalent(as.list(mydtu$Transcripts["LC1", list(elig_xp, elig, elig_fx)]),
                    list(FALSE, FALSE, TRUE))
  expect_equivalent(as.list(mydtu$Transcripts["LC2", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, TRUE))
  expect_equivalent(as.list(mydtu$Transcripts["MIX6.c1", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, TRUE))
  expect_equivalent(as.list(mydtu$Transcripts["MIX6.c2", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, TRUE))
  expect_equivalent(as.list(mydtu$Transcripts["MIX6.c3", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, FALSE))
  expect_equivalent(as.list(mydtu$Transcripts["MIX6.c4", .(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, TRUE))
  expect_equivalent(as.list(mydtu$Transcripts["MIX6.d", .(elig_xp, elig, elig_fx)]),
                    list(FALSE, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Transcripts["MIX6.nc", .(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, FALSE))
  expect_equivalent(as.list(mydtu$Transcripts["NIB.1", .(elig_xp, elig, elig_fx)]),
                    list(FALSE, FALSE, NA))
  expect_equivalent(as.list(mydtu$Transcripts["1NN", .(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, FALSE))
  expect_equivalent(as.list(mydtu$Transcripts["2NN", .(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, FALSE))
})


#==============================================================================
#==============================================================================
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
                              testmode= "transc", correction= "bonferroni", verbose= FALSE, boots= "genes",
                              bootnum= 2, COUNTS_COL= "counts", TARGET_COL= "target", 
                              PARENT_COL= "parent", BS_TARGET_COL= "id"))
  # No false alarms with defaults.
  sim <- sim_sleuth_data(cnames=c(name_A, name_B))
  expect_silent(call_DTU(sim$slo, sim$annot, name_A, name_B, verbose = FALSE))
  
  # Annottaion is not a dataframe.
  expect_error(call_DTU(sim$slo, c("not", "a", "dataframe"), name_A, name_B, verbose = FALSE), "annot is not a data.frame.")
  # Annotation field names.
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, TARGET_COL= wrong_name, verbose = FALSE),
               "target and/or parent IDs field-names do not exist in annot", fixed= TRUE)
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, PARENT_COL= wrong_name, verbose = FALSE),
               "target and/or parent IDs field-names do not exist in annot", fixed= TRUE)
  
  # Bootstrap field names.
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, BS_TARGET_COL= wrong_name, verbose = FALSE),
               "target IDs field-name does not exist in the bootstraps", fixed= TRUE)
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, COUNTS_COL= wrong_name, verbose = FALSE),
               "counts field-name does not exist", fixed= TRUE)
  
  # Correction method.
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, correction= wrong_name, verbose = FALSE),
               "Invalid p-value correction method name", fixed= TRUE)
  
  # Covariate name.
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, varname= wrong_name, verbose = FALSE),
               "covariate name does not exist", fixed= TRUE)
  
  # Condition names.
  expect_error(call_DTU(sim$slo, sim$annot, wrong_name, name_B, verbose = FALSE),
               "conditions do not exist", fixed= TRUE)
  expect_error(call_DTU(sim$slo, sim$annot, name_A, wrong_name, verbose = FALSE),
               "conditions do not exist", fixed= TRUE)
  
  # Verbose is bool.
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, verbose="yes"),
               "not interpretable as logical", fixed= TRUE)
  
  # Probability threshold.
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, p_thresh = 666, verbose = FALSE),
               "Invalid p-value threshold", fixed= TRUE)
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, p_thresh = -0.05, verbose = FALSE),
               "Invalid p-value threshold", fixed= TRUE)
  
  # Read counts threshold.
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, count_thresh = -5, verbose = FALSE),
               "Invalid read-count threshold", fixed= TRUE)
  
  # Tests.
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, testmode="GCSE", verbose = FALSE),
               "Unrecognized value for testmode", fixed= TRUE)
  expect_silent(call_DTU(sim$slo, sim$annot, name_A, name_B, testmode="genes", verbose = FALSE))
  expect_silent(call_DTU(sim$slo, sim$annot, name_A, name_B, testmode="transc", verbose = FALSE))
  
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, boots="GCSE", verbose = FALSE),
               "Unrecognized value for boots", fixed= TRUE)
  expect_silent(call_DTU(sim$slo, sim$annot, name_A, name_B, boots="genes", bootnum = 2, verbose = FALSE))
  expect_silent(call_DTU(sim$slo, sim$annot, name_A, name_B, boots="transc", bootnum = 2, verbose = FALSE))
  
  # Number of bootstraps.
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, bootnum = -5, verbose = FALSE),
               "Invalid number of bootstraps", fixed= TRUE)
  
  # Proportion change threshold.
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, dprop_thresh = -2, verbose = FALSE),
               "Invalid proportion difference threshold", fixed= TRUE)
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, dprop_thresh = 2, verbose = FALSE),
               "Invalid proportion difference threshold", fixed= TRUE)
  
  # Inconsistent annotation.
  sim <- sim_sleuth_data(errannot_inconsistent= TRUE, cnames= c(name_A, name_B))
  expect_error(call_DTU(sim$slo, sim$annot, name_A, name_B, verbose = FALSE),
               "Inconsistent set of transcript IDs", fixed= TRUE)
})


#==============================================================================
#==============================================================================
context("DTU Output")

#==============================================================================
test_that("The output structure is correct", {
  sim <- sim_sleuth_data(cnames=c("ONE","TWO"))
  mydtu <- call_DTU(sim$slo, sim$annot, "ONE", "TWO", boots="both", bootnum=2, verbose = FALSE)
  
  expect_type(mydtu, "list")
  expect_equal(length(mydtu), 3)
  expect_named(mydtu, c("Parameters", "Genes", "Transcripts"))
  
  expect_type(mydtu$Parameters, "list")
  expect_length(mydtu$Parameters, 11)
  expect_named(mydtu$Parameters, c("var_name", "cond_A", "cond_B", "num_replic_A", "num_replic_B", "p_thresh", 
                                  "count_thresh", "dprop_thresh", "tests", "bootstrap", "bootnum"))
  
  expect_true(is.data.frame(mydtu$Genes))
  expect_equal(dim(mydtu$Genes)[2], 23)
  expect_named(mydtu$Genes, c("parent_id", "DTU", "transc_DTU", "known_transc", "detect_transc", "elig_transc",  
                             "elig", "elig_fx", "pvalAB", "pvalBA", "pvalAB_corr", "pvalBA_corr", "sig", "boot_freq", 
                             "boot_meanAB", "boot_meanBA", "boot_stdevAB", "boot_stdevBA", "boot_minAB", "boot_minBA", 
                             "boot_maxAB", "boot_maxBA", "boot_na"))
  expect_true(is.numeric(mydtu$Genes[["known_transc"]]))
  expect_true(is.numeric(mydtu$Genes[["detect_transc"]]))
  expect_true(is.numeric(mydtu$Genes[["pvalAB"]]))
  expect_true(is.numeric(mydtu$Genes[["pvalBA"]]))
  expect_true(is.numeric(mydtu$Genes[["pvalAB_corr"]]))
  expect_true(is.numeric(mydtu$Genes[["pvalBA_corr"]]))
  expect_true(is.numeric(mydtu$Genes[["boot_freq"]]))
  expect_true(is.numeric(mydtu$Genes[["boot_meanAB"]]))
  expect_true(is.numeric(mydtu$Genes[["boot_meanBA"]]))
  expect_true(is.numeric(mydtu$Genes[["boot_stdevAB"]]))
  expect_true(is.numeric(mydtu$Genes[["boot_stdevBA"]]))
  expect_true(is.numeric(mydtu$Genes[["boot_minAB"]]))
  expect_true(is.numeric(mydtu$Genes[["boot_minBA"]]))
  expect_true(is.numeric(mydtu$Genes[["boot_maxAB"]]))
  expect_true(is.numeric(mydtu$Genes[["boot_maxBA"]]))
  expect_true(is.numeric(mydtu$Genes[["boot_na"]]))
  expect_true(is.logical(mydtu$Genes[["elig"]]))
  expect_true(is.logical(mydtu$Genes[["elig_fx"]]))
  expect_true(is.logical(mydtu$Genes[["sig"]]))
  expect_true(is.logical(mydtu$Genes[["DTU"]]))
  expect_true(is.logical(mydtu$Genes[["transc_DTU"]]))
  
  expect_true(is.data.frame(mydtu$Transcripts))
  expect_equal(dim(mydtu$Transcripts)[2], 27)
  expect_named(mydtu$Transcripts, c("target_id", "parent_id", "DTU", "gene_DTU", "meanA", "meanB", "stdevA", "stdevB",
                                   "sumA", "sumB", "totalA", "totalB", "elig_xp", "elig", "propA", "propB", "Dprop", 
                                   "elig_fx", "pval", "pval_corr", "sig", "boot_freq", "boot_mean", "boot_stdev", 
                                   "boot_min", "boot_max", "boot_na"))
  expect_true(is.logical(mydtu$Transcripts[["elig_xp"]]))
  expect_true(is.logical(mydtu$Transcripts[["elig"]]))
  expect_true(is.logical(mydtu$Transcripts[["elig_fx"]]))
  expect_true(is.logical(mydtu$Transcripts[["sig"]]))
  expect_true(is.logical(mydtu$Transcripts[["DTU"]]))
  expect_true(is.logical(mydtu$Transcripts[["gene_DTU"]]))
  expect_true(is.numeric(mydtu$Transcripts[["propA"]]))
  expect_true(is.numeric(mydtu$Transcripts[["propB"]]))
  expect_true(is.numeric(mydtu$Transcripts[["Dprop"]]))
  expect_true(is.numeric(mydtu$Transcripts[["sumA"]]))
  expect_true(is.numeric(mydtu$Transcripts[["sumB"]]))
  expect_true(is.numeric(mydtu$Transcripts[["meanA"]]))
  expect_true(is.numeric(mydtu$Transcripts[["meanB"]]))
  expect_true(is.numeric(mydtu$Transcripts[["stdevA"]]))
  expect_true(is.numeric(mydtu$Transcripts[["stdevB"]]))
  expect_true(is.numeric(mydtu$Transcripts[["pval"]]))
  expect_true(is.numeric(mydtu$Transcripts[["pval_corr"]]))
  expect_true(is.numeric(mydtu$Transcripts[["boot_freq"]]))
  expect_true(is.numeric(mydtu$Transcripts[["boot_mean"]]))
  expect_true(is.numeric(mydtu$Transcripts[["boot_stdev"]]))
  expect_true(is.numeric(mydtu$Transcripts[["boot_min"]]))
  expect_true(is.numeric(mydtu$Transcripts[["boot_max"]]))
  expect_true(is.numeric(mydtu$Transcripts[["boot_na"]]))

  mydtu <- call_DTU(sim$slo, sim$annot, "ONE", "TWO", verbose = FALSE)  
  expect_false(any(c("boot_freq", "boot_mean", "boot_stdev", "boot_min", "boot_max", "boot_na") %in% names(mydtu$Genes)))
  expect_false(any(c("boot_freq", "boot_mean", "boot_stdev", "boot_min", "boot_max", "boot_na") %in% names(mydtu$Transcripts)))
  
  mydtu <- call_DTU(sim$slo, sim$annot, "ONE", "TWO", boots="transc", bootnum=2, verbose = FALSE)  
  expect_false( any(c("boot_freq", "boot_mean", "boot_stdev", "boot_min", "boot_max", "boot_na") %in% names(mydtu$Genes)))
  
  mydtu <- call_DTU(sim$slo, sim$annot, "ONE", "TWO", boots="genes", bootnum=2, verbose = FALSE)  
  expect_false(any(c("boot_freq", "boot_mean", "boot_stdev", "boot_min", "boot_max", "boot_na") %in% names(mydtu$Transcripts)))
})

#==============================================================================
test_that("The output content is complete", {
  sim <- sim_sleuth_data(cnames=c("ONE","TWO"))
  mydtu <- call_DTU(sim$slo, sim$annot, "ONE", "TWO", boots="both", bootnum=2, verbose = FALSE)
  
  expect_false(is.na(mydtu$Parameters$"var_name"))
  expect_false(is.na(mydtu$Parameters$"cond_A"))
  expect_false(is.na(mydtu$Parameters$"cond_B"))
  expect_false(is.na(mydtu$Parameters$"num_replic_A"))
  expect_false(is.na(mydtu$Parameters$"num_replic_B"))
  expect_false(is.na(mydtu$Parameters$"p_thresh"))
  expect_false(is.na(mydtu$Parameters$"dprop_thresh"))
  expect_false(is.na(mydtu$Parameters$"count_thresh"))
  expect_false(is.na(mydtu$Parameters$"tests"))
  expect_false(is.na(mydtu$Parameters$"bootstrap"))
  expect_false(is.na(mydtu$Parameters$"bootnum"))

  expect_false(all(is.na(mydtu$Genes[["parent_id"]])))
  expect_false(all(is.na(mydtu$Genes[["known_transc"]])))
  expect_false(all(is.na(mydtu$Genes[["detect_transc"]])))
  expect_false(all(is.na(mydtu$Genes[["pvalAB"]])))
  expect_false(all(is.na(mydtu$Genes[["pvalBA"]])))
  expect_false(all(is.na(mydtu$Genes[["pvalAB_corr"]])))
  expect_false(all(is.na(mydtu$Genes[["pvalBA_corr"]])))
  expect_false(all(is.na(mydtu$Genes[["boot_freq"]])))
  expect_false(all(is.na(mydtu$Genes[["boot_meanAB"]])))
  expect_false(all(is.na(mydtu$Genes[["boot_meanBA"]])))
  expect_false(all(is.na(mydtu$Genes[["boot_stdevAB"]])))
  expect_false(all(is.na(mydtu$Genes[["boot_stdevBA"]])))
  expect_false(all(is.na(mydtu$Genes[["boot_minAB"]])))
  expect_false(all(is.na(mydtu$Genes[["boot_minBA"]])))
  expect_false(all(is.na(mydtu$Genes[["boot_maxAB"]])))
  expect_false(all(is.na(mydtu$Genes[["boot_maxBA"]])))
  expect_false(all(is.na(mydtu$Genes[["boot_na"]])))
  expect_false(all(is.na(mydtu$Genes[["elig"]])))
  expect_false(all(is.na(mydtu$Genes[["elig_fx"]])))
  expect_false(all(is.na(mydtu$Genes[["sig"]])))
  expect_false(all(is.na(mydtu$Genes[["DTU"]])))
  expect_false(all(is.na(mydtu$Genes[["transc_DTU"]])))
  
  expect_false(all(is.na(mydtu$Transcripts[["target_id"]])))
  expect_false(all(is.na(mydtu$Transcripts[["parent_id"]])))
  expect_false(all(is.na(mydtu$Transcripts[["DTU"]])))
  expect_false(all(is.na(mydtu$Transcripts[["gene_DTU"]])))
  expect_false(all(is.na(mydtu$Transcripts[["meanA"]])))
  expect_false(all(is.na(mydtu$Transcripts[["meanB"]])))
  expect_false(all(is.na(mydtu$Transcripts[["stdevA"]])))
  expect_false(all(is.na(mydtu$Transcripts[["stdevB"]])))
  expect_false(all(is.na(mydtu$Transcripts[["sumA"]])))
  expect_false(all(is.na(mydtu$Transcripts[["sumB"]])))
  expect_false(all(is.na(mydtu$Transcripts[["elig_xp"]])))
  expect_false(all(is.na(mydtu$Transcripts[["elig"]])))
  expect_false(all(is.na(mydtu$Transcripts[["propA"]])))
  expect_false(all(is.na(mydtu$Transcripts[["propB"]])))
  expect_false(all(is.na(mydtu$Transcripts[["Dprop"]])))
  expect_false(all(is.na(mydtu$Transcripts[["elig_fx"]])))
  expect_false(all(is.na(mydtu$Transcripts[["pval"]])))
  expect_false(all(is.na(mydtu$Transcripts[["pval_corr"]])))
  expect_false(all(is.na(mydtu$Transcripts[["sig"]])))
  expect_false(all(is.na(mydtu$Transcripts[["boot_freq"]])))
  expect_false(all(is.na(mydtu$Transcripts[["boot_mean"]])))
  expect_false(all(is.na(mydtu$Transcripts[["boot_stdev"]])))
  expect_false(all(is.na(mydtu$Transcripts[["boot_min"]])))
  expect_false(all(is.na(mydtu$Transcripts[["boot_max"]])))
  expect_false(all(is.na(mydtu$Transcripts[["boot_na"]])))
})


#==============================================================================
test_that("The summaries work", {
  sim <- sim_sleuth_data(cnames=c("ONE","TWO"))
  mydtu <- call_DTU(sim$slo, sim$annot, "ONE", "TWO", boots="both", bootnum=2, verbose = FALSE)
  
  expect_silent(tally <- dtu_summary(mydtu))
  expect_true(is.numeric(tally))
  expect_named(tally, c("DTU genes", "non-DTU genes", "NA genes", "DTU transcripts", "non-DTU transcripts", "NA transcripts"))
  expect_false(any(is.na(tally)))
  
  ids <- get_dtu_ids(mydtu)
  expect_type(ids, "list")
  expect_named(ids, c("dtu-genes", "dtu-transc", "ndtu-genes", "ndtu-transc", "na-genes", "na-transc"))
  for (v in ids) {
    expect_false(any(is.na(v)))
  }
})
