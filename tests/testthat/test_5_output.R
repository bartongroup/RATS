#==============================================================================
#==============================================================================
context("DTU Output")

#==============================================================================
test_that("The output structure is correct", {
  sim <- sim_boot_data(clean=TRUE)
  mydtu <- call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= "ONE", name_B= "TWO", qbootnum=2, qboot=TRUE, rboot=TRUE, verbose = FALSE, lean=FALSE)
  
  expect_type(mydtu, "list")
  expect_equal(length(mydtu), 4)
  expect_named(mydtu, c("Parameters", "Genes", "Transcripts", "Abundances"))
  
  expect_type(mydtu$Parameters, "list")
  expect_length(mydtu$Parameters, 27)
  expect_named(mydtu$Parameters, c("description", "time", "rats_version", "R_version",
                                   "var_name", "cond_A", "cond_B", "data_type", "num_replic_A", "num_replic_B", "num_genes", "num_transc",
                                   "tests", "p_thresh", "abund_thresh", "dprop_thresh", "correction", "abund_scaling",
                                   "quant_boot", "quant_reprod_thresh", "quant_bootnum",
                                   "rep_boot", "rep_reprod_thresh", "rep_bootnum", "seed", "reckless", "lean"))
  
  expect_true(is.data.frame(mydtu$Genes))
  expect_equal(dim(mydtu$Genes)[2], 24)
  expect_named(mydtu$Genes, c("parent_id", "elig", "sig", "elig_fx", "quant_reprod", "rep_reprod", "DTU", "transc_DTU",
                              "known_transc", "detect_transc", "elig_transc", "maxDprop", "pval", "pval_corr", 
                              "quant_p_median", "quant_p_min", "quant_p_max", "quant_na_freq", "quant_dtu_freq",
                              "rep_p_median", "rep_p_min", "rep_p_max", "rep_na_freq", "rep_dtu_freq") )
  expect_true(is.numeric(mydtu$Genes[["known_transc"]]))
  expect_true(is.numeric(mydtu$Genes[["detect_transc"]]))
  expect_true(is.numeric(mydtu$Genes[["maxDprop"]]))
  expect_true(is.numeric(mydtu$Genes[["pval"]]))
  expect_true(is.numeric(mydtu$Genes[["pval_corr"]]))
  expect_true(is.numeric(mydtu$Genes[["quant_dtu_freq"]]))
  expect_true(is.numeric(mydtu$Genes[["quant_p_median"]]))
  expect_true(is.numeric(mydtu$Genes[["quant_p_min"]]))
  expect_true(is.numeric(mydtu$Genes[["quant_p_max"]]))
  expect_true(is.numeric(mydtu$Genes[["quant_na_freq"]]))
  expect_true(is.logical(mydtu$Genes[["quant_reprod"]]))
  expect_true(is.numeric(mydtu$Genes[["rep_dtu_freq"]]))
  expect_true(is.numeric(mydtu$Genes[["rep_p_median"]]))
  expect_true(is.numeric(mydtu$Genes[["rep_p_min"]]))
  expect_true(is.numeric(mydtu$Genes[["rep_p_max"]]))
  expect_true(is.numeric(mydtu$Genes[["rep_na_freq"]]))
  expect_true(is.logical(mydtu$Genes[["rep_reprod"]]))
  expect_true(is.logical(mydtu$Genes[["elig"]]))
  expect_true(is.logical(mydtu$Genes[["elig_fx"]]))
  expect_true(is.logical(mydtu$Genes[["sig"]]))
  expect_true(is.logical(mydtu$Genes[["DTU"]]))
  expect_true(is.logical(mydtu$Genes[["transc_DTU"]]))
  
  expect_true(is.data.frame(mydtu$Transcripts))
  expect_equal(dim(mydtu$Transcripts)[2], 42)
  expect_named(mydtu$Transcripts, c("target_id", "parent_id", "elig_xp", "elig", "sig", "elig_fx", "quant_reprod", "rep_reprod", "DTU", "gene_DTU", 
                                    "meanA", "meanB", "stdevA", "stdevB", "sumA", "sumB", "log2FC", "totalA", "totalB", "propA", "propB", "Dprop", "pval", "pval_corr", 
                                    "quant_p_median", "quant_p_min","quant_p_max", "quant_Dprop_mean", "quant_Dprop_stdev", "quant_Dprop_min","quant_Dprop_max", "quant_na_freq", "quant_dtu_freq",
                                    "rep_p_median", "rep_p_min","rep_p_max", "rep_Dprop_mean", "rep_Dprop_stdev", "rep_Dprop_min","rep_Dprop_max", "rep_na_freq", "rep_dtu_freq") )
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
  expect_true(is.numeric(mydtu$Transcripts[["log2FC"]]))
  expect_true(is.numeric(mydtu$Transcripts[["meanA"]]))
  expect_true(is.numeric(mydtu$Transcripts[["meanB"]]))
  expect_true(is.numeric(mydtu$Transcripts[["stdevA"]]))
  expect_true(is.numeric(mydtu$Transcripts[["stdevB"]]))
  expect_true(is.numeric(mydtu$Transcripts[["pval"]]))
  expect_true(is.numeric(mydtu$Transcripts[["pval_corr"]]))
  expect_true(is.numeric(mydtu$Transcripts[["quant_dtu_freq"]]))
  expect_true(is.numeric(mydtu$Transcripts[["quant_p_median"]]))
  expect_true(is.numeric(mydtu$Transcripts[["quant_p_min"]]))
  expect_true(is.numeric(mydtu$Transcripts[["quant_p_max"]]))
  expect_true(is.numeric(mydtu$Transcripts[["quant_Dprop_mean"]]))
  expect_true(is.numeric(mydtu$Transcripts[["quant_Dprop_stdev"]]))
  expect_true(is.numeric(mydtu$Transcripts[["quant_Dprop_min"]]))
  expect_true(is.numeric(mydtu$Transcripts[["quant_Dprop_max"]]))
  expect_true(is.numeric(mydtu$Transcripts[["quant_na_freq"]]))
  expect_true(is.logical(mydtu$Transcripts[["quant_reprod"]]))
  expect_true(is.numeric(mydtu$Transcripts[["rep_dtu_freq"]]))
  expect_true(is.numeric(mydtu$Transcripts[["rep_p_median"]]))
  expect_true(is.numeric(mydtu$Transcripts[["rep_p_min"]]))
  expect_true(is.numeric(mydtu$Transcripts[["rep_p_max"]]))
  expect_true(is.numeric(mydtu$Transcripts[["rep_Dprop_mean"]]))
  expect_true(is.numeric(mydtu$Transcripts[["rep_Dprop_stdev"]]))
  expect_true(is.numeric(mydtu$Transcripts[["rep_Dprop_min"]]))
  expect_true(is.numeric(mydtu$Transcripts[["rep_Dprop_max"]]))
  expect_true(is.numeric(mydtu$Transcripts[["rep_na_freq"]]))
  expect_true(is.logical(mydtu$Transcripts[["rep_reprod"]]))
  expect_true(is.list(mydtu$Abundances))
  expect_true(is.data.frame(mydtu$Abundances[[1]]))
  expect_true(is.data.frame(mydtu$Abundances[[2]]))
  
  # No quantification bootstraps.
  mydtu <- call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= "ONE", name_B= "TWO", verbose = FALSE, qboot = FALSE, lean=FALSE)  
  expect_false(any(c("quant_dtu_freq", "quant_p_median", "quant_p_min", "quant_p_max", "quant_na_freq", "quant_reprod") %in% names(mydtu$Genes)))
  expect_false(any(c("quant_dtu_freq", "quant_p_median", "quant_p_min", "quant_p_max", "quant_Dprop_mean", "quant_Dprop_stdev", "quant_Dprop_min", "quant_Dprop_max", "quant_na_freq", "quant_reprod") %in% names(mydtu$Transcripts)))
  
  # No cross-replicate bootstraps.
  mydtu <- call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= "ONE", name_B= "TWO", verbose = FALSE, rboot = FALSE, lean=FALSE)  
  expect_false(any(c("rep_dtu_freq", "rep_p_median", "rep_p_min", "rep_p_max", "rep_na_freq", "rep_reprod") %in% names(mydtu$Genes)))
  expect_false(any(c("rep_dtu_freq", "rep_p_median", "rep_p_min", "rep_p_max", "rep_Dprop_mean", "rep_Dprop_stdev", "rep_Dprop_min", "rep_Dprop_max", "rep_na_freq", "rep_reprod") %in% names(mydtu$Transcripts)))
  
  # No gene tests.
  mydtu <- call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= "ONE", name_B= "TWO", testmode="transc", qbootnum=2, rboot=TRUE, qboot=TRUE, verbose = FALSE, lean=FALSE)
  expect_false(any(c("quant_dtu_freq", "quant_p_median", "quant_p_min", "quant_p_max", "quant_na_freq", "quant_reprod",
                     "rep_dtu_freq", "rep_p_median", "rep_p_min", "rep_p_max", "rep_na_freq", "rep_reprod") %in% names(mydtu$Genes)))
  expect_true(all(c("quant_dtu_freq", "quant_p_median", "quant_p_min", "quant_p_max", "quant_Dprop_mean", "quant_Dprop_stdev", "quant_Dprop_min", "quant_Dprop_max", "quant_na_freq", "quant_reprod",
                    "rep_dtu_freq", "rep_p_median", "rep_p_min", "rep_p_max", "rep_Dprop_mean", "rep_Dprop_stdev", "rep_Dprop_min", "rep_Dprop_max", "rep_na_freq", "rep_reprod") %in% names(mydtu$Transcripts)))
  
  # No transcript tests.
  mydtu <- call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= "ONE", name_B= "TWO", testmode="genes", qbootnum=2, rboot=TRUE, qboot=TRUE, verbose = FALSE, lean=FALSE)
  expect_false(any(c("quant_dtu_freq", "quant_p_median", "quant_p_min", "quant_p_max", "quant_Dprop_mean", "quant_Dprop_stdev", "quant_Dprop_min", "quant_Dprop_max", "quant_na_freq", "quant_reprod",
                     "rep_dtu_freq", "rep_p_median", "rep_p_min", "rep_p_max", "rep_Dprop_mean", "rep_Dprop_stdev", "rep_Dprop_min", "rep_Dprop_max", "rep_na_freq", "rep_reprod") %in% names(mydtu$Transcripts)))
  expect_true(all(c("quant_dtu_freq", "quant_p_median", "quant_p_min", "quant_p_max", "quant_na_freq", "quant_reprod",
                    "rep_dtu_freq", "rep_p_median", "rep_p_min", "rep_p_max", "rep_na_freq", "rep_reprod") %in% names(mydtu$Genes)))
  
  # No bootstrap variance statstics.
  mydtu <- call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= "ONE", name_B= "TWO", qbootnum=2, rboot=TRUE, qboot=TRUE, verbose = FALSE, lean=TRUE)
  expect_false(any(c("quant_p_median", "quant_p_min", "quant_p_max", 
                     "rep_p_median", "rep_p_min", "rep_p_max") %in% names(mydtu$Genes)))
  expect_true(all(c("quant_dtu_freq", "quant_reprod", "quant_na_freq",
                    "rep_dtu_freq", "rep_reprod", "rep_na_freq") %in% names(mydtu$Genes)))
  expect_false(any(c("quant_p_median", "quant_p_min", "quant_p_max", "quant_Dprop_mean", "quant_Dprop_stdev", "quant_Dprop_min", "quant_Dprop_max",
                     "rep_p_median", "rep_p_min", "rep_p_max", "rep_Dprop_mean", "rep_Dprop_stdev", "rep_Dprop_min", "rep_Dprop_max") %in% names(mydtu$Transcripts)))
  expect_true(all(c("quant_dtu_freq", "quant_reprod", "quant_na_freq",
                    "rep_dtu_freq", "rep_reprod", "rep_na_freq") %in% names(mydtu$Transcripts)))
})

#==============================================================================
test_that("The output content is complete", {
  sim <- sim_boot_data()
  # Emulate non-sleuth bootstrap data.
  data_A <- sim$boots_A
  data_B <- sim$boots_B
  # Emulate non-bootstrap data.
  counts_A <- data_A[[2]]
  counts_B <- data_B[[1]]
  
  mydtu <- list(call_DTU(annot= sim$annot, count_data_A = counts_A, count_data_B = counts_B, rboot=FALSE, qboot=FALSE, verbose = FALSE, description="test", reckless=TRUE, lean=FALSE),
                call_DTU(annot= sim$annot, boot_data_A = data_A, boot_data_B = data_B, rboot=TRUE, qboot=TRUE, qbootnum=2, verbose = FALSE, description="test", reckless=TRUE, lean=FALSE)
                )
  
  for (x in length(mydtu)) {
    # Parameters.
    expect_false(is.na(mydtu[[x]]$Parameters$"var_name"))
    expect_false(is.na(mydtu[[x]]$Parameters$"cond_A"))
    expect_false(is.na(mydtu[[x]]$Parameters$"cond_B"))
    expect_false(is.na(mydtu[[x]]$Parameters$"num_replic_A"))
    expect_false(is.na(mydtu[[x]]$Parameters$"num_replic_B"))
    expect_false(is.na(mydtu[[x]]$Parameters$"p_thresh"))
    expect_false(is.na(mydtu[[x]]$Parameters$"dprop_thresh"))
    expect_false(is.na(mydtu[[x]]$Parameters$"abund_thresh"))
    expect_false(is.na(mydtu[[x]]$Parameters$"tests"))
    expect_false(is.na(mydtu[[x]]$Parameters$"data_type"))
    expect_false(is.na(mydtu[[x]]$Parameters$"num_genes"))
    expect_false(is.na(mydtu[[x]]$Parameters$"num_transc"))
    expect_false(is.na(mydtu[[x]]$Parameters$"rats_version"))
    expect_false(any(is.na(mydtu[[x]]$Parameters$"R_version")))
    expect_false(is.na(mydtu[[x]]$Parameters$"description"))
    expect_false(is.na(mydtu[[x]]$Parameters$"time"))
    expect_false(is.na(mydtu[[x]]$Parameters$"rep_boot"))
    expect_false(is.na(mydtu[[x]]$Parameters$"quant_boot"))
    # expect_false(is.na(mydtu[[x]]$Parameters$"seed"))  # Actually it may be NA, if no specific seed is given.
    expect_false(is.na(mydtu[[x]]$Parameters$"reckless"))
    if (x>1) {
      expect_false(is.na(mydtu[[x]]$Parameters$"rep_bootnum"))
      expect_false(is.na(mydtu[[x]]$Parameters$"quant_bootnum"))
      expect_false(is.na(mydtu[[x]]$Parameters$"rep_reprod_thresh"))
      expect_false(is.na(mydtu[[x]]$Parameters$"quant_reprod_thresh"))
    }
    # Genes.
    expect_false(all(is.na(mydtu[[x]]$Genes[["parent_id"]])))
    expect_false(all(is.na(mydtu[[x]]$Genes[["known_transc"]])))
    expect_false(all(is.na(mydtu[[x]]$Genes[["detect_transc"]])))
    expect_false(all(is.na(mydtu[[x]]$Genes[["pval"]])))
    expect_false(all(is.na(mydtu[[x]]$Genes[["pval_corr"]])))
    expect_false(all(is.na(mydtu[[x]]$Genes[["elig"]])))
    expect_false(all(is.na(mydtu[[x]]$Genes[["elig_fx"]])))
    expect_false(all(is.na(mydtu[[x]]$Genes[["sig"]])))
    expect_false(all(is.na(mydtu[[x]]$Genes[["DTU"]])))
    expect_false(all(is.na(mydtu[[x]]$Genes[["transc_DTU"]])))
    if (x>1) {
      expect_false(all(is.na(mydtu[[x]]$Genes[["quant_dtu_freq"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["quant_p_median"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["quant_p_min"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["quant_p_max"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["quant_na_freq"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["quant_reprod"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["rep_dtu_freq"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["rep_p_median"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["rep_p_min"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["rep_p_max"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["rep_na_freq"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["rep_reprod"]])))
    }
    # Transcripts.
    expect_false(all(is.na(mydtu[[x]]$Transcripts[["target_id"]])))
    expect_false(all(is.na(mydtu[[x]]$Transcripts[["parent_id"]])))
    expect_false(all(is.na(mydtu[[x]]$Transcripts[["DTU"]])))
    expect_false(all(is.na(mydtu[[x]]$Transcripts[["gene_DTU"]])))
    expect_false(all(is.na(mydtu[[x]]$Transcripts[["meanA"]])))
    expect_false(all(is.na(mydtu[[x]]$Transcripts[["meanB"]])))
    expect_false(all(is.na(mydtu[[x]]$Transcripts[["stdevA"]])))
    expect_false(all(is.na(mydtu[[x]]$Transcripts[["stdevB"]])))
    expect_false(all(is.na(mydtu[[x]]$Transcripts[["sumA"]])))
    expect_false(all(is.na(mydtu[[x]]$Transcripts[["sumB"]])))
    expect_false(all(is.na(mydtu[[x]]$Transcripts[["elig_xp"]])))
    expect_false(all(is.na(mydtu[[x]]$Transcripts[["elig"]])))
    expect_false(all(is.na(mydtu[[x]]$Transcripts[["propA"]])))
    expect_false(all(is.na(mydtu[[x]]$Transcripts[["propB"]])))
    expect_false(all(is.na(mydtu[[x]]$Transcripts[["Dprop"]])))
    expect_false(all(is.na(mydtu[[x]]$Transcripts[["elig_fx"]])))
    expect_false(all(is.na(mydtu[[x]]$Transcripts[["pval"]])))
    expect_false(all(is.na(mydtu[[x]]$Transcripts[["pval_corr"]])))
    expect_false(all(is.na(mydtu[[x]]$Transcripts[["sig"]])))
    if (x>1) {
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["quant_dtu_freq"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["quant_p_median"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["quant_p_min"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["quant_p_max"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["quant_Dprop_mean"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["quant_Dprop_stdev"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["quant_Dprop_min"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["quant_Dprop_max"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["quant_na_freq"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["quant_reprod"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["rep_dtu_freq"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["rep_p_median"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["rep_p_min"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["rep_p_max"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["rep_Dprop_mean"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["rep_Dprop_stdev"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["rep_Dprop_min"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["rep_Dprop_max"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["rep_na_freq"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["rep_reprod"]])))
    }
    # Paranoid check that all the annotation entries are present in the output.
    expect_equal(dim(mydtu[[x]]$Genes)[1], length(unique(sim$annot$parent_id)))
    expect_equal(dim(mydtu[[x]]$Transcripts)[1], length(sim$annot$target_id))
  }
})

#==============================================================================
test_that("The result is consistent across input data formats", {
  sim <- sim_boot_data()
  data_A <- sim$boots_A
  data_B <- sim$boots_B
  # Emulate non-bootstrap data.
  counts_A <- data_A[[2]]
  counts_B <- data_B[[1]]
  
  mydtu <- list(call_DTU(annot= sim$annot, boot_data_A = data_A, boot_data_B = data_B, rboot=FALSE, qboot=FALSE, verbose = FALSE, reckless=TRUE),
                call_DTU(annot= sim$annot, count_data_A = counts_A, count_data_B = counts_B, rboot=FALSE, qboot=FALSE, verbose = FALSE, reckless=TRUE))
  
  expect_equal(mydtu[[1]][[2]][, seq(1,7), with=FALSE], mydtu[[2]][[2]][, seq(1,7), with=FALSE])
  expect_equal(mydtu[[1]][[3]][, seq(1,4), with=FALSE], mydtu[[2]][[3]][, seq(1,4), with=FALSE])
})

#==============================================================================
test_that("Filters work correctly", {
  
  # !!! This test is tightly dependent on the data used for the test and the parameter values, in order to
  # !!! ensure correct response to specific scenarios.
  
  sim <- sim_boot_data()
  mydtu <- call_DTU(annot= sim$annot, boot_data_A=sim$boots_A, boot_data_B=sim$boots_B, name_A= "ONE", name_B= "TWO", abund_thresh=8, dprop_thresh=0.1, verbose = FALSE, rboot=FALSE, qboot = FALSE, seed=666, reckless=TRUE)
  
  expect_equivalent(as.list(mydtu$Genes["1A1B", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(2, 2, 2, TRUE, TRUE))
  expect_equivalent(as.list(mydtu$Genes["1A1N", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(1, 1, 0, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Genes["1B1C", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(2, 1, 0, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Genes["1D1C", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(2, 1, 0, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Genes["ALLA", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(1, 1, 0, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Genes["ALLB", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(2, 2, 0, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Genes["CC", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(2, 2, 2, TRUE, TRUE))
  expect_equivalent(as.list(mydtu$Genes["LC", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(2, 2, 1, FALSE, TRUE))
  expect_equivalent(as.list(mydtu$Genes["MIX6", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(6, 5, 5, TRUE, TRUE))
  expect_equivalent(as.list(mydtu$Genes["NIB", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(1, 0, 0, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Genes["NN", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(2, 2, 2, TRUE, FALSE))
  
  setkey(mydtu$Transcripts, target_id)
  expect_equivalent(as.list(mydtu$Transcripts["1A1B.a", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, TRUE))
  expect_equivalent(as.list(mydtu$Transcripts["1A1B.b", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, TRUE))
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
                    list(TRUE, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Transcripts["ALLB1", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Transcripts["ALLB2", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, FALSE, FALSE))
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
                    list(FALSE, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Transcripts["1NN", .(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, FALSE))
  expect_equivalent(as.list(mydtu$Transcripts["2NN", .(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, FALSE))
})


