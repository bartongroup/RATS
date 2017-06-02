#==============================================================================
#==============================================================================
context("DTU Output")

#==============================================================================
test_that("The output structure is correct", {
  sim <- sim_sleuth_data(cnames=c("ONE","TWO"))
  mydtu <- call_DTU(annot= sim$annot, slo= sim$slo, name_A= "ONE", name_B= "TWO", qbootnum=2, qboot=TRUE, rboot=TRUE, verbose = FALSE)
  
  expect_type(mydtu, "list")
  expect_equal(length(mydtu), 4)
  expect_named(mydtu, c("Parameters", "Genes", "Transcripts", "Abundances"))
  
  expect_type(mydtu$Parameters, "list")
  expect_length(mydtu$Parameters, 24)
  expect_named(mydtu$Parameters, c("description", "time", "rats_version", "R_version",
                                   "var_name", "cond_A", "cond_B", "data_type", "num_replic_A", "num_replic_B", "num_genes", "num_transc",
                                   "tests", "p_thresh", "abund_thresh", "dprop_thresh", "abund_scaling",
                                   "quant_reprod_thresh", "quant_boot", "quant_bootnum",
                                   "rep_reprod_thresh", "rep_boot", "rep_bootnum", "rep_reprod_as_crit"))
  
  expect_true(is.data.frame(mydtu$Genes))
  expect_equal(dim(mydtu$Genes)[2], 26)
  expect_named(mydtu$Genes, c("parent_id", "elig", "sig", "elig_fx", "quant_reprod", "rep_reprod", "DTU", "transc_DTU",
                              "known_transc", "detect_transc", "elig_transc", "maxDprop", "pval", "pval_corr", 
                              "quant_p_mean", "quant_p_stdev", "quant_p_min", "quant_p_max", "quant_na_freq", "quant_dtu_freq",
                              "rep_p_mean",  "rep_p_stdev", "rep_p_min", "rep_p_max", "rep_na_freq", "rep_dtu_freq") )
  expect_true(is.numeric(mydtu$Genes[["known_transc"]]))
  expect_true(is.numeric(mydtu$Genes[["detect_transc"]]))
  expect_true(is.numeric(mydtu$Genes[["maxDprop"]]))
  expect_true(is.numeric(mydtu$Genes[["pval"]]))
  expect_true(is.numeric(mydtu$Genes[["pval_corr"]]))
  expect_true(is.numeric(mydtu$Genes[["quant_dtu_freq"]]))
  expect_true(is.numeric(mydtu$Genes[["quant_p_mean"]]))
  expect_true(is.numeric(mydtu$Genes[["quant_p_stdev"]]))
  expect_true(is.numeric(mydtu$Genes[["quant_p_min"]]))
  expect_true(is.numeric(mydtu$Genes[["quant_p_max"]]))
  expect_true(is.numeric(mydtu$Genes[["quant_na_freq"]]))
  expect_true(is.logical(mydtu$Genes[["quant_reprod"]]))
  expect_true(is.numeric(mydtu$Genes[["rep_dtu_freq"]]))
  expect_true(is.numeric(mydtu$Genes[["rep_p_mean"]]))
  expect_true(is.numeric(mydtu$Genes[["rep_p_stdev"]]))
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
  expect_equal(dim(mydtu$Transcripts)[2], 36)
  expect_named(mydtu$Transcripts, c("target_id", "parent_id", "elig_xp", "elig", "sig", "elig_fx", "quant_reprod", "rep_reprod", "DTU", "gene_DTU", 
                                    "meanA", "meanB", "stdevA", "stdevB", "sumA", "sumB", "log2FC", "totalA", "totalB", "propA", "propB", "Dprop", "pval", "pval_corr", 
                                    "quant_p_mean", "quant_p_stdev", "quant_p_min","quant_p_max", "quant_na_freq", "quant_dtu_freq",
                                    "rep_p_mean", "rep_p_stdev", "rep_p_min","rep_p_max", "rep_na_freq", "rep_dtu_freq") )
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
  expect_true(is.numeric(mydtu$Transcripts[["quant_p_mean"]]))
  expect_true(is.numeric(mydtu$Transcripts[["quant_p_stdev"]]))
  expect_true(is.numeric(mydtu$Transcripts[["quant_p_min"]]))
  expect_true(is.numeric(mydtu$Transcripts[["quant_p_max"]]))
  expect_true(is.numeric(mydtu$Transcripts[["quant_na_freq"]]))
  expect_true(is.logical(mydtu$Transcripts[["quant_reprod"]]))
  expect_true(is.numeric(mydtu$Transcripts[["rep_dtu_freq"]]))
  expect_true(is.numeric(mydtu$Transcripts[["rep_p_mean"]]))
  expect_true(is.numeric(mydtu$Transcripts[["rep_p_stdev"]]))
  expect_true(is.numeric(mydtu$Transcripts[["rep_p_min"]]))
  expect_true(is.numeric(mydtu$Transcripts[["rep_p_max"]]))
  expect_true(is.numeric(mydtu$Transcripts[["rep_na_freq"]]))
  expect_true(is.logical(mydtu$Transcripts[["rep_reprod"]]))
  
  expect_true(is.list(mydtu$Abundances))
  expect_true(is.data.frame(mydtu$Abundances[[1]]))
  expect_true(is.data.frame(mydtu$Abundances[[2]]))
  
  mydtu <- call_DTU(annot= sim$annot, slo= sim$slo, name_A= "ONE", name_B= "TWO", verbose = FALSE, qboot = FALSE)  
  expect_false(any(c("quant_dtu_freq", "quant_p_mean", "quant_p_stdev", "quant_p_min", "quant_p_max", "quant_na_freq", "quant_reprod") %in% names(mydtu$Genes)))
  expect_false(any(c("quant_dtu_freq", "quant_p_mean", "quant_p_stdev", "quant_p_min", "quant_p_max", "quant_na_freq", "quant_reprod") %in% names(mydtu$Transcripts)))
  
  mydtu <- call_DTU(annot= sim$annot, slo= sim$slo, name_A= "ONE", name_B= "TWO", verbose = FALSE, rboot = FALSE)  
  expect_false(any(c("rep_dtu_freq", "rep_p_mean", "rep_p_stdev", "rep_p_min", "rep_p_max", "rep_na_freq", "rep_reprod") %in% names(mydtu$Genes)))
  expect_false(any(c("rep_dtu_freq", "rep_p_mean", "rep_p_stdev", "rep_p_min", "rep_p_max", "rep_na_freq", "rep_reprod") %in% names(mydtu$Transcripts)))
  
  
  mydtu <- call_DTU(annot= sim$annot, slo= sim$slo, name_A= "ONE", name_B= "TWO", testmode="transc", qbootnum=2, rboot=TRUE, qboot=TRUE, verbose = FALSE)
  expect_false(any(c("quant_dtu_freq", "quant_p_mean", "quant_p_stdev", "quant_p_min", "quant_p_max", "quant_na_freq", "quant_reprod",
                     "rep_dtu_freq", "rep_p_mean", "rep_p_stdev", "rep_p_min", "rep_p_max", "rep_na_freq", "rep_reprod") %in% names(mydtu$Genes)))
  expect_true(all(c("quant_dtu_freq", "quant_p_mean", "quant_p_stdev", "quant_p_min", "quant_p_max", "quant_na_freq", "quant_reprod",
                    "rep_dtu_freq", "rep_p_mean", "rep_p_stdev", "rep_p_min", "rep_p_max", "rep_na_freq", "rep_reprod") %in% names(mydtu$Transcripts)))
  
  mydtu <- call_DTU(annot= sim$annot, slo= sim$slo, name_A= "ONE", name_B= "TWO", testmode="genes", qbootnum=2, rboot=TRUE, qboot=TRUE, verbose = FALSE)
  expect_false(any(c("quant_dtu_freq", "quant_p_mean", "quant_p_stdev", "quant_p_min", "quant_p_max", "quant_na_freq", "quant_reprod",
                     "rep_dtu_freq", "rep_p_mean", "rep_p_stdev", "rep_p_min", "rep_p_max", "rep_na_freq", "rep_reprod") %in% names(mydtu$Transcripts)))
  expect_true(all(c("quant_dtu_freq", "quant_p_mean", "quant_p_stdev", "quant_p_min", "quant_p_max", "quant_na_freq", "quant_reprod",
                    "rep_dtu_freq", "rep_p_mean", "rep_p_stdev", "rep_p_min", "rep_p_max", "rep_na_freq", "rep_reprod") %in% names(mydtu$Genes)))
  
})

#==============================================================================
test_that("The output content is complete", {
  sim <- sim_sleuth_data(cnames=c("ONE","TWO"))
  # Emulate non-sleuth bootstrap data.
  data_A <- denest_sleuth_boots(sim$slo, sim$annot, c(1,3), "est_counts", "target_id")
  data_B <- denest_sleuth_boots(sim$slo, sim$annot, c(2,4), "est_counts", "target_id")
  # Emulate non-bootstrap data.
  counts_A <- data_A[[1]]
  counts_B <- data_B[[2]]
  
  mydtu <- list(call_DTU(annot= sim$annot, slo= sim$slo, name_A= "ONE", name_B= "TWO", rboot=TRUE, qboot=TRUE, qbootnum=2, verbose = FALSE, description="test"),
                call_DTU(annot= sim$annot, boot_data_A = data_A, boot_data_B = data_B, rboot=TRUE, qboot=TRUE, qbootnum=2, verbose = FALSE, description="test"),
                call_DTU(annot= sim$annot, count_data_A = counts_A, count_data_B = counts_B, rboot=FALSE, qboot=FALSE, verbose = FALSE, description="test"))
  
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
    if (x<3) {
      expect_false(is.na(mydtu[[x]]$Parameters$"rep_bootnum"))
      expect_false(is.na(mydtu[[x]]$Parameters$"rep_reprod_as_crit"))
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
    if (x <3) {
      expect_false(all(is.na(mydtu[[x]]$Genes[["quant_dtu_freq"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["quant_p_mean"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["quant_p_stdev"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["quant_p_min"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["quant_p_max"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["quant_na_freq"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["quant_reprod"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["rep_dtu_freq"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["rep_p_mean"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["rep_p_stdev"]])))
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
    if (x <3) {
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["quant_dtu_freq"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["quant_p_mean"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["quant_p_stdev"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["quant_p_min"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["quant_p_max"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["quant_na_freq"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["quant_reprod"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["rep_dtu_freq"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["rep_p_mean"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["rep_p_stdev"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["rep_p_min"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["rep_p_max"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["rep_na_freq"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["rep_reprod"]])))
    }
    # Paranoid check that all the annotation entries are present in the output.
    expect_equal(dim(mydtu[[x]]$Genes)[1], length(levels(sim$annot$parent_id)))
    expect_equal(dim(mydtu[[x]]$Transcripts)[1], length(levels(sim$annot$target_id)))
  }
})

#==============================================================================
test_that("The result is consistent across input data formats", {
  sim <- sim_sleuth_data(cnames=c("ONE","TWO"))
  # Emulate non-sleuth bootstrap data.
  data_A <- denest_sleuth_boots(sim$slo, sim$annot, c(1,3), "est_counts", "target_id")
  data_B <- denest_sleuth_boots(sim$slo, sim$annot, c(2,4), "est_counts", "target_id")
  # Emulate non-bootstrap data.
  counts_A <- data_A[[1]]
  counts_B <- data_B[[2]]
  
  mydtu <- list(call_DTU(annot= sim$annot, slo= sim$slo, name_A= "ONE", name_B= "TWO", rboot=FALSE, qboot=FALSE, verbose = FALSE),
                call_DTU(annot= sim$annot, boot_data_A = data_A, boot_data_B = data_B, rboot=FALSE, qboot=FALSE, verbose = FALSE),
                call_DTU(annot= sim$annot, count_data_A = counts_A, count_data_B = counts_B, rboot=FALSE, qboot=FALSE, verbose = FALSE))
  
  expect_equal(mydtu[[1]][[2]], mydtu[[2]][[2]])
  expect_equal(mydtu[[1]][[2]][, seq(1,7), with=FALSE], mydtu[[3]][[2]][, seq(1,7), with=FALSE])
  expect_equal(mydtu[[1]][[3]], mydtu[[2]][[3]])
  expect_equal(mydtu[[1]][[3]][, seq(1,4), with=FALSE], mydtu[[3]][[3]][, seq(1,4), with=FALSE])
})



