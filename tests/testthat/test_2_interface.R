#==============================================================================
#==============================================================================
context("DTU Input checks.")

#==============================================================================
test_that("The input checks work", {
  name_A <- "one"
  name_B <- "two"
  wrong_name <- "RUBBISH_COLUMN_NAME"
  
  # Emulate sleuth and annotation.
  sim1 <- sim_sleuth_data(cnames=c(name_A, name_B))
  sim2 <- sim_sleuth_data(varname= "waffles", COUNTS_COL= "counts", TARGET_COL= "target" , PARENT_COL= "parent", 
                         BS_TARGET_COL= "id", cnames= c("AAAA","BBBB"))
  # Emulate non-sleuth bootstrap data.
  data_A <- denest_sleuth_boots(sim1$slo, sim1$annot$target_id, c(1,3), "est_counts", "target_id")
  data_B <- denest_sleuth_boots(sim1$slo, sim1$annot$target_id, c(2,4), "est_counts", "target_id")
  # Emulate non-bootstrap data.
  counts_A <- data_A[[1]]
  counts_B <- data_B[[1]]
  
  # No false alarms with valid calls.
  expect_silent(call_DTU(annot= sim2$annot, slo= sim2$slo, name_A = "AAAA", name_B = "BBBB", varname= "waffles", p_thresh= 0.01, count_thresh= 10,
                         conf_thresh = 0.8, testmode= "transc", correction= "bonferroni", verbose= FALSE, boots= "genes",
                         bootnum= 2, COUNTS_COL= "counts", TARGET_COL= "target", 
                         PARENT_COL= "parent", BS_TARGET_COL= "id"))
  expect_silent(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= name_B, verbose = FALSE))
  expect_silent(call_DTU(annot= sim1$annot, boot_data_A= data_A, boot_data_B= data_B, verbose = FALSE))
  expect_silent(call_DTU(annot= sim1$annot, count_data_A= counts_A, count_data_B= counts_B, boots= "none", verbose = FALSE))
  expect_silent(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= name_B, verbose = FALSE,
                         boot_data_A= data_A, boot_data_B= data_B, count_data_A= counts_A, count_data_B= counts_B))
  
  # Annotation is not a dataframe.
  expect_error(call_DTU(annot= list("not", "a", "dataframe"), slo= sim1$slo, name_A= name_A, name_B= name_B, verbose = FALSE), "annot is not a data.frame")
  # Annotation field names.
  expect_error(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= name_B, TARGET_COL= wrong_name, verbose = FALSE),
               "target and/or parent IDs field names do not exist in annot", fixed= TRUE)
  expect_error(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= name_B, PARENT_COL= wrong_name, verbose = FALSE),
               "target and/or parent IDs field names do not exist in annot", fixed= TRUE)
  # Inconsistent annotation.
  sim3 <- sim_sleuth_data(errannot_inconsistent= TRUE, cnames= c(name_A, name_B))
  expect_error(call_DTU(annot= sim3$annot, slo= sim3$slo, name_A= name_A, name_B= name_B, verbose = FALSE),
               "Inconsistent set of transcript IDs", fixed= TRUE)
  # Non unique IDs.
  a <- copy(sim1$annot)
  a[1, "target_id"] <- a[2, "target_id"]
  expect_error(call_DTU(annot= a, slo= sim1$slo, name_A= name_A, name_B= name_B, verbose = FALSE), "transcript identifiers are not unique")
  
  # Boot data is not a list of datatables.
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= c("not", "a", "list"), boot_data_B= data_B, verbose = FALSE), "bootstrap data are not lists")
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= data_A, boot_data_B= c("not", "a", "list"), verbose = FALSE), "bootstrap data are not lists")
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= list( list("not"), list("a"), list("table")), boot_data_B= data_B, verbose = FALSE), "bootstrap data are not lists")
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= data_A, boot_data_B= list( list("not"), list("a"), list("table")), verbose = FALSE), "bootstrap data are not lists")
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= list(data.frame(a="not", b="a", c="table")), boot_data_B= data_B, verbose = FALSE), "not lists of data.tables")
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= data_A, boot_data_B= list(data.frame(a="not", b="a", c="table")), verbose = FALSE), "not lists of data.tables")
  
  # Counts data is not datatables.
  expect_error(call_DTU(annot= sim1$annot, count_data_A= c("not", "a", "list"), count_data_B= counts_B, verbose = FALSE), "counts data are not data.tables")
  expect_error(call_DTU(annot= sim1$annot, count_data_A= counts_A, count_data_B= c("not", "a", "list"), verbose = FALSE), "counts data are not data.tables")
  expect_error(call_DTU(annot= sim1$annot, count_data_A= data.frame(a="not", b="a", c="table"), count_data_B= counts_B, verbose = FALSE), "counts data are not data.tables")
  expect_error(call_DTU(annot= sim1$annot, count_data_A= counts_A, count_data_B= data.frame(a="not", b="a", c="table"), verbose = FALSE), "counts data are not data.tables")
  
  # Bootstrap field names.
  expect_error(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= name_B, BS_TARGET_COL= wrong_name, verbose = FALSE),
               "target IDs field name does not exist in the bootstraps", fixed= TRUE)
  expect_error(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= name_B, COUNTS_COL= wrong_name, verbose = FALSE),
               "counts field name does not exist", fixed= TRUE)
  # Number of bootstraps.
  expect_error(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= name_B, bootnum = -5, verbose = FALSE),
               "Invalid number of bootstraps", fixed= TRUE)
  # Bootstraps without boot data.
  expect_error(call_DTU(annot= sim1$annot, count_data_A= counts_A, count_data_B= counts_B, boots="both", bootnum=2, verbose= FALSE ), 
               "No bootstrapped estimates", fixed= TRUE)
  
  # Correction method.
  expect_error(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= name_B, correction= wrong_name, verbose = FALSE),
               "Invalid p-value correction method name", fixed= TRUE)
  
  # Covariate name.
  expect_error(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= name_B, varname= wrong_name, verbose = FALSE),
               "covariate name does not exist", fixed= TRUE)
  
  # Condition names.
  expect_error(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= wrong_name, name_B= name_B, verbose = FALSE),
               "conditions do not exist", fixed= TRUE)
  expect_error(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= wrong_name, verbose = FALSE),
               "conditions do not exist", fixed= TRUE)
  
  # Verbose is bool.
  expect_error(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= name_B, verbose="yes"),
               "not interpretable as logical", fixed= TRUE)
  
  # Tests.
  expect_error(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= name_B, testmode="GCSE", verbose = FALSE),
               "Unrecognized value for testmode", fixed= TRUE)
  expect_silent(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= name_B, testmode="genes", verbose = FALSE, boots = "none"))
  expect_silent(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= name_B, testmode="transc", verbose = FALSE, boots = "none"))
  
  expect_error(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= name_B, boots="GCSE", verbose = FALSE),
               "Unrecognized value for boots", fixed= TRUE)
  expect_silent(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= name_B, boots="genes", bootnum = 2, verbose = FALSE))
  expect_silent(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= name_B, boots="transc", bootnum = 2, verbose = FALSE))
  
  # Probability threshold.
  expect_error(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= name_B, p_thresh = 666, verbose = FALSE),
               "Invalid p-value threshold", fixed= TRUE)
  expect_error(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= name_B, p_thresh = -0.05, verbose = FALSE),
               "Invalid p-value threshold", fixed= TRUE)
  
  # Read counts threshold.
  expect_error(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= name_B, count_thresh = -5, verbose = FALSE),
               "Invalid read-count threshold", fixed= TRUE)
  
  # Proportion change threshold.
  expect_error(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= name_B, dprop_thresh = -2, verbose = FALSE),
               "Invalid proportion difference threshold", fixed= TRUE)
  expect_error(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= name_B, dprop_thresh = 2, verbose = FALSE),
               "Invalid proportion difference threshold", fixed= TRUE)
  
  # Confidence threshold.
  expect_error(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= name_B, conf_thresh = -2, verbose = FALSE),
               "Invalid confidence threshold", fixed= TRUE)
  expect_error(call_DTU(annot= sim1$annot, slo= sim1$slo, name_A= name_A, name_B= name_B, conf_thresh = 2, verbose = FALSE),
               "Invalid confidence threshold", fixed= TRUE)
  
})


#==============================================================================
#==============================================================================
context("DTU Output")

#==============================================================================
test_that("The output structure is correct", {
  sim <- sim_sleuth_data(cnames=c("ONE","TWO"))
  mydtu <- call_DTU(annot= sim$annot, slo= sim$slo, name_A= "ONE", name_B= "TWO", bootnum=2, verbose = FALSE)
  
  expect_type(mydtu, "list")
  expect_equal(length(mydtu), 3)
  expect_named(mydtu, c("Parameters", "Genes", "Transcripts"))
  
  expect_type(mydtu$Parameters, "list")
  expect_length(mydtu$Parameters, 18)
  expect_named(mydtu$Parameters, c("var_name", "cond_A", "cond_B", "num_replic_A", "num_replic_B", "p_thresh", 
                                   "count_thresh", "dprop_thresh", "conf_thresh", "tests", "bootstrap", "bootnum",
                                   "data_type", "num_genes", "num_transc", "description", 
                                   "rats_version", "R_version"))
  
  expect_true(is.data.frame(mydtu$Genes))
  expect_equal(dim(mydtu$Genes)[2], 24)
  expect_named(mydtu$Genes, c("parent_id", "DTU", "transc_DTU", "known_transc", "detect_transc", "elig_transc",  
                              "elig", "elig_fx", "pvalAB", "pvalBA", "pvalAB_corr", "pvalBA_corr", "sig", "boot_dtu_freq", "conf",
                              "boot_p_meanAB", "boot_p_meanBA", "boot_p_stdevAB", "boot_p_stdevBA", "boot_p_minAB", "boot_p_minBA", 
                              "boot_p_maxAB", "boot_p_maxBA", "boot_na"))
  expect_true(is.numeric(mydtu$Genes[["known_transc"]]))
  expect_true(is.numeric(mydtu$Genes[["detect_transc"]]))
  expect_true(is.numeric(mydtu$Genes[["pvalAB"]]))
  expect_true(is.numeric(mydtu$Genes[["pvalBA"]]))
  expect_true(is.numeric(mydtu$Genes[["pvalAB_corr"]]))
  expect_true(is.numeric(mydtu$Genes[["pvalBA_corr"]]))
  expect_true(is.numeric(mydtu$Genes[["boot_dtu_freq"]]))
  expect_true(is.numeric(mydtu$Genes[["boot_p_meanAB"]]))
  expect_true(is.numeric(mydtu$Genes[["boot_p_meanBA"]]))
  expect_true(is.numeric(mydtu$Genes[["boot_p_stdevAB"]]))
  expect_true(is.numeric(mydtu$Genes[["boot_p_stdevBA"]]))
  expect_true(is.numeric(mydtu$Genes[["boot_p_minAB"]]))
  expect_true(is.numeric(mydtu$Genes[["boot_p_minBA"]]))
  expect_true(is.numeric(mydtu$Genes[["boot_p_maxAB"]]))
  expect_true(is.numeric(mydtu$Genes[["boot_p_maxBA"]]))
  expect_true(is.numeric(mydtu$Genes[["boot_na"]]))
  expect_true(is.logical(mydtu$Genes[["elig"]]))
  expect_true(is.logical(mydtu$Genes[["elig_fx"]]))
  expect_true(is.logical(mydtu$Genes[["sig"]]))
  expect_true(is.logical(mydtu$Genes[["DTU"]]))
  expect_true(is.logical(mydtu$Genes[["conf"]]))
  expect_true(is.logical(mydtu$Genes[["transc_DTU"]]))
  
  expect_true(is.data.frame(mydtu$Transcripts))
  expect_equal(dim(mydtu$Transcripts)[2], 28)
  expect_named(mydtu$Transcripts, c("target_id", "parent_id", "DTU", "gene_DTU", "meanA", "meanB", "stdevA", "stdevB",
                                    "sumA", "sumB", "totalA", "totalB", "elig_xp", "elig", "propA", "propB", "Dprop", 
                                    "elig_fx", "pval", "pval_corr", "sig", "boot_dtu_freq", "conf", "boot_p_mean", "boot_p_stdev", 
                                    "boot_p_min", "boot_p_max", "boot_na"))
  expect_true(is.logical(mydtu$Transcripts[["elig_xp"]]))
  expect_true(is.logical(mydtu$Transcripts[["elig"]]))
  expect_true(is.logical(mydtu$Transcripts[["elig_fx"]]))
  expect_true(is.logical(mydtu$Transcripts[["sig"]]))
  expect_true(is.logical(mydtu$Transcripts[["DTU"]]))
  expect_true(is.logical(mydtu$Transcripts[["conf"]]))
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
  expect_true(is.numeric(mydtu$Transcripts[["boot_dtu_freq"]]))
  expect_true(is.numeric(mydtu$Transcripts[["boot_p_mean"]]))
  expect_true(is.numeric(mydtu$Transcripts[["boot_p_stdev"]]))
  expect_true(is.numeric(mydtu$Transcripts[["boot_p_min"]]))
  expect_true(is.numeric(mydtu$Transcripts[["boot_p_max"]]))
  expect_true(is.numeric(mydtu$Transcripts[["boot_na"]]))
  
  mydtu <- call_DTU(annot= sim$annot, slo= sim$slo, name_A= "ONE", name_B= "TWO", verbose = FALSE, boots = "none")  
  expect_false(any(c("boot_dtu_freq", "boot_p_mean", "boot_p_stdev", "boot_p_min", "boot_p_max", "boot_na") %in% names(mydtu$Genes)))
  expect_false(any(c("boot_dtu_freq", "boot_p_mean", "boot_p_stdev", "boot_p_min", "boot_p_max", "boot_na") %in% names(mydtu$Transcripts)))
  
  mydtu <- call_DTU(annot= sim$annot, slo= sim$slo, name_A= "ONE", name_B= "TWO", boots="transc", bootnum=2, verbose = FALSE)  
  expect_false( any(c("boot_dtu_freq", "boot_p_mean", "boot_p_stdev", "boot_p_min", "boot_p_max", "boot_na") %in% names(mydtu$Genes)))
  
  mydtu <- call_DTU(annot= sim$annot, slo= sim$slo, name_A= "ONE", name_B= "TWO", boots="genes", bootnum=2, verbose = FALSE)  
  expect_false(any(c("boot_dtu_freq", "boot_p_mean", "boot_p_stdev", "boot_p_min", "boot_p_max", "boot_na") %in% names(mydtu$Transcripts)))
})

#==============================================================================
test_that("The output content is complete", {
  sim <- sim_sleuth_data(cnames=c("ONE","TWO"))
  # Emulate non-sleuth bootstrap data.
  data_A <- denest_sleuth_boots(sim$slo, sim$annot$target_id, c(1,3), "est_counts", "target_id")
  data_B <- denest_sleuth_boots(sim$slo, sim$annot$target_id, c(2,4), "est_counts", "target_id")
  # Emulate non-bootstrap data.
  counts_A <- data_A[[1]]
  counts_B <- data_B[[2]]
  
  mydtu <- list(call_DTU(annot= sim$annot, slo= sim$slo, name_A= "ONE", name_B= "TWO", boots="both", bootnum=2, verbose = FALSE),
                call_DTU(annot= sim$annot, boot_data_A = data_A, boot_data_B = data_B, boots="both", bootnum=2, verbose = FALSE),
                call_DTU(annot= sim$annot, count_data_A = counts_A, count_data_B = counts_B, boots="none", verbose = FALSE))
  
  for (x in length(mydtu)) {
    # Parameters.
    expect_false(is.na(mydtu[[x]]$Parameters$"var_name"))
    expect_false(is.na(mydtu[[x]]$Parameters$"cond_A"))
    expect_false(is.na(mydtu[[x]]$Parameters$"cond_B"))
    expect_false(is.na(mydtu[[x]]$Parameters$"num_replic_A"))
    expect_false(is.na(mydtu[[x]]$Parameters$"num_replic_B"))
    expect_false(is.na(mydtu[[x]]$Parameters$"p_thresh"))
    expect_false(is.na(mydtu[[x]]$Parameters$"dprop_thresh"))
    expect_false(is.na(mydtu[[x]]$Parameters$"count_thresh"))
    expect_false(is.na(mydtu[[x]]$Parameters$"conf_thresh"))
    expect_false(is.na(mydtu[[x]]$Parameters$"tests"))
    expect_false(is.na(mydtu[[x]]$Parameters$"bootstrap"))
    expect_false(is.na(mydtu[[x]]$Parameters$"bootnum"))
    expect_false(is.na(mydtu[[x]]$Parameters$"data_type"))
    expect_false(is.na(mydtu[[x]]$Parameters$"num_genes"))
    expect_false(is.na(mydtu[[x]]$Parameters$"num_transc"))
    expect_false(is.na(mydtu[[x]]$Parameters$"rats_version"))
    expect_false(any(is.na(mydtu[[x]]$Parameters$"R_version")))
    # Genes.
    expect_false(all(is.na(mydtu[[x]]$Genes[["parent_id"]])))
    expect_false(all(is.na(mydtu[[x]]$Genes[["known_transc"]])))
    expect_false(all(is.na(mydtu[[x]]$Genes[["detect_transc"]])))
    expect_false(all(is.na(mydtu[[x]]$Genes[["pvalAB"]])))
    expect_false(all(is.na(mydtu[[x]]$Genes[["pvalBA"]])))
    expect_false(all(is.na(mydtu[[x]]$Genes[["pvalAB_corr"]])))
    expect_false(all(is.na(mydtu[[x]]$Genes[["pvalBA_corr"]])))
    expect_false(all(is.na(mydtu[[x]]$Genes[["elig"]])))
    expect_false(all(is.na(mydtu[[x]]$Genes[["elig_fx"]])))
    expect_false(all(is.na(mydtu[[x]]$Genes[["sig"]])))
    expect_false(all(is.na(mydtu[[x]]$Genes[["DTU"]])))
    expect_false(all(is.na(mydtu[[x]]$Genes[["transc_DTU"]])))
    if (x <3) {
      expect_false(all(is.na(mydtu[[x]]$Genes[["boot_dtu_freq"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["boot_p_meanAB"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["boot_p_meanBA"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["boot_p_stdevAB"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["boot_p_stdevBA"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["boot_p_minAB"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["boot_p_minBA"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["boot_p_maxAB"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["boot_p_maxBA"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["boot_na"]])))
      expect_false(all(is.na(mydtu[[x]]$Genes[["conf"]])))
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
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["boot_dtu_freq"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["boot_p_mean"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["boot_p_stdev"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["boot_p_min"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["boot_p_max"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["boot_na"]])))
      expect_false(all(is.na(mydtu[[x]]$Transcripts[["conf"]])))
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
  data_A <- denest_sleuth_boots(sim$slo, sim$annot$target_id, c(1,3), "est_counts", "target_id")
  data_B <- denest_sleuth_boots(sim$slo, sim$annot$target_id, c(2,4), "est_counts", "target_id")
  # Emulate non-bootstrap data.
  counts_A <- data_A[[1]]
  counts_B <- data_B[[1]]
  
  mydtu <- list(call_DTU(annot= sim$annot, slo= sim$slo, name_A= "ONE", name_B= "TWO", boots="none", verbose = FALSE),
                call_DTU(annot= sim$annot, boot_data_A = data_A, boot_data_B = data_B, boots="none", verbose = FALSE),
                call_DTU(annot= sim$annot, count_data_A = counts_A, count_data_B = counts_B, boots="none", verbose = FALSE))
  
  expect_equal(mydtu[[1]][[2]], mydtu[[2]][[2]])
  expect_equal(mydtu[[1]][[2]][, seq(1,7), with=FALSE], mydtu[[3]][[2]][, seq(1,7), with=FALSE])
  expect_equal(mydtu[[1]][[3]], mydtu[[2]][[3]])
  expect_equal(mydtu[[1]][[3]][, seq(1,4), with=FALSE], mydtu[[3]][[3]][, seq(1,4), with=FALSE])
})

#==============================================================================
test_that("The summaries work", {
  sim <- sim_sleuth_data(cnames=c("ONE","TWO"))
  mydtu <- call_DTU(annot= sim$annot, slo= sim$slo, name_A= "ONE", name_B= "TWO", boots="both", bootnum=2, verbose = FALSE)
  
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
