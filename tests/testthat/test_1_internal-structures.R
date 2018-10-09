#==============================================================================
#==============================================================================
context("DTU Internal structures")

#==============================================================================
test_that("The reporting structures are created correctly", {
  sim <- sim_boot_data()
  full <- alloc_out(sim$annot, "full")
  short <- alloc_out(sim$annot, "short")
  
  expect_type(full, "list")
  expect_type(short, "list")
  expect_equal(length(full), 4)
  expect_equal(length(full), length(short))
  expect_named(full, c("Parameters", "Genes", "Transcripts", "Abundances"))
  expect_true(all(names(full) == names(short)))
  
  expect_type(full$Parameters, "list")
  expect_true(typeof(full$Parameters) == typeof(short$Parameters))
  expect_length(full$Parameters, 27)
  expect_length(short$Parameters, 2)
  expect_named(full$Parameters, c("description", "time", "rats_version", "R_version",
                                  "var_name", "cond_A", "cond_B", "data_type", "num_replic_A", "num_replic_B", "num_genes", "num_transc",
                                  "tests", "p_thresh", "abund_thresh", "dprop_thresh", "correction", "abund_scaling",
                                  "quant_boot", "quant_reprod_thresh", "quant_bootnum",
                                  "rep_boot", "rep_reprod_thresh", "rep_bootnum", "seed", "reckless", "lean"))
  expect_named(short$Parameters, c("num_replic_A", "num_replic_B"))
  
  expect_true(is.data.frame(full$Genes))
  expect_true(is.data.frame(short$Genes))
  expect_equal(dim(full$Genes)[2], 24)
  expect_equal(dim(short$Genes)[2], 8)
  expect_named(full$Genes, c("parent_id", "elig", "sig", "elig_fx", "quant_reprod", "rep_reprod", "DTU", "transc_DTU",
                             "known_transc", "detect_transc", "elig_transc", "maxDprop", "pval", "pval_corr", 
                             "quant_p_median", "quant_p_min", "quant_p_max", "quant_na_freq", "quant_dtu_freq",
                             "rep_p_median", "rep_p_min", "rep_p_max", "rep_na_freq", "rep_dtu_freq") )
  expect_named(short$Genes, c("parent_id", "elig_transc", "elig", "elig_fx", "pval", "pval_corr", "sig", "DTU") )
  
  expect_true(is.data.frame(full$Transcripts))
  expect_true(is.data.frame(short$Transcripts))
  expect_equal(dim(full$Transcripts)[2], 42)
  expect_equal(dim(short$Transcripts)[2], 17)
  expect_named(full$Transcripts, c("target_id", "parent_id", "elig_xp", "elig", "sig", "elig_fx", "quant_reprod", "rep_reprod", "DTU", "gene_DTU", 
                                   "meanA", "meanB", "stdevA", "stdevB", "sumA", "sumB", "log2FC", "totalA", "totalB", "propA", "propB", "Dprop", "pval", "pval_corr", 
                                   "quant_p_median", "quant_p_min", "quant_p_max", "quant_Dprop_mean", "quant_Dprop_stdev", "quant_Dprop_min", "quant_Dprop_max", "quant_na_freq", "quant_dtu_freq", 
                                   "rep_p_median", "rep_p_min", "rep_p_max", "rep_Dprop_mean", "rep_Dprop_stdev", "rep_Dprop_min", "rep_Dprop_max", "rep_na_freq", "rep_dtu_freq") )
  expect_named(short$Transcripts, c("target_id", "parent_id", "sumA", "sumB", "log2FC", "totalA", "totalB", "elig_xp", "elig",
                                    "propA", "propB", "Dprop", "elig_fx", "pval", "pval_corr", "sig", "DTU") )
  
  expect_type(full$Abundances, "list")
  expect_type(short$Abundances, "NULL")
})


