#==============================================================================
#==============================================================================
context("DTU Internal structures")

#==============================================================================
test_that("The reporting structures are created correctly", {
  sim <- sim_sleuth_data()
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
  expect_length(full$Parameters, 23)
  expect_length(short$Parameters, 2)
  expect_named(full$Parameters, c("description", "time", "rats_version", "R_version",
                                  "var_name", "cond_A", "cond_B", "data_type", "num_replic_A", "num_replic_B", "num_genes", "num_transc",
                                  "tests", "p_thresh", "abund_thresh", "dprop_thresh",
                                  "quant_reprod_thresh", "quant_boot", "quant_bootnum",
                                  "rep_reprod_thresh", "rep_boot", "rep_bootnum", "rep_reprod_as_crit"))
  expect_named(short$Parameters, c("num_replic_A", "num_replic_B"))
  
  expect_true(is.data.frame(full$Genes))
  expect_true(is.data.frame(short$Genes))
  expect_equal(dim(full$Genes)[2], 35)
  expect_equal(dim(short$Genes)[2], 10)
  expect_named(full$Genes, c("parent_id", "elig", "sig", "elig_fx", "quant_reprod", "rep_reprod", "DTU", "transc_DTU",
                             "known_transc", "detect_transc", "elig_transc", "pvalAB", "pvalBA", "pvalAB_corr", "pvalBA_corr", 
                             "quant_p_meanAB", "quant_p_meanBA", "quant_p_stdevAB", "quant_p_stdevBA", 
                             "quant_p_minAB", "quant_p_minBA", "quant_p_maxAB", "quant_p_maxBA", "quant_na_freq", "quant_dtu_freq",
                             "rep_p_meanAB", "rep_p_meanBA", "rep_p_stdevAB", "rep_p_stdevBA",
                             "rep_p_minAB", "rep_p_minBA", "rep_p_maxAB", "rep_p_maxBA", "rep_na_freq", "rep_dtu_freq") )
  expect_named(short$Genes, c("parent_id", "DTU", "elig_transc", "elig", "elig_fx", "pvalAB", 
                             "pvalBA", "pvalAB_corr", "pvalBA_corr", "sig") )
  
  expect_true(is.data.frame(full$Transcripts))
  expect_true(is.data.frame(short$Transcripts))
  expect_equal(dim(full$Transcripts)[2], 35)
  expect_equal(dim(short$Transcripts)[2], 16)
  expect_named(full$Transcripts, c("target_id", "parent_id", "elig_xp", "elig", "sig", "elig_fx", "quant_reprod", "rep_reprod", "DTU", "gene_DTU", 
                                   "meanA", "meanB", "stdevA", "stdevB", "sumA", "sumB", "totalA", "totalB", "propA", "propB", "Dprop", "pval", "pval_corr", 
                                   "quant_p_mean", "quant_p_stdev", "quant_p_min","quant_p_max", "quant_na_freq", "quant_dtu_freq",
                                   "rep_p_mean", "rep_p_stdev", "rep_p_min","rep_p_max", "rep_na_freq", "rep_dtu_freq") )
  
  expect_type(full$Abundances, "list")
})

