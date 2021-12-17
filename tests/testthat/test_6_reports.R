#==============================================================================
#==============================================================================
context("DTU Reports")

#==============================================================================
test_that("The summaries work", {
  sim <- sim_boot_data()
  mydtu <- call_DTU(annot= sim$annot, boot_data_A=sim$boots_A, boot_data_B=sim$boots_B, name_A= "ONE", name_B= "TWO", abund_thresh=20, dprop_thresh=0.2, verbose = FALSE, rboot=FALSE, qboot = FALSE, seed=666, reckless=TRUE, use_sums=TRUE)
  
  expect_silent(ids <- get_dtu_ids(mydtu))
  expect_type(ids, "list")
  expect_named(ids, c("DTU genes (gene test)", "non-DTU genes (gene test)", "ineligible genes (gene test)", 
                      "DTU genes (transc. test)", "non-DTU genes (transc. test)", "ineligible genes (transc. test)", 
                      "DTU genes (both tests)", "non-DTU genes (both tests)", "ineligible genes (both tests)", 
                      "DTU transcripts", "non-DTU transcripts", "ineligible transcripts"))
  for (v in ids) {
    expect_false(any(is.na(v)))
  }
  expect_equal(ids[[1]], c("XSW", "MIX", "SW"))
  expect_equal(ids[[2]], c("D2TU", "D2TE", "SAME"))
  expect_equal(ids[[3]], c("LC", "ALLA", "FAKE", "LONE",  "SOLO",  "NN"))
  expect_equal(ids[[4]], c("XSW", "MIX", "SW"))
  expect_equal(ids[[5]], c("D2TU", "D2TE", "SAME", "FAKE"))
  expect_equal(ids[[6]], c("LC", "ALLA", "LONE", "SOLO", "NN"))
  expect_equal(ids[[7]], c("XSW", "MIX", "SW"))
  expect_equal(ids[[8]], c("D2TU", "D2TE", "SAME"))
  expect_equal(ids[[9]], c("LC", "ALLA", "FAKE", "LONE", "SOLO", "NN"))
  expect_equal(ids[[10]], c("XSW:one", "XSW:two", "MIX.b", "SW1", "SW2", "MIX.a"))
  expect_equal(ids[[11]], c("MIX.ab", "2D2TU", "1D2TU", "MIX.l2", "D2TE_b", "D2TE_a", "SAME_1", "SAME_2", "FAKE-2"))
  expect_equal(ids[[12]], c("LC1", "LC2", "ALLA:1", "ALLA:2", "FAKE-1", "MIX.l1", "LONE.a", "MIX.n", "SOLO.1", "NNa", "NNb"))
  
  expect_silent(tally <- dtu_summary(mydtu))
  expect_true(is.data.frame(tally))
  expect_equal(tally[[1]], c("DTU genes (gene test)", "non-DTU genes (gene test)", "ineligible genes (gene test)", 
                        "DTU genes (transc. test)", "non-DTU genes (transc. test)", "ineligible genes (transc. test)", 
                        "DTU genes (both tests)", "non-DTU genes (both tests)", "ineligible genes (both tests)", 
                        "DTU transcripts", "non-DTU transcripts", "ineligible transcripts"))
  expect_equivalent(tally[[2]], unlist(lapply(ids, length)))
  expect_equal(tally[[2]],  c(3, 3, 6, 3, 4, 5, 3, 3, 6, 6, 9, 11)) 
  expect_false(any(is.na(tally)))
  
  # Lower threshold to verify that DTU changes accordingly.
  mydtu <- call_DTU(annot= sim$annot, boot_data_A=sim$boots_A, boot_data_B=sim$boots_B, name_A= "ONE", name_B= "TWO", abund_thresh=10, dprop_thresh=0.1, verbose = FALSE, rboot=FALSE, qboot = FALSE, seed=666, reckless=TRUE, use_sums=TRUE)
  
  ids <- get_dtu_ids(mydtu)
  expect_equal(ids[[1]], c("XSW", "MIX", "SW", "LC", "ALLA", "D2TU"))
  expect_equal(ids[[2]], c("D2TE", "SAME"))
  expect_equal(ids[[3]], c("FAKE", "LONE",  "SOLO",  "NN"))
  expect_equal(ids[[4]], c("XSW", "MIX", "SW", "LC", "ALLA", "D2TU"))
  expect_equal(ids[[5]], c("D2TE", "SAME", "FAKE"))
  expect_equal(ids[[6]], c("LONE", "SOLO", "NN"))
  expect_equal(ids[[7]], c("XSW", "MIX", "SW", "LC", "ALLA", "D2TU"))
  expect_equal(ids[[8]], c("D2TE", "SAME"))
  expect_equal(ids[[9]], c("FAKE", "LONE", "SOLO", "NN"))
  expect_equal(ids[[10]], c("XSW:one", "XSW:two", "MIX.b", "SW1", "SW2", "MIX.a", "LC1", "LC2", "ALLA:1", "ALLA:2", "MIX.ab", "2D2TU", "1D2TU", "MIX.l2"))
  expect_equal(ids[[11]], c("D2TE_b", "D2TE_a", "SAME_1", "SAME_2", "FAKE-2", "MIX.l1"))
  expect_equal(ids[[12]], c("FAKE-1", "LONE.a", "MIX.n", "SOLO.1", "NNa", "NNb"))
  
  expect_silent(tally <- dtu_summary(mydtu))
  expect_equivalent(tally[[2]], unlist(lapply(ids, length)))
  expect_equal(tally[[2]],  c(6, 2, 4, 6, 3, 3, 6, 2, 4, 14, 6, 6)) 
})

#==============================================================================
test_that("The isoform switching summaries work", {
  sim <- sim_boot_data()
  mydtu <- call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= "ONE", name_B= "TWO", qboot=FALSE, rboot=FALSE, verbose = FALSE, dprop_thresh=0.1, reckless=TRUE, use_sums=TRUE)
  
  expect_silent(ids <- get_switch_ids(mydtu))
  expect_type(ids, "list")
  expect_named(ids, c("Primary switch (gene test)", "Non-primary switch (gene test)", 
                      "Primary switch (transc. test)", "Non-primary switch (transc. test)", 
                      "Primary switch (both tests)", "Non-primary switch (both tests)"))
  for (v in ids) {
    expect_false(any(is.na(v)))
  }
  expect_equivalent(unlist(ids), c("LC", "MIX", "SW", "XSW", "MIX", "LC", "MIX", "SW", "XSW", "MIX", "LC", "MIX", "SW", "XSW", "MIX"))
  
  expect_silent(tally <- dtu_switch_summary(mydtu))
  expect_true(is.data.frame(tally))
  expect_equal(tally[[1]], c("Primary switch (gene test)", "Non-primary switch (gene test)", 
                        "Primary switch (transc. test)", "Non-primary switch (transc. test)", 
                        "Primary switch (both tests)", "Non-primary switch (both tests)"))
  expect_false(any(is.na(tally)))
  expect_equal(tally[[2]], c(4, 1, 4, 1, 4, 1))
  expect_equivalent(tally[[2]], unlist(lapply(ids, length)))
 })


#==============================================================================
test_that("The plurality summaries work", {
  sim <- sim_boot_data()
  mydtu <- call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= "ONE", name_B= "TWO", qboot=FALSE, rboot=FALSE, verbose = FALSE, dprop_thresh=0.1, reckless=TRUE, use_sums=TRUE)
  
  expect_silent(ids <- get_plurality_ids(mydtu))
  expect_type(ids, "list")
  expect_named(ids, c("2", "4"))
  for (v in ids) {
    expect_false(any(is.na(v)))
  }
  expect_equivalent(unlist(ids), c("ALLA", "D2TU", "LC", "SW", "XSW", "MIX"))
  
  expect_silent(tally <- dtu_plurality_summary(mydtu))
  expect_true(is.data.frame(tally))
  expect_equal(tally[[1]], c("2", "4"))
  expect_false(any(is.na(tally)))
  expect_equal(tally[[2]], c(5, 1))
  expect_equivalent(tally[[2]], unlist(lapply(ids, length)))
  
  mydtu <- call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= "ONE", name_B= "TWO", qboot=FALSE, rboot=FALSE, verbose = FALSE, dprop_thresh=0.3, reckless=TRUE, use_sums=TRUE)
  expect_silent(ids <- get_plurality_ids(mydtu))
  expect_silent(tally <- dtu_plurality_summary(mydtu))
  expect_named(ids, c("2"))
  expect_equivalent(unlist(ids), c("LC", "MIX", "SW", "XSW"))
  expect_equal(tally[[2]], c(4))
})


#==============================================================================
test_that("The gene plotting commands work", {
  sim <- sim_boot_data()
  mydtu <- call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= "ONE", name_B= "TWO", qbootnum=2, verbose = FALSE, reckless=TRUE)
  
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX", style="bycondition"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX", style="byisoform"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX", style="lines"))
  
  expect_error(plot_gene(dtuo=mydtu, pid="MIX", fillby="replicate"))
  expect_error(plot_gene(dtuo=mydtu, pid="MIX", style="bycondition", fillby="garbage"))
  expect_error(plot_gene(dtuo=mydtu, pid="MIX", style="bycondition", colourby="garbage"))
  expect_error(plot_gene(dtuo=mydtu, pid="MIX", style="bycondition", shapeby="garbage"))
  expect_error(plot_gene(dtuo=mydtu, pid="MIX", style="byisoform", fillby="none", colourby="none"))
  
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX", style="bycondition", fillby="none"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX", style="bycondition", fillby="isoform"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX", style="bycondition", fillby="condition"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX", style="bycondition", fillby="DTU"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX", style="bycondition", colourby="replicate"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX", style="bycondition", shapeby="replicate"))
  
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX", style="byisoform", fillby="isoform"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX", style="byisoform", fillby="DTU"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX", style="byisoform", fillby="none"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX", style="byisoform", fillby="condition"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX", style="byisoform", colourby="isoform"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX", style="byisoform", colourby="condition"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX", style="byisoform", colourby="DTU"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX", style="byisoform", colourby="replicate"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX", style="byisoform", colourby="none"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX", style="byisoform", shapeby="isoform"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX", style="byisoform", shapeby="condition"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX", style="byisoform", shapeby="DTU"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX", style="byisoform", shapeby="replicate"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX", style="byisoform", shapeby="none"))
})

#==============================================================================
test_that("All gene scenarios can be plotted", {
  sim <- sim_boot_data(clean=FALSE)
  mydtu <- call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, qbootnum=2, verbose = FALSE, reckless=TRUE)
  
  expect_silent(plot_gene(dtuo=mydtu, pid="ALLA"))
  expect_silent(plot_gene(dtuo=mydtu, pid="D2TE"))
  expect_silent(plot_gene(dtuo=mydtu, pid="D2TU"))
  expect_silent(plot_gene(dtuo=mydtu, pid="FAKE"))
  expect_silent(plot_gene(dtuo=mydtu, pid="LC"))
  expect_silent(plot_gene(dtuo=mydtu, pid="LONE"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX"))
  expect_silent(plot_gene(dtuo=mydtu, pid="NID")) #
  expect_silent(plot_gene(dtuo=mydtu, pid="NN"))  #
  expect_silent(plot_gene(dtuo=mydtu, pid="SAME"))
  expect_silent(plot_gene(dtuo=mydtu, pid="SOLO"))
  expect_silent(plot_gene(dtuo=mydtu, pid="SW"))
  expect_silent(plot_gene(dtuo=mydtu, pid="XSW"))
})

#==============================================================================
test_that("The overview plotting commands work", {
  sim <- sim_boot_data()
  mydtu <- call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= "ONE", name_B= "TWO", qbootnum=2, verbose = FALSE, reckless=TRUE)
  
  expect_silent(plot_overview(mydtu))
  expect_silent(plot_overview(dtuo=mydtu, type="tvolcano"))
  expect_silent(plot_overview(dtuo=mydtu, type="gvolcano"))
  expect_silent(plot_overview(dtuo=mydtu, type="fcvolcano"))
  expect_silent(plot_overview(dtuo=mydtu, type="dprop"))
  expect_silent(plot_overview(dtuo=mydtu, type="maxdprop"))
  expect_silent(plot_overview(dtuo=mydtu, type="fcVSdprop"))
  expect_silent(plot_overview(dtuo=mydtu, type="reprod"))
  expect_silent(plot_overview(dtuo=mydtu, type="reprodVSdprop"))
})

#==============================================================================
test_that("The diagnostic plots commands work", {
  sim <- sim_boot_data()
  mydtu <- call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= "ONE", name_B= "TWO", qbootnum=2, verbose = FALSE, reckless=TRUE)
  
  expect_silent(plot_diagnostics(mydtu))
  expect_silent(plot_diagnostics(dtuo=mydtu, type="cormat"))
})
