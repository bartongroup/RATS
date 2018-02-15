#==============================================================================
#==============================================================================
context("DTU Reports")

#==============================================================================
test_that("The summaries work", {
  sim <- sim_boot_data()
  mydtu <- call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= "ONE", name_B= "TWO", qboot=FALSE, rboot=FALSE, verbose = FALSE)
  
  expect_silent(ids <- get_dtu_ids(mydtu))
  expect_type(ids, "list")
  expect_named(ids, c("DTU genes (gene test)", "non-DTU genes (gene test)", "ineligible genes (gene test)", 
                      "DTU genes (transc. test)", "non-DTU genes (transc. test)", "ineligible genes (transc. test)", 
                      "DTU genes (both tests)", "non-DTU genes (both tests)", "ineligible genes (both tests)", 
                      "DTU transcripts", "non-DTU transcripts", "ineligible transcripts"))
  for (v in ids) {
    expect_false(any(is.na(v)))
  }
  expect_equal(ids[[1]], c("MIX6"))
  expect_equal(ids[[2]], c("CC", "NN"))
  expect_equal(ids[[3]], c("LC", "1A1N", "1B1C",  "1D1C",  "ALLA", "ALLB", "NIB"))
  expect_equal(ids[[4]], c("MIX6"))
  expect_equal(ids[[5]], c("LC", "CC", "NN"))
  expect_equal(ids[[6]], c("1A1N", "1B1C", "1D1C", "ALLA", "ALLB", "NIB"))
  expect_equal(ids[[7]], c("MIX6"))
  expect_equal(ids[[8]], c("CC", "NN"))
  expect_equal(ids[[9]], c("LC", "1A1N", "1B1C", "1D1C", "ALLA", "ALLB", "NIB"))
  expect_equal(ids[[10]], c("MIX6.c1", "MIX6.c2"))
  expect_equal(ids[[11]], c("MIX6.c4", "LC2", "CC_a", "CC_b", "MIX6.c3", "2NN", "1NN", "MIX6.nc"))
  expect_equal(ids[[12]], c("LC1", "1A1N-2", "1B1C.1", "1B1C.2", "1D1C:one", "1D1C:two", "MIX6.d", "ALLA1", "ALLB1", "ALLB2", "NIB.1"))
  
  expect_silent(tally <- dtu_summary(mydtu))
  expect_true(is.data.frame(tally))
  expect_equal(tally[[1]], c("DTU genes (gene test)", "non-DTU genes (gene test)", "ineligible genes (gene test)", 
                        "DTU genes (transc. test)", "non-DTU genes (transc. test)", "ineligible genes (transc. test)", 
                        "DTU genes (both tests)", "non-DTU genes (both tests)", "ineligible genes (both tests)", 
                        "DTU transcripts", "non-DTU transcripts", "ineligible transcripts"))
  expect_equal(tally[[2]], as.vector(sapply(ids, length)))
  expect_equal(tally[[2]],  c(1, 2, 7, 1, 3, 6, 1, 2, 7, 2, 8, 11)) 
  expect_false(any(is.na(tally)))
  
  # Lower effect size threshold to verify that DTU changes accordingly.
  mydtu <- call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= "ONE", name_B= "TWO", qboot=FALSE, rboot=FALSE, verbose = FALSE, dprop_thresh=0.1)
  
  expect_silent(ids <- get_dtu_ids(mydtu))
  expect_equal(ids[[1]], c("MIX6", "CC"))
  expect_equal(ids[[2]], c("NN"))
  expect_equal(ids[[3]], c("LC", "1A1N", "1B1C", "1D1C", "ALLA", "ALLB", "NIB"))
  expect_equal(ids[[4]], c("MIX6", "CC"))
  expect_equal(ids[[5]], c("LC", "NN"))
  expect_equal(ids[[6]], c("1A1N", "1B1C", "1D1C", "ALLA", "ALLB", "NIB"))
  expect_equal(ids[[7]], c("MIX6", "CC"))
  expect_equal(ids[[8]], c("NN"))
  expect_equal(ids[[9]], c("LC", "1A1N", "1B1C", "1D1C", "ALLA", "ALLB", "NIB"))
  expect_equal(ids[[10]], c("MIX6.c1",  "MIX6.c2", "MIX6.c4", "CC_a", "CC_b"))
  expect_equal(ids[[11]], c("LC2", "MIX6.c3", "2NN", "1NN", "MIX6.nc"))
  expect_equal(ids[[12]], c("LC1", "1A1N-2", "1B1C.1", "1B1C.2", "1D1C:one", "1D1C:two", "MIX6.d", "ALLA1", "ALLB1", "ALLB2", "NIB.1"))
  
  expect_silent(tally <- dtu_summary(mydtu))
  expect_equal(tally[[2]], as.vector(sapply(ids, length)))
  expect_equal(tally[[2]],  c(2, 1, 7, 2, 2, 6, 2, 1, 7, 5, 5, 11)) 
})

#==============================================================================
test_that("The isoform switching summaries work", {
  sim <- sim_boot_data()
  mydtu <- call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= "ONE", name_B= "TWO", qboot=FALSE, rboot=FALSE, verbose = FALSE, dprop_thresh=0.1)
  
  expect_silent(ids <- get_switch_ids(mydtu))
  expect_type(ids, "list")
  expect_named(ids, c("Primary switch (gene test)", "Non-primary switch (gene test)", 
                      "Primary switch (transc. test)", "Non-primary switch (transc. test)", 
                      "Primary switch (both tests)", "Non-primary switch (both tests)"))
  for (v in ids) {
    expect_false(any(is.na(v)))
  }
  expect_equal(as.vector(unlist(ids)), c("MIX6", "MIX6", "MIX6", "MIX6", "MIX6", "MIX6"))
  
  expect_silent(tally <- dtu_switch_summary(mydtu))
  expect_true(is.data.frame(tally))
  expect_equal(tally[[1]], c("Primary switch (gene test)", "Non-primary switch (gene test)", 
                        "Primary switch (transc. test)", "Non-primary switch (transc. test)", 
                        "Primary switch (both tests)", "Non-primary switch (both tests)"))
  expect_false(any(is.na(tally)))
  expect_equal(tally[[2]], c(1, 1, 1, 1, 1, 1))
  expect_equal(tally[[2]], as.vector(sapply(ids, length)))
 })


#==============================================================================
test_that("The plurality summaries work", {
  sim <- sim_boot_data()
  mydtu <- call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= "ONE", name_B= "TWO", qboot=FALSE, rboot=FALSE, verbose = FALSE, dprop_thresh=0.1)
  
  expect_silent(ids <- get_plurality_ids(mydtu))
  expect_type(ids, "list")
  expect_named(ids, c("2", "3"))
  for (v in ids) {
    expect_false(any(is.na(v)))
  }
  expect_equal(as.vector(unlist(ids)), c("CC", "MIX6"))
  
  expect_silent(tally <- dtu_plurality_summary(mydtu))
  expect_true(is.data.frame(tally))
  expect_equal(tally[[1]], c("2", "3"))
  expect_false(any(is.na(tally)))
  expect_equal(tally[[2]], c(1, 1))
  expect_equal(tally[[2]], as.vector(sapply(ids, length)))
  
  mydtu <- call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= "ONE", name_B= "TWO", qboot=FALSE, rboot=FALSE, verbose = FALSE, dprop_thresh=0.3)
  expect_silent(ids <- get_plurality_ids(mydtu))
  expect_silent(tally <- dtu_plurality_summary(mydtu))
  expect_named(ids, c("1"))
  expect_equal(as.vector(unlist(ids)), c("MIX6"))
  expect_equal(tally[[2]], c(1))
})


#==============================================================================
test_that("The gene plotting commands work", {
  sim <- sim_boot_data()
  mydtu <- call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= "ONE", name_B= "TWO", qbootnum=2, verbose = FALSE)
  
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="bycondition"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="byisoform"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="lines"))
  
  expect_error(plot_gene(dtuo=mydtu, pid="MIX6", fillby="replicate"))
  expect_error(plot_gene(dtuo=mydtu, pid="MIX6", style="bycondition", fillby="garbage"))
  expect_error(plot_gene(dtuo=mydtu, pid="MIX6", style="bycondition", colourby="garbage"))
  expect_error(plot_gene(dtuo=mydtu, pid="MIX6", style="bycondition", shapeby="garbage"))
  expect_error(plot_gene(dtuo=mydtu, pid="MIX6", style="byisoform", fillby="none", colourby="none"))
  
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="bycondition", fillby="none"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="bycondition", fillby="isoform"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="bycondition", fillby="condition"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="bycondition", fillby="DTU"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="bycondition", colourby="replicate"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="bycondition", shapeby="replicate"))
  
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="byisoform", fillby="isoform"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="byisoform", fillby="DTU"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="byisoform", fillby="none"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="byisoform", fillby="condition"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="byisoform", colourby="isoform"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="byisoform", colourby="condition"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="byisoform", colourby="DTU"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="byisoform", colourby="replicate"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="byisoform", colourby="none"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="byisoform", shapeby="isoform"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="byisoform", shapeby="condition"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="byisoform", shapeby="DTU"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="byisoform", shapeby="replicate"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="byisoform", shapeby="none"))
})

#==============================================================================
test_that("The overview plotting commands work", {
  sim <- sim_boot_data()
  mydtu <- call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= "ONE", name_B= "TWO", qbootnum=2, verbose = FALSE)
  
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
