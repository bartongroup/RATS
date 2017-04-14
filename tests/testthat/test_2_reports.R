#==============================================================================
#==============================================================================
context("DTU Reports")

#==============================================================================
test_that("The summaries work", {
  sim <- sim_sleuth_data(cnames=c("ONE","TWO"))
  mydtu <- call_DTU(annot= sim$annot, slo= sim$slo, name_A= "ONE", name_B= "TWO", qboot=FALSE, rboot=FALSE, verbose = FALSE)
  
  expect_silent(tally <- dtu_summary(mydtu))
  expect_true(is.numeric(tally))
  expect_named(tally, c("DTU genes (gene test)", "non-DTU genes (gene test)", "NA genes (gene test)", 
                        "DTU genes (transc. test)", "non-DTU genes (transc. test)", "NA genes (transc. test)", 
                        "DTU genes (both tests)", "non-DTU genes (both tests)", "NA genes (both tests)", 
                        "DTU transcripts", "non-DTU transcripts", "NA transcripts"))
  expect_false(any(is.na(tally)))
  expect_equal(as.vector(tally),  c(1, 2, 7, 1, 9, 0, 1, 2, 0, 2, 8, 11)) 
  
  expect_silent(ids <- get_dtu_ids(mydtu))
  expect_type(ids, "list")
  expect_named(ids, c("DTU genes (gene test)", "non-DTU genes (gene test)", "NA genes (gene test)", 
                      "DTU genes (transc. test)", "non-DTU genes (transc. test)", "NA genes (transc. test)", 
                      "DTU genes (both tests)", "non-DTU genes (both tests)", "NA genes (both tests)", 
                      "DTU transcripts", "non-DTU transcripts", "NA transcripts"))
  for (v in ids) {
    expect_false(any(is.na(v)))
  }
  expect_equal(as.vector(tally), as.vector(sapply(ids, length)))
  expect_equal(as.vector(unlist(ids)), c("MIX6", "CC", "NN", "LC", "1A1N", "1B1C", "1D1C", "ALLA", "ALLB", "NIB", 
                                         "MIX6", "LC", "CC", "NN", "1A1N", "1B1C", "1D1C", "ALLA", "ALLB", "NIB", 
                                         "MIX6", "CC", "NN", 
                                         "MIX6.c1", "MIX6.c2", "MIX6.c4", "LC2", "CC_a", "CC_b", "MIX6.c3", "2NN", "1NN", "MIX6.nc", "LC1", "1A1N-2", "1B1C.1", "1B1C.2", "1D1C:one", "1D1C:two", "MIX6.d", "ALLA1", "ALLB1", "ALLB2", "NIB.1"))
  
  mydtu <- call_DTU(annot= sim$annot, slo= sim$slo, name_A= "ONE", name_B= "TWO", qboot=FALSE, rboot=FALSE, verbose = FALSE, dprop_thresh=0.1)
  expect_silent(tally <- dtu_summary(mydtu))
  expect_silent(ids <- get_dtu_ids(mydtu))
  expect_equal(as.vector(tally),  c(2, 1, 7, 2, 8, 0, 2, 1, 0, 5, 5, 11)) 
  expect_equal(as.vector(unlist(ids)), c("MIX6", "CC", "NN", "LC", "1A1N", "1B1C", "1D1C", "ALLA", "ALLB", "NIB", 
                                         "MIX6", "CC", "LC", "NN", "1A1N", "1B1C", "1D1C", "ALLA", "ALLB", "NIB", 
                                         "MIX6", "CC", "NN",
                                         "MIX6.c1", "MIX6.c2", "MIX6.c4", "CC_a", "CC_b", "LC2", "MIX6.c3", "2NN", "1NN", "MIX6.nc", "LC1", "1A1N-2", "1B1C.1", "1B1C.2", "1D1C:one", "1D1C:two", "MIX6.d", "ALLA1", "ALLB1", "ALLB2", "NIB.1"))
  expect_equal(as.vector(tally), as.vector(sapply(ids, length)))
})

#==============================================================================
test_that("The isoform switching summaries work", {
  sim <- sim_sleuth_data(cnames=c("ONE","TWO"))
  mydtu <- call_DTU(annot= sim$annot, slo= sim$slo, name_A= "ONE", name_B= "TWO", qboot=FALSE, rboot=FALSE, verbose = FALSE, dprop_thresh=0.1)
  
  expect_silent(tally <- dtu_switch_summary(mydtu))
  expect_true(is.numeric(tally))
  expect_named(tally, c("Primary switch (gene test)", "Non-primary switch (gene test)", 
                        "Primary switch (transc. test)", "Non-primary switch (transc. test)", 
                        "Primary switch (both tests)", "Non-primary switch (both tests)"))
  expect_false(any(is.na(tally)))
  expect_equal(as.vector(tally), c(1, 1, 1, 1, 1, 1))
  
  expect_silent(ids <- get_switch_ids(mydtu))
  expect_type(ids, "list")
  expect_named(ids, c("Primary switch (gene test)", "Non-primary switch (gene test)", 
                      "Primary switch (transc. test)", "Non-primary switch (transc. test)", 
                      "Primary switch (both tests)", "Non-primary switch (both tests)"))
  for (v in ids) {
    expect_false(any(is.na(v)))
  }
  expect_equal(as.vector(tally), as.vector(sapply(ids, length)))
  expect_equal(as.vector(unlist(ids)), c("MIX6", "MIX6", "MIX6", "MIX6", "MIX6", "MIX6"))
 })

#==============================================================================
test_that("The gene plotting commands work", {
  sim <- sim_sleuth_data(cnames=c("ONE","TWO"))
  mydtu <- call_DTU(annot= sim$annot, slo= sim$slo, name_A= "ONE", name_B= "TWO", qbootnum=2, verbose = FALSE)
  
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="lines"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="bycondition"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="merged"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="byisoform"))
  expect_silent(plot_gene(dtuo=mydtu, pid="MIX6", style="linesonly"))
  
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
test_that("The gene plotting commands work", {
  sim <- sim_sleuth_data(cnames=c("ONE","TWO"))
  mydtu <- call_DTU(annot= sim$annot, slo= sim$slo, name_A= "ONE", name_B= "TWO", qbootnum=2, verbose = FALSE)
  
  expect_silent(plot_overview(mydtu))
  expect_silent(plot_overview(dtuo=mydtu, type="volcano"))
  expect_silent(plot_overview(dtuo=mydtu, type="maxdprop"))
})
