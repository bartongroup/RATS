#==============================================================================
#==============================================================================
context("DTU Reports")

#==============================================================================
test_that("The summaries work", {
  sim <- sim_sleuth_data(cnames=c("ONE","TWO"))
  mydtu <- call_DTU(annot= sim$annot, slo= sim$slo, name_A= "ONE", name_B= "TWO", qbootnum=2, verbose = FALSE)
  
  expect_silent(tally <- dtu_summary(mydtu))
  expect_true(is.numeric(tally))
  expect_named(tally, c("DTU genes (gene test)", "non-DTU genes (gene test)", "NA genes (gene test)", 
                        "DTU genes (transc. test)", "non-DTU genes (transc. test)", "NA genes (transc. test)", 
                        "DTU genes (both tests)", "non-DTU genes (both tests)", "NA genes (both tests)", 
                        "DTU transcripts", "non-DTU transcripts", "NA transcripts"))
  expect_false(any(is.na(tally)))
  
  ids <- get_dtu_ids(mydtu)
  expect_type(ids, "list")
  expect_named(ids, c("DTU genes (gene test)", "non-DTU genes (gene test)", "NA genes (gene test)", 
                      "DTU genes (transc. test)", "non-DTU genes (transc. test)", "NA genes (transc. test)", 
                      "DTU genes (both tests)", "non-DTU genes (both tests)", "NA genes (both tests)", 
                      "DTU transcripts", "non-DTU transcripts", "NA transcripts"))
  for (v in ids) {
    expect_false(any(is.na(v)))
  }
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
