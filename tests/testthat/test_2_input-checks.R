#==============================================================================
#==============================================================================
context("DTU Input checks")

#==============================================================================
test_that("The input checks work", {
  name_A <- "one"
  name_B <- "two"
  wrong_name <- "RUBBISH_COLUMN_NAME"
  
  # Emulate bootstrap data.
  sim1 <- sim_boot_data()
  sim2 <- sim_boot_data(varname="waffles", cnames=c("AAAA", "BBBB"), errannot_inconsistent=FALSE, PARENT_COL="parent", TARGET_COL="target")
  # Emulate non-bootstrap data.
  counts_A <- sim1$boots_A[[1]]
  counts_B <- sim1$boots_B[[1]]
  
  # No false alarms with valid calls.
  expect_silent(call_DTU(annot= sim2$annot, boot_data_A= sim2$boots_A, boot_data_B= sim2$boots_B, name_A = "AAAA", name_B = "BBBB", varname= "waffles", p_thresh= 0.01, abund_thresh= 10,
                         rrep_thresh = 0.6, qrep_thresh = 0.8, testmode= "transc", correction= "bonferroni", verbose= FALSE, rboot= FALSE, qboot=TRUE,
                         qbootnum= 100, TARGET_COL= "target", PARENT_COL= "parent", threads= 2, dbg= "prep"))
  expect_silent(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, verbose = FALSE, dbg= "prep"))
  expect_silent(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, rboot=FALSE, qboot= FALSE, verbose = FALSE, dbg= "prep"))
  
  # Mixed input formats.
  expect_silent(call_DTU(annot= sim1$annot, name_A= name_A, name_B= name_B, verbose = FALSE,
                         boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, count_data_A= counts_A, count_data_B= counts_B, dbg= "prep"))
  expect_warning(call_DTU(annot= sim1$annot, verbose = TRUE, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, 
                          count_data_A= counts_A, count_data_B= counts_B, qbootnum= 100, dbg= "prep"),
                 "Only the bootstrapped data will be used", fixed= TRUE)
  
  # Annotation is not a dataframe.
  expect_error(call_DTU(annot= list("not", "a", "dataframe"), boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, verbose = FALSE, dbg= "prep"), "annot is not a data.frame")
  # Annotation field names.
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, TARGET_COL= wrong_name, verbose = FALSE, dbg= "prep"),
               "target and/or parent IDs field names do not exist in annot", fixed= TRUE)
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, PARENT_COL= wrong_name, verbose = FALSE, dbg= "prep"),
               "target and/or parent IDs field names do not exist in annot", fixed= TRUE)
  # Inconsistent annotation.
  sim3 <- sim_boot_data(errannot_inconsistent= TRUE, cnames= c(name_A, name_B))
  expect_error(call_DTU(annot= sim3$annot, boot_data_A= sim3$boots_A, boot_data_B= sim3$boots_B, name_A= name_A, name_B= name_B, verbose = FALSE, dbg= "prep"),
               "Inconsistent set of transcript IDs", fixed= TRUE)
  # Non unique IDs.
  a <- copy(sim1$annot)
  a[1, "target_id"] <- a[2, "target_id"]
  expect_error(call_DTU(annot= a, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, verbose = FALSE, dbg= "prep"), "transcript identifiers are not unique")
  
  # Boot data is not a list of datatables.
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= c("not", "a", "list"), boot_data_B= sim1$boots_B, verbose = FALSE, dbg= "prep"), "bootstrap data are not lists")
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= c("not", "a", "list"), verbose = FALSE, dbg= "prep"), "bootstrap data are not lists")
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= list( list("not"), list("a"), list("table")), boot_data_B= sim1$boots_B, verbose = FALSE, dbg= "prep"), "bootstrap data are not lists")
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= list( list("not"), list("a"), list("table")), verbose = FALSE, dbg= "prep"), "bootstrap data are not lists")
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= list(data.frame(a="not", b="a", c="table")), boot_data_B= sim1$boots_B, verbose = FALSE, dbg= "prep"), "not lists of data.tables")
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= list(data.frame(a="not", b="a", c="table")), verbose = FALSE, dbg= "prep"), "not lists of data.tables")
  
  # Counts data is not datatables.
  expect_error(call_DTU(annot= sim1$annot, count_data_A= c("not", "a", "list"), count_data_B= counts_B, verbose = FALSE, dbg= "prep"), "counts data are not data.tables")
  expect_error(call_DTU(annot= sim1$annot, count_data_A= counts_A, count_data_B= c("not", "a", "list"), verbose = FALSE, dbg= "prep"), "counts data are not data.tables")
  expect_error(call_DTU(annot= sim1$annot, count_data_A= data.frame(a="not", b="a", c="table"), count_data_B= counts_B, verbose = FALSE, dbg= "prep"), "counts data are not data.tables")
  expect_error(call_DTU(annot= sim1$annot, count_data_A= counts_A, count_data_B= data.frame(a="not", b="a", c="table"), verbose = FALSE, dbg= "prep"), "counts data are not data.tables")
  
  # Number of bootstraps.
  expect_warning(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, verbose = TRUE, dbg= "prep"),
                 "few bootstrap iterations", fixed= TRUE)
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, qbootnum = -5, verbose = FALSE, dbg= "prep"),
               "Invalid number of bootstraps", fixed= TRUE)
  expect_warning(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, qbootnum = 5, verbose = TRUE, dbg= "prep"),
                 "qbootnum is low", fixed= TRUE)
  expect_warning(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, qbootnum = 5000, verbose = TRUE, dbg= "prep"),
                 "number of quantification bootstraps is very high", fixed= TRUE)
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, qbootnum = 5000000000000, verbose = FALSE, dbg= "prep"),
               "number of quantification bootstraps would exceed the maximum capacity", fixed= TRUE)

  # Bootstraps without boot data.
  expect_warning(call_DTU(annot= sim1$annot, count_data_A= counts_A, count_data_B= counts_B, qboot=TRUE, qbootnum=2, verbose= TRUE, dbg= "prep"), 
               "qboot is TRUE but no bootstrapped estimates", fixed= TRUE)
  
  # Correction method.
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, correction= wrong_name, verbose= FALSE, dbg= "prep"),
               "Invalid p-value correction method name", fixed= TRUE)
  
  # Verbose is bool.
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, verbose="yes", dbg= "prep"),
               "not interpretable as logical", fixed= TRUE)
  
  # Tests.
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, testmode="GCSE", verbose = FALSE, dbg= "prep"),
               "Unrecognized value for testmode", fixed= TRUE)
  expect_silent(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, testmode="genes", verbose = FALSE, rboot=FALSE, qboot = FALSE, dbg= "prep"))
  expect_silent(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, testmode="transc", verbose = FALSE, rboot=FALSE, qboot = FALSE, dbg= "prep"))
  
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, qboot="GCSE", verbose = FALSE, dbg= "prep"),
               "Unrecognized value for qboot", fixed= TRUE)
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, rboot="GCSE", verbose = FALSE, dbg= "prep"),
               "Unrecognized value for rboot", fixed= TRUE)
  
  # Probability threshold.
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, p_thresh = 666, verbose = FALSE, dbg= "prep"),
               "Invalid p-value threshold", fixed= TRUE)
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, p_thresh = -0.05, verbose = FALSE, dbg= "prep"),
               "Invalid p-value threshold", fixed= TRUE)
  
  # Read counts threshold.
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, abund_thresh = -5, verbose = FALSE, dbg= "prep"),
               "Invalid abundance threshold", fixed= TRUE)
  
  # Proportion change threshold.
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, dprop_thresh = -2, verbose = FALSE, dbg= "prep"),
               "Invalid proportion difference threshold", fixed= TRUE)
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, dprop_thresh = 2, verbose = FALSE, dbg= "prep"),
               "Invalid proportion difference threshold", fixed= TRUE)
  
  # Confidence threshold.
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, rrep_thresh = -2, verbose = FALSE, dbg= "prep"),
               "Invalid reproducibility threshold", fixed= TRUE)
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, rrep_thresh = 2, verbose = FALSE, dbg= "prep"),
               "Invalid reproducibility threshold", fixed= TRUE)
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, qrep_thresh = -2, verbose = FALSE, dbg= "prep"),
               "Invalid reproducibility threshold", fixed= TRUE)
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, qrep_thresh = 2, verbose = FALSE, dbg= "prep"),
               "Invalid reproducibility threshold", fixed= TRUE)
  
  # Number of threads.
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, threads = 0, verbose= FALSE, dbg= "prep"),
               "Invalid number of threads", fixed= TRUE)
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, threads = -4, verbose= FALSE, dbg= "prep"),
               "Invalid number of threads", fixed= TRUE)
  expect_silent(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, threads = parallel::detectCores(logical=FALSE), verbose= FALSE, dbg= "prep"))
  expect_silent(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, threads = parallel::detectCores(logical=TRUE), verbose= FALSE, dbg= "prep"))
  expect_error(call_DTU(annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, name_A= name_A, name_B= name_B, threads = parallel::detectCores(logical=TRUE) + 1, verbose= FALSE, dbg= "prep"),
               "threads exceed", fixed= TRUE)
})


