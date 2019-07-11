#==============================================================================
#==============================================================================
context("DTU Input checks")

#==============================================================================
test_that("No false alarms", {
  name_A <- "one"
  name_B <- "two"
  
  sim <- sim_boot_data(clean=TRUE, errannot_inconsistent=FALSE, TARGET_COL= "target", PARENT_COL= "parent")
  expect_silent(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A = "AAAA", name_B = "BBBB", varname= "waffles", p_thresh= 0.01, abund_thresh= 10,
                         rrep_thresh = 0.6, qrep_thresh = 0.8, testmode= "transc", correction= "bonferroni", verbose= FALSE, rboot= FALSE, qboot=TRUE,
                         qbootnum= 100, TARGET_COL= "target", PARENT_COL= "parent", threads= 2, dbg= "prep", scaling=c(10, 11, 20, 21), lean=FALSE, use_sums=TRUE))
  
  sim <- sim_boot_data(clean=TRUE, errannot_inconsistent=FALSE)
  expect_silent(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, verbose = FALSE, dbg= "prep", scaling=30))
  expect_silent(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, rboot=FALSE, qboot= FALSE, verbose = FALSE, dbg= "prep"))
  expect_silent(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, verbose = FALSE, dbg= "prep"))
  expect_silent(call_DTU(annot= sim$annot, count_data_A= sim$boots_A[[1]], count_data_B= sim$boots_B[[1]], verbose = FALSE, dbg= "prep"))
})
  
#==============================================================================
test_that("Feature ID inconsistencies abort or warn", {
  sim <- sim_boot_data(clean=TRUE, errannot_inconsistent=TRUE)
  
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, verbose = FALSE, dbg= "prep"),
                 "Inconsistent set of transcript IDs across samples", fixed= TRUE)
  expect_error(call_DTU(annot= sim$annot, count_data_A= sim$boots_A[[1]], count_data_B= sim$boots_B[[1]], verbose = FALSE, dbg= "prep"),
                 "Transcript IDs do not match completely between conditions", fixed= TRUE)
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_B, boot_data_B= sim$boots_A, verbose = FALSE, dbg= "prep"),
               "The transcript IDs in the quantifications and the annotation do not match", fixed= TRUE)
  expect_error(call_DTU(annot= sim$annot, count_data_A= sim$boots_B[[1]], count_data_B= sim$boots_A[[1]], verbose = FALSE, dbg= "prep"),
               "The transcript IDs in the quantifications and the annotation do not match", fixed= TRUE)
  
  expect_warning(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, verbose = TRUE, dbg= "prep", reckless=TRUE),
               "Inconsistent set of transcript IDs across samples", fixed= TRUE)
  expect_warning(call_DTU(annot= sim$annot, count_data_A= sim$boots_A[[1]], count_data_B= sim$boots_B[[1]], verbose = TRUE, dbg= "prep", reckless=TRUE),
               "Transcript IDs do not match completely between conditions", fixed= TRUE)
  expect_warning(call_DTU(annot= sim$annot, boot_data_A= sim$boots_B, boot_data_B= sim$boots_A, verbose = TRUE, dbg= "prep", reckless=TRUE),
               "The transcript IDs in the quantifications and the annotation do not match", fixed= TRUE)
  expect_warning(call_DTU(annot= sim$annot, count_data_A= sim$boots_B[[1]], count_data_B= sim$boots_A[[1]], verbose = TRUE, dbg= "prep", reckless=TRUE),
               "The transcript IDs in the quantifications and the annotation do not match", fixed= TRUE)
  
  expect_silent(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, verbose = FALSE, dbg= "prep", reckless=TRUE))
  expect_silent(call_DTU(annot= sim$annot, count_data_A= sim$boots_A[[1]], count_data_B= sim$boots_B[[1]], verbose = FALSE, dbg= "prep", reckless=TRUE))
  expect_silent(call_DTU(annot= sim$annot, boot_data_A= sim$boots_B, boot_data_B= sim$boots_A, verbose = FALSE, dbg= "prep", reckless=TRUE))
  expect_silent(call_DTU(annot= sim$annot, count_data_A= sim$boots_B[[1]], count_data_B= sim$boots_A[[1]], verbose = FALSE, dbg= "prep", reckless=TRUE))
})
  
#==============================================================================
test_that("Annotation format is checked", {
  sim <- sim_boot_data(clean=TRUE)
  
  # Annotation is not a dataframe.
  expect_error(call_DTU(annot= list("not", "a", "dataframe"), boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, verbose = FALSE, dbg= "prep"), 
               "annot is not a data.frame")
  # Annotation field names.
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, TARGET_COL= "wrong_name", verbose = FALSE, dbg= "prep"),
               "target and/or parent IDs field names do not exist in annot", fixed= TRUE)
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, PARENT_COL= "wrong_name", verbose = FALSE, dbg= "prep"),
               "target and/or parent IDs field names do not exist in annot", fixed= TRUE)
  
  # Non unique transcript IDs.
  a <- copy(sim$annot)
  a[1, "target_id"] <- a[2, "target_id"]
  expect_error(call_DTU(annot= a, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, verbose = FALSE, dbg= "prep"), 
               "transcript identifiers are not unique")
})
  
#==============================================================================
test_that("Quantification formats are checked", {
  sim <- sim_boot_data(clean=TRUE)
  
  # Redundant quantification formats.
  expect_silent(call_DTU(annot= sim$annot, verbose = FALSE, dbg= "prep", boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, count_data_A= sim$boots_A[[1]], count_data_B= sim$boots_B[[1]]))
  expect_warning(call_DTU(annot= sim$annot, verbose = TRUE, dbg= "prep", boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, count_data_A= sim$boots_A[[1]], count_data_B= sim$boots_B[[1]]),
                 "Only the bootstrapped data will be used", fixed= TRUE)
  
  # Boot data is not a list of datatables.
  expect_error(call_DTU(annot= sim$annot, boot_data_A= c("not", "a", "list"), boot_data_B= sim$boots_B, verbose = FALSE, dbg= "prep"), "bootstrap data are not lists")
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= c("not", "a", "list"), verbose = FALSE, dbg= "prep"), "bootstrap data are not lists")
  expect_error(call_DTU(annot= sim$annot, boot_data_A= list( list("not"), list("a"), list("table")), boot_data_B= sim$boots_B, verbose = FALSE, dbg= "prep"), "bootstrap data are not lists")
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= list( list("not"), list("a"), list("table")), verbose = FALSE, dbg= "prep"), "bootstrap data are not lists")
  expect_error(call_DTU(annot= sim$annot, boot_data_A= list(data.frame(a="not", b="a", c="table")), boot_data_B= sim$boots_B, verbose = FALSE, dbg= "prep"), "not lists of data.tables")
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= list(data.frame(a="not", b="a", c="table")), verbose = FALSE, dbg= "prep"), "not lists of data.tables")
  
  # Counts data is not datatables.
  expect_error(call_DTU(annot= sim$annot, count_data_A= c("not", "a", "list"), count_data_B= sim$boots_B[[1]], verbose = FALSE, dbg= "prep"), "counts data are not data.tables")
  expect_error(call_DTU(annot= sim$annot, count_data_A= sim$boots_A[[1]], count_data_B= c("not", "a", "list"), verbose = FALSE, dbg= "prep"), "counts data are not data.tables")
  expect_error(call_DTU(annot= sim$annot, count_data_A= data.frame(a="not", b="a", c="table"), count_data_B= sim$boots_B[[1]], verbose = FALSE, dbg= "prep"), "counts data are not data.tables")
  expect_error(call_DTU(annot= sim$annot, count_data_A= sim$boots_A[[1]], count_data_B= data.frame(a="not", b="a", c="table"), verbose = FALSE, dbg= "prep"), "counts data are not data.tables")
})

#==============================================================================
test_that("Bootstrap parameters are checked", {
  sim <- sim_boot_data(clean=TRUE)
  
  # Number of bootstraps.
  expect_warning(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, verbose = TRUE, dbg= "prep"),
                 "few bootstrap iterations", fixed= TRUE)
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, qbootnum = -5, verbose = FALSE, dbg= "prep"),
               "Invalid number of bootstraps", fixed= TRUE)
  expect_warning(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, qbootnum = 5, verbose = TRUE, dbg= "prep"),
                 "qbootnum is low", fixed= TRUE)
  expect_warning(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, qbootnum = 5000, verbose = TRUE, dbg= "prep"),
                 "number of quantification bootstraps is very high", fixed= TRUE)
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, qbootnum = 5000000000000, verbose = FALSE, dbg= "prep"),
               "number of quantification bootstraps would exceed the maximum capacity", fixed= TRUE)

  # Bootstraps without boot data.
  expect_warning(call_DTU(annot= sim$annot, count_data_A= sim$boots_A[[1]], count_data_B= sim$boots_B[[1]], qboot=TRUE, qbootnum=2, verbose= TRUE, dbg= "prep"), 
               "qboot is TRUE but no bootstrapped estimates", fixed= TRUE)
  
  # Confidence threshold.
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, rrep_thresh = -2, verbose = FALSE, dbg= "prep"),
               "Invalid reproducibility threshold", fixed= TRUE)
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, rrep_thresh = 2, verbose = FALSE, dbg= "prep"),
               "Invalid reproducibility threshold", fixed= TRUE)
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, qrep_thresh = -2, verbose = FALSE, dbg= "prep"),
               "Invalid reproducibility threshold", fixed= TRUE)
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, qrep_thresh = 2, verbose = FALSE, dbg= "prep"),
               "Invalid reproducibility threshold", fixed= TRUE)
})

#==============================================================================
test_that("Testing parameters are checked", {
  sim <- sim_boot_data(clean=TRUE)
  
  # Correction method.
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, correction= "wrong_name", verbose= FALSE, dbg= "prep"),
               "Invalid p-value correction method name", fixed= TRUE)
  
  # Tests.
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, testmode="GCSE", verbose = FALSE, dbg= "prep"),
               "Unrecognized value for testmode", fixed= TRUE)
  expect_silent(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, testmode="genes", verbose = FALSE, rboot=FALSE, qboot = FALSE, dbg= "prep"))
  expect_silent(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, testmode="transc", verbose = FALSE, rboot=FALSE, qboot = FALSE, dbg= "prep"))
  
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, qboot="GCSE", verbose = FALSE, dbg= "prep"),
               "Unrecognized value for qboot", fixed= TRUE)
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, rboot="GCSE", verbose = FALSE, dbg= "prep"),
               "Unrecognized value for rboot", fixed= TRUE)
  
  # Sums or means
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, verbose = FALSE, use_sums="sums", dbg= "prep"),
               "Invalid value for use_sums", fixed= TRUE)
})
  
#==============================================================================
test_that("Thresholds are checked", {
  sim <- sim_boot_data(clean=TRUE)

    # Probability threshold.
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, p_thresh = 666, verbose = FALSE, dbg= "prep"),
               "Invalid p-value threshold", fixed= TRUE)
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, p_thresh = -0.05, verbose = FALSE, dbg= "prep"),
               "Invalid p-value threshold", fixed= TRUE)
  
  # Read counts threshold.
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, abund_thresh = -5, verbose = FALSE, dbg= "prep"),
               "Invalid abundance threshold", fixed= TRUE)
  
  # Proportion change threshold.
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, dprop_thresh = -2, verbose = FALSE, dbg= "prep"),
               "Invalid proportion difference threshold", fixed= TRUE)
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, dprop_thresh = 2, verbose = FALSE, dbg= "prep"),
               "Invalid proportion difference threshold", fixed= TRUE)
})  
  
#==============================================================================
test_that("Scaling parameters are checked", {
  sim <- sim_boot_data(clean=TRUE)
    
  # Scaling of abundances.
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, verbose = FALSE, dbg= "prep", scaling=0),
               "Invalid scaling factor", fixed= TRUE)
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, verbose = FALSE, dbg= "prep", scaling=-33),
               "Invalid scaling factor", fixed= TRUE)
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, verbose = FALSE, dbg= "prep", scaling=c(2,3,4)),
               "Invalid scaling factor", fixed= TRUE)
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, verbose = FALSE, dbg= "prep", scaling=c(1,2,3,4,5)),
               "Invalid scaling factor", fixed= TRUE)
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, verbose = FALSE, dbg= "prep", scaling=c(1,0,2,3,4)),
               "Invalid scaling factor", fixed= TRUE)
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, verbose = FALSE, dbg= "prep", scaling=c(1,2,-3,4,5)),
               "Invalid scaling factor", fixed= TRUE)
})

#==============================================================================
test_that("General behaviour parameters are checked", {
  sim <- sim_boot_data(clean=TRUE)
  
  # Number of threads.
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, threads = 0, verbose= FALSE, dbg= "prep"),
               "Invalid number of threads", fixed= TRUE)
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, threads = -4, verbose= FALSE, dbg= "prep"),
               "Invalid number of threads", fixed= TRUE)
  expect_silent(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, threads = parallel::detectCores(logical=FALSE), verbose= FALSE, dbg= "prep"))
  expect_silent(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, threads = parallel::detectCores(logical=TRUE), verbose= FALSE, dbg= "prep"))
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, threads = parallel::detectCores(logical=TRUE) + 1, verbose= FALSE, dbg= "prep"),
               "threads exceed", fixed= TRUE)
  
  # Verbose is bool.
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, verbose="yes", dbg= "prep"),
               "Must be TRUE/FALSE", fixed= TRUE)
  
  # Lean is bool.
  expect_error(call_DTU(annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, name_A= name_A, name_B= name_B, lean="yes", dbg= "prep"),
               "Must be TRUE/FALSE", fixed= TRUE)
})

#EOF
