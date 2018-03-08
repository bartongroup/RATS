#==============================================================================
#==============================================================================
context("RATs main")

#==============================================================================
test_that("The paramater interpretations are correct", {
  sim1 <- sim_boot_data()
  sim2 <- sim_count_data()
  
  # Bootstraps or counts?
  expect_equal(2, call_DTU(dbg='prep', annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, verbose = FALSE)[[1]])
  expect_equal(1, call_DTU(dbg='prep', annot= sim2$annot, count_data_A= sim2$counts_A, count_data_B= sim2$counts_B, verbose = FALSE)[[1]])
  
  # How many bootstraps? The guess function is already tested separately. This confirms that the value finds its way through flow control.
  expect_equal(infer_bootnum(sim1$boots_A, sim1$boots_B), 
               call_DTU(dbg='prep', qboot=TRUE, qbootnum=0, annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, verbose = FALSE)[[2]])
  expect_equal(100, 
               call_DTU(dbg='prep', qboot=TRUE, qbootnum=100, annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, verbose = FALSE)[[2]])
  
  # Which tests?
  expect_true(call_DTU(dbg='prep', testmode= "both", annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, verbose = FALSE)[[3]])
  expect_true(call_DTU(dbg='prep', testmode= "both", annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, verbose = FALSE)[[4]])
  
  expect_true(call_DTU(dbg='prep', testmode= "transc", annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, verbose = FALSE)[[3]])
  expect_false(call_DTU(dbg='prep', testmode= "transc", annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, verbose = FALSE)[[4]])
  
  expect_false(call_DTU(dbg='prep', testmode= "genes", annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, verbose = FALSE)[[3]])
  expect_true(call_DTU(dbg='prep', testmode= "genes", annot= sim1$annot, boot_data_A= sim1$boots_A, boot_data_B= sim1$boots_B, verbose = FALSE)[[4]])
})

#==============================================================================
test_that("The index is received", {
  sim <- sim_boot_data()
  # tidy_annot() should be tested separately.
  expect_equal(tidy_annot(sim$annot), 
               call_DTU(dbg='indx', annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, verbose = FALSE))
  
})

#==============================================================================
test_that("Scaling is applied", {
  sim <- sim_boot_data()
  
  # Single factor.
  S <- 10
  bdt <- call_DTU(dbg='bootin', annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, verbose = FALSE, scaling=S)
  cdt <- call_DTU(dbg='countin', annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, verbose = FALSE, scaling=S)
  sdt <- call_DTU(dbg='scale', annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, verbose = FALSE, scaling=S)
  expect_equal(S * cdt[[1]], sdt[[1]])
  expect_equal(S * cdt[[2]], sdt[[2]])
  for (b in 1:length(bdt$bootA)) {
    expect_equal(S * bdt$bootA[[b]][,-1], sdt$bootA[[b]])
  }
  for (b in 1:length(bdt$bootB)) {
    expect_equal(S * bdt$bootB[[b]][,-1], sdt$bootB[[b]])
  }
  
  # Vector.
  S <- c(10, 20, 30, 40)
  bdt <- call_DTU(dbg='bootin', annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, verbose = FALSE, scaling=S)
  cdt <- call_DTU(dbg='countin', annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, verbose = FALSE, scaling=S)
  sdt <- call_DTU(dbg='scale', annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, verbose = FALSE, scaling=S)
  
  expect_equal(S[1] * cdt[[1]]$V1, sdt[[1]]$V1)
  expect_equal(S[2] * cdt[[1]]$V2, sdt[[1]]$V2)
  expect_equal(S[3] * cdt[[2]]$V1, sdt[[2]]$V1)
  expect_equal(S[4] * cdt[[2]]$V2, sdt[[2]]$V2)
  
  expect_equal(S[1] * bdt$bootA[[1]][,-1], sdt[[3]][[1]])
  expect_equal(S[2] * bdt$bootA[[2]][,-1], sdt[[3]][[2]])
  expect_equal(S[3] * bdt$bootB[[1]][,-1], sdt[[4]][[1]])
  expect_equal(S[4] * bdt$bootB[[2]][,-1], sdt[[4]][[2]])
})

test_that("Parameters are recorded", {
  sim <- sim_boot_data()
  param <- call_DTU(dbg='info', annot= sim$annot, boot_data_A= sim$boots_A, boot_data_B= sim$boots_B, scaling=c(1,2,3,4),
                    description='Test.',name_A='cA', name_B='cB', varname='testing',
                    p_thresh=0.001, abund_thresh=1, dprop_thresh=0.15, correction='bonferroni',
                    qboot=TRUE, qbootnum=99, qrep_thresh=0.8, rboot=TRUE, rrep_thresh=0.6, 
                    testmode="genes", seed=666, verbose = FALSE)$Parameters
  
  expect_equal(param$description, 'Test.')
  expect_true(!is.na(param$time))
  expect_equal(param$rats_version, packageVersion('rats'))
  expect_equal(param$R_version, R.Version()[c("platform", "version.string")])
  expect_equal(c(param$var_name, param$cond_A, param$cond_B), c('testing', 'cA', 'cB'))
  expect_equal(c(param$data_type, param$tests), c("bootstrapped abundance estimates", 'genes'))
  expect_equal(c(param$num_replic_A, param$num_replic_B), c(2, 2))
  expect_equal(c(param$num_genes, param$num_transc), c(11, 23))
  expect_equal(list(param$p_thresh, param$abund_thresh, param$dprop_thresh, param$correction, param$abund_scaling), list(0.001, 1, 0.15, 'bonferroni', c(1,2,3,4)))
  expect_equal(list(param$quant_boot, param$quant_bootnum, param$quant_reprod_thresh), list(TRUE, 99, 0.8))
  expect_equal(list(param$rep_boot, param$rep_bootnum, param$rep_reprod_thresh), list(TRUE, NA_integer_, 0.6))
})
