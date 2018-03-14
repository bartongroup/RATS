#==============================================================================
#==============================================================================
context("DTU Internal data munging")


#==============================================================================
test_that("The number of iterations is detected correctly", {
  sim <- sim_boot_data()
  expect_equal(infer_bootnum(sim$boots_A, sim$boots_B), 2)
})

