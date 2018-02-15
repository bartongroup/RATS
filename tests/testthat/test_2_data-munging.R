#==============================================================================
#==============================================================================
context("DTU Internal data munging")

#==============================================================================
test_that("Samples are grouped correctly", {
  sim <- sim_sleuth_data(cv_dt=FALSE)
  expect_silent(r <- group_samples(sim$slo$sample_to_covariates))
  # number of covariates
  expect_equal(length(r), length(sim$slo$sample_to_covariates))
  # names of covariates
  expect_named(r, names(sim$slo$sample_to_covariates))
  # number of values of each covariate
  expect_equal(length(r[[1]]), length(levels(as.factor(sim$slo$sample_to_covariates[[1]]))))
  expect_equal(length(r[[2]]), length(levels(as.factor(sim$slo$sample_to_covariates[[2]]))))
  # total number of samples
  expect_equal(sum(sapply(r[[1]],length)), length(sim$slo$sample_to_covariates[[1]]))
  expect_equal(sum(sapply(r[[2]],length)), length(sim$slo$sample_to_covariates[[2]]))

  # Repeat tests with data.table covariate
  sim <- sim_sleuth_data(cv_dt=TRUE)
  expect_silent(r <- group_samples(sim$slo$sample_to_covariates))
  # number of covariates
  expect_equal(length(r), length(sim$slo$sample_to_covariates))
  # names of covariates
  expect_named(r, names(sim$slo$sample_to_covariates))
  # number of values of each covariate
  expect_equal(length(r[[1]]), length(levels(as.factor(sim$slo$sample_to_covariates[[1]]))))
  expect_equal(length(r[[2]]), length(levels(as.factor(sim$slo$sample_to_covariates[[2]]))))
  # total number of samples
  expect_equal(sum(sapply(r[[1]],length)), length(sim$slo$sample_to_covariates[[1]]))
  expect_equal(sum(sapply(r[[2]],length)), length(sim$slo$sample_to_covariates[[2]]))

})


#==============================================================================
test_that("Bootstrapped counts are extracted correctly", {
  samples <- c(1,3)
  bst <- "id"
  cnt <- "counts"
  sim <- sim_sleuth_data(COUNTS_COL = cnt, BS_TARGET_COL = bst)
  lr <- denest_sleuth_boots(sim$slo, sim$annot, samples, cnt, bst)

  for (i in 1:length(lr)) {
    # The transcripts supposed to be there are there.
    expect_true(all(sim$isx %in% as.character(lr[[i]]$target_id)))
    # No NA.
    expect_false(any(is.na(lr[[i]])))

    # Number of bootstraps per sample.
    expect_equal(length(lr[[i]]) - 1, length(sim$slo$kal[[samples[i]]]$bootstrap))  #  one column in lr[[i]] is IDs.

    # All target counts pulled from the correct bootstraps and the correct transcripts.
    fltr2 <- match(sim$isx, lr[[i]][["target_id"]])                           # Where in the extracted counts are the expected IDs.
    counts_ok <- sapply(2:(length(lr[[i]])), function(j) {
      fltr1 <- match(sim$isx, sim$slo$kal[[samples[i]]]$bootstrap[[j-1]][[bst]])  # Where in the boot are the expected IDs.
      all(sim$slo$kal[[samples[i]]]$bootstrap[[j-1]][[cnt]][fltr1] == lr[[i]][[j]][fltr2])  # Both vectors' elements are ordered by the same IDs.
    })
    expect_true(all(counts_ok))

    # IDs in annotation, but not in bootstraps, should be 0.
    missing_from_boots_ok <- sapply(2:(length(lr[[i]])), function(j) {
      nib <- setdiff(sim$annot$target_id, sim$isx)
      nib <- match(nib, lr[[i]][["target_id"]])
      all(lr[[i]][[j]][nib] == 0)
    })
    expect_true(all(missing_from_boots_ok))

    # IDs in bootstrap, but not in annotation, should be absent.
    missing_from_annot_ok <- sapply(1:(length(lr[[i]]) - 1), function(j) {
      nia <- setdiff(sim$slo$kal[[i]]$bootstrap[[j]], sim$isx)
      any(nia %in% lr[[i]][["target_id"]])
    })
    expect_false(any(missing_from_annot_ok))
  }
})


#==============================================================================
test_that("The number of iterations is detected correctly", {
  sim <- sim_boot_data()
  expect_equal(infer_bootnum(sim$boots_A, sim$boots_B), 2)
})

