#==============================================================================
#==============================================================================
context("DTU internal data munging")

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
  lr <- denest_sleuth_boots(sim$slo, sim$annot[[1]], samples, cnt, bst)
  
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
test_that("Filters work correctly", {
  
  # !!! This test is tightly dependent on the data used for the test and the default parameter values, in order to
  # !!! ensure correct response to specific scenarios.
  
  sim <- sim_sleuth_data(cnames=c("ONE","TWO"))
  mydtu <- call_DTU(annot= sim$annot, slo= sim$slo, name_A= "ONE", name_B= "TWO", abund_thresh=10, dprop_thresh=0.1, verbose = FALSE, rboot=FALSE, qboot = FALSE)
  
  expect_equivalent(as.list(mydtu$Genes["1A1N", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(1, 1, 0, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Genes["1B1C", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(2, 1, 0, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Genes["1D1C", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(2, 1, 0, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Genes["ALLA", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(1, 1, 0, FALSE, NA))
  expect_equivalent(as.list(mydtu$Genes["ALLB", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(2, 2, 0, FALSE, NA))
  expect_equivalent(as.list(mydtu$Genes["CC", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(2, 2, 2, TRUE, TRUE))
  expect_equivalent(as.list(mydtu$Genes["LC", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(2, 2, 1, FALSE, TRUE))
  expect_equivalent(as.list(mydtu$Genes["MIX6", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(6, 5, 5, TRUE, TRUE))
  expect_equivalent(as.list(mydtu$Genes["NIB", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(1, 0, 0, FALSE, NA))
  expect_equivalent(as.list(mydtu$Genes["NN", list(known_transc, detect_transc, elig_transc, elig, elig_fx)]), 
                    list(2, 2, 2, TRUE, FALSE))
  
  setkey(mydtu$Transcripts, target_id)
  expect_equivalent(as.list(mydtu$Transcripts["1A1N-2", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Transcripts["1B1C.1", list(elig_xp, elig, elig_fx)]),
                    list(FALSE, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Transcripts["1B1C.2", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Transcripts["1D1C:one", list(elig_xp, elig, elig_fx)]),
                    list(FALSE, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Transcripts["1D1C:two", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Transcripts["ALLA1", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, FALSE, NA))
  expect_equivalent(as.list(mydtu$Transcripts["ALLB1", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, FALSE, NA))
  expect_equivalent(as.list(mydtu$Transcripts["ALLB2", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, FALSE, NA))
  expect_equivalent(as.list(mydtu$Transcripts["CC_a", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, TRUE))
  expect_equivalent(as.list(mydtu$Transcripts["CC_b", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, TRUE))
  expect_equivalent(as.list(mydtu$Transcripts["LC1", list(elig_xp, elig, elig_fx)]),
                    list(FALSE, FALSE, TRUE))
  expect_equivalent(as.list(mydtu$Transcripts["LC2", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, TRUE))
  expect_equivalent(as.list(mydtu$Transcripts["MIX6.c1", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, TRUE))
  expect_equivalent(as.list(mydtu$Transcripts["MIX6.c2", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, TRUE))
  expect_equivalent(as.list(mydtu$Transcripts["MIX6.c3", list(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, FALSE))
  expect_equivalent(as.list(mydtu$Transcripts["MIX6.c4", .(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, TRUE))
  expect_equivalent(as.list(mydtu$Transcripts["MIX6.d", .(elig_xp, elig, elig_fx)]),
                    list(FALSE, FALSE, FALSE))
  expect_equivalent(as.list(mydtu$Transcripts["MIX6.nc", .(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, FALSE))
  expect_equivalent(as.list(mydtu$Transcripts["NIB.1", .(elig_xp, elig, elig_fx)]),
                    list(FALSE, FALSE, NA))
  expect_equivalent(as.list(mydtu$Transcripts["1NN", .(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, FALSE))
  expect_equivalent(as.list(mydtu$Transcripts["2NN", .(elig_xp, elig, elig_fx)]),
                    list(TRUE, TRUE, FALSE))
})

#==============================================================================
test_that("The number of iterations is detected correctly", {
  sim <- sim_sleuth_data()
  expect_equal(infer_bootnum(sim$slo, NULL, NULL), 2)
  
  sim <- sim_boot_data()
  expect_equal(infer_bootnum(NULL,sim$boots_A, sim$boots_B), 2)
})

