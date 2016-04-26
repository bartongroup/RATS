library(rats)
context("Transcript proportions")

#==============================================================================
test_that("Stats are correct", {

  # call calculate_tx_proportions and check result
  # in theory uses lazy loading to get the data, but not sure if it's working
  stats <- calculate_DTU(pseudo_sleuth, mini_anno)

  ids <- c("AT1G03180.1", "AT1G03180.2", "AT1G01020.2", "AT1G01020.1")
  m <- c(55.08233, 12.53521, 41.24901, 458.99861)
  v <- c(43.95496, 40.20190, 144.82521, 476.89568)
  p <- c(0.81461, 0.18538, 0.08245, 0.91754)
  test <- data.frame(ids,m,v,p)

  mapping <- match(test$ids, stats$Col[[TARGET_ID]])

  expect_equal(stats$Col$mean[mapping], test$m, tolerance = 1e-6)
  expect_equal(stats$Col$variance[mapping], test$v, tolerance = 1e-6)
  expect_equal(stats$Col$proportion[mapping], test$p, tolerance = 1e-5)

  m <- c(46.807238, 6.068015, 34.354291, 411.072110)
  v <- c(56.38875, 11.68474, 62.14951, 710.51429)
  p <- c(0.88523, 0.11476, 0.07712, 0.92287)
  test <- data.frame(ids,m,v,p)

  mapping <- match(test$ids, stats$Vir[[TARGET_ID]])

  expect_equal(stats$Vir$mean[mapping], test$m, tolerance = 1e-6)
  expect_equal(stats$Vir$variance[mapping], test$v, tolerance = 1e-6)
  expect_equal(stats$Vir$proportion[mapping], test$p, tolerance = 1e-5)
})

#==============================================================================
test_that("Mixed order bootstraps give same results as unmixed", {

  # make a pseudo sleuth object with mixed up bootstrap rows
  mixed_psuedo_sleuth <- pseudo_sleuth
  mixed_pseudo_sleuth$kal[[1]]$bootstrap[[3]] <- mixed_pseudo_sleuth$kal[[1]]$bootstrap[[3]][c(1,3,2,5,4),]
  mixed_pseudo_sleuth$kal[[1]]$bootstrap[[4]] <- mixed_pseudo_sleuth$kal[[1]]$bootstrap[[4]][c(5,4,3,2,1),]
  mixed_pseudo_sleuth$kal[[3]]$bootstrap[[1]] <- mixed_pseudo_sleuth$kal[[3]]$bootstrap[[1]][c(3,5,1,4,2),]

  mixed_stats <- calculate_DTU(mixed_pseudo_sleuth, mini_anno)
  unmixed_stats <- calculate_DTU(pseudo_sleuth, mini_anno)

  expect_equal(sort(mixed_stats$Col$mean), sort(unmixed_stats$Col$mean), tolerance = 1e-6)
  expect_equal(sort(mixed_stats$Col$variance), sort(unmixed_stats$Col$variance), tolerance = 1e-6)
  expect_equal(sort(mixed_stats$Col$proportion), sort(unmixed_stats$Col$proportion), tolerance = 1e-6)
  expect_equal(sort(mixed_stats$Vir$mean), sort(unmixed_stats$Vir$mean), tolerance = 1e-6)
  expect_equal(sort(mixed_stats$Vir$variance), sort(unmixed_stats$Vir$variance), tolerance = 1e-6)
  expect_equal(sort(mixed_stats$Vir$proportion), sort(unmixed_stats$Vir$proportion), tolerance = 1e-6)
})

#================================================================================
context("Annotation Preparation")

#================================================================================
test_that("Sibling IDs identification (Kimon's version)", {
  result <- mark_sibling_targets2(mini_anno)
  reference <- mini_anno
  reference["has_siblings"] <- c(FALSE,TRUE,TRUE,TRUE,TRUE)
  expect_identical(result, reference)
})

test_that("Sibling IDs identification (Nick's version)", {
  generalized_anno <- mini_anno
  generalized_anno$target_id <- c("t1","t2","t3","t4","t5")
  result <- mark_sibling_targets(generalized_anno)
  reference <- mini_anno
  reference["has_siblings"] <- c(FALSE,TRUE,TRUE,TRUE,TRUE)
  reference$target_id <- c("t1","t2","t3","t4","t5")
  expect_identical(result, reference)
})

test_that("Reverse covariates look-up tables creation", {
  result <- group_samples(pseudo_sleuth$sample_to_covariates)
  expect_length(result, 3)
  expect_equal(result$sample, list("Col-1"=1L, "Col-2"=2L, "Vir-1"=3L))
  expect_equal(result$condition, list("Col"=c(1L, 2L), "Vir"=c(3L)))
  expect_equal(result$batch, list("batch1"=c(1L, 2L, 3L)))
})

#================================================================================
context("Counts-based calculations")

#================================================================================
test_that("Row-wise means and variances (Kimon's version)", {
  df <- data.frame("a"=c(1,20,300), "b"=c(2,21,301), "c"=c(2,40,600), "d"=c(3, 30, 450))
  rownames(df) <- c("A", "B", "C")
  result <- calculate_stats2(df)
  reference <- data.frame("mean"=c(mean(c(1,2,2,3)), mean(c(20,21,40,30)), mean(c(300,301,600,450))),
                          "variance"=c(var(c(1,2,2,3)), var(c(20,21,40,30)), var(c(300,301,600,450))))
  rownames(reference) <- c("A", "B", "C")
  expect_identical(result, reference)
})

test_that("Row-wise means and variances (Kira's version)", {
  df <- data.frame("a"=c(1,20,300), "b"=c(2,21,301), "c"=c(2,40,600), "d"=c(3, 30, 450))
  rownames(df) <- c("A", "B", "C")
  result <- calculate_stats(df)
  reference <- data.frame("mean"=c(mean(c(1,2,2,3)), mean(c(20,21,40,30)), mean(c(300,301,600,450))),
                          "variance"=c(var(c(1,2,2,3)), var(c(20,21,40,30)), var(c(300,301,600,450))))
  rownames(reference) <- c("A", "B", "C")
  expect_identical(result, reference)
})

