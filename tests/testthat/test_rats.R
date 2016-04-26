library(rats)
context("Transcript proportions")

#==============================================================================
test_that("Stats are correct", {

  # call calculate_tx_proportions and check result
  # in theory uses lazy loading to get the data, but not sure if it's working
  stats <- calculate_DTU(pseudo_sleuth, mini_anno)

  m <- c(458.99861,41.24901,55.08233,12.53521)
  expect_equal(sort(stats$Col$mean), sort(m), tolerance = 1e-6)

  v <- c(476.89568,144.82521,43.95496,40.20190)
  expect_equal(sort(stats$Col$variance), sort(v), tolerance = 1e-6)

  m <- c(411.072110,34.354291,46.807238,6.068015)
  expect_equal(sort(stats$Vir$mean), sort(m), tolerance = 1e-6)

  v <- c(710.51429,62.14951,56.38875,11.68474)
  expect_equal(sort(stats$Vir$variance), sort(v), tolerance = 1e-6)
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
  expect_equal(sort(mixed_stats$Vir$mean), sort(unmixed_stats$Vir$mean), tolerance = 1e-6)
  expect_equal(sort(mixed_stats$Vir$variance), sort(unmixed_stats$Vir$variance), tolerance = 1e-6)
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
test_that("Row-wise means and variances", {
  df <- data.frame("a"=c(1,20,300), "b"=c(2,21,301), "c"=c(2,40,600), "d"=c(3, 30, 450))
  rownames(df) <- c("A", "B", "C")
  result <- calculate_stats2(df)
  reference <- data.frame("mean"=c(mean(c(1,2,2,3)), mean(c(20,21,40,30)), mean(c(300,301,600,450))),
                          "variance"=c(var(c(1,2,2,3)), var(c(20,21,40,30)), var(c(300,301,600,450))))
  rownames(reference) <- c("A", "B", "C")
  expect_identical(result, reference)
})


