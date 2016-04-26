library(rats)
context("Transcript proportions")

test_that("Stats are correct", {

  # call calculate_tx_proportions and check result
  # in theory uses lazy loading to get the data, but in practice not sure if it's working
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

#================================================================================
context("Annotation Filtering")

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

#================================================================================
context("Covariates to samples")

#================================================================================
test_that("Look-up tables creation", {
  result <- group_samples(pseudo_sleuth$sample_to_covariates)
  expect_length(result, 3)
  expect_equal(result$sample, list("Col-1"=1L, "Col-2"=2L, "Vir-1"=3L))
  expect_equal(result$condition, list("Col"=c(1L, 2L), "Vir"=c(3L)))
  expect_equal(result$batch, list("batch1"=c(1L, 2L, 3L)))
})


