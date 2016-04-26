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
