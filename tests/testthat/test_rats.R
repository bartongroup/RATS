library(rats)
context("Transcript proportions")

test_that("Stats are correct", {

  # call calculate_tx_proportions and check result
  # in theory uses lazy loading to get the data, but in practice not sure if it's working
  stats <- calculate_tx_proportions(pseudo_sleuth, mini_anno)

  m <- c(458.99861,41.24901,55.08233,12.53521)
  expect_equal(stats$Col$mean, m, tolerance = 1e-6)

  v <- c(476.89568,144.82521,43.95496,40.20190)
  expect_equal(stats$Col$variance, v, tolerance = 1e-6)

  m <- c(411.072110,34.354291,46.807238,6.068015)
  expect_equal(stats$Vir$mean, m, tolerance = 1e-6)

  v <- c(710.51429,62.14951,56.38875,11.68474)
  expect_equal(stats$Vir$variance, v, tolerance = 1e-6)

})
