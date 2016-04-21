library(rats)
context("Transcript proportions")

test_that("Stats are correct", {

  # call calculate_tx_proportions and check result
  stats <- calculate_tx_proportions(psuedo_sleuth, mini_anno)

  m <- c(458.998606,41.249013,55.082330,12.535210,192.501678,411.072110,34.354291,46.807238,6.068015,221.341819)
  expect_equal(stats$means, m, tolerance = 1e-6)

  v <- c(476.89568,144.82521,43.95496,40.20190,772.10266,710.51429,62.14951,56.38875,11.68474,752.44403)
  expect_equal(stats$variance, v, tolerance = 1e-6)

})
