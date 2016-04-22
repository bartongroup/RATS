library(rats)
context("Transcript proportions")

test_that("Stats are correct", {

  # call calculate_tx_proportions and check result
  # in theory uses lazy loading to get the data, but in practice not sure if it's working
  stats <- calculate_tx_proportions(pseudo_sleuth, mini_anno)

  m <- c(41.249013,55.082330,12.535210,192.501678,34.354291,46.807238,6.068015,221.341819)
  expect_equal(stats$means, m, tolerance = 1e-6)

  v <- c(144.82521,43.95496,40.20190,772.10266,62.14951,56.38875,11.68474,752.44403)
  expect_equal(stats$variance, v, tolerance = 1e-6)

})
