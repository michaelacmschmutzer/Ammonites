# Perform unit tests on statistical helper functions

library(testthat)

source('../../code/statistical-functions.R')

test_x <- c(1, 1, 5, 4, 6, 8, 6)
test_xna <- c(1, 1, 5, 4, 6, 8, 6, NA)
test_cor1 <- c(1, 2, 4, 3, 5, 7, 6)
test_cor2 <- c(2, 3, 5, 4, 1, 6, 7)

test_that(
  'Scaled median absolute deviations are calculated correctly', {
    expect_equal(scaled.median.deviation(test_x), 0.2)
    expect_true(is.na(scaled.median.deviation(test_xna, na.rm = FALSE)))
    expect_equal(scaled.median.deviation(test_xna, na.rm = TRUE), 0.2)
  }
)

test_that(
  'Spearman correlation returns correct rho', {
    expect_equal(round(get.spearmans.rho(test_cor1, test_cor2), digits = 7),
      0.6071429)
  }
)

