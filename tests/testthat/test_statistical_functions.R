# Perform unit tests on statistical helper functions

library(testthat)

source('../../code/statistical-functions.R')
source('../../code/power_analysis.R')

test_x <- c(1, 1, 5, 4, 6, 8, 6)
test_xna <- c(1, 1, 5, 4, 6, 8, 6, NA)
test_cor1 <- c(1, 2, 4, 3, 5, 7, 6)
test_cor2 <- c(2, 3, 5, 4, 1, 6, 7)
test_df <- data.frame(
  values = c(1, 1, 3, 7, 2, 6, 8, 9, 4),
  survival = c(FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
)

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

test_that(
  'Common language effect size and CI are calculated correctly', {
    x <- wilmanwhit.effect.size(test_df, 'values')
    expect_equal(round(x[1], digits = 7), 0.05)
    expect_equal(round(x[2], digits = 7), 0.0020122)
    expect_equal(round(x[3], digits = 7), 0.5787490)
  }
)

test_that(
  'Normal distributions are parameterised correctly', {
    y.mean <- find_mean_shift(1-0.3273604, 1, 2, 3)
    expect_equal(round(y.mean, digits = 6), 2)
  }
)

test_that(
  'Statistical power is calculated correctly', {
    set.seed(1)
    power <- determine_power(0.8, 100, 10, 2, 3, 100)
    expect_equal(round(power, digits = 7), 0.92)
  }
)
