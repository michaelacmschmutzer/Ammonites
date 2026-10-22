# Helper functions for statistical analyses and figure plotting

library(stats)
library(scales)
library(auRoc)
library(spatstat)
library(MASS)

##################### STATISTICS FUNCTIONS ##################### 

scaled.median.deviation <- function(x, na.rm = TRUE) {
  # Calculates a robust non-parametric variant of the coefficient of variation.
  # It is based on the median absolute deviation (MAD), which calculates the 
  # absolute deviation of a data point from the median. Dividing the MAD by
  # the median results in a scale-free, unit-less estimator of dispersion.
  
  # Calculate the median
  x_med <- median(x, na.rm = na.rm)
  
  # Calculate the median absolute deviance.
  mad <- median(abs(x - x_med), na.rm = na.rm)
  
  # Calculate and return the scaled MAD
  return(mad / x_med)
}

get.spearmans.rho <- function(x, y) {
  # Wrapper function in order to obtain only Spearman's rho and nothing else
  rho <- cor.test(x, y, method = 'spearman')$estimate
  return(unname(rho))
}

wilmanwhit.effect.size <- function(data.fr, x_col, y_col = 'survival',
  ci.method = 'DL.corr', alpha = 0.05, nboot = 1000) {
  # Function to calculate the common language effect size with its associated
  # confidence interval.
  # Variables:
  #  data.fr    = Data frame containing values in x_col and grouping factor in 
  #               y_col
  #  x_col      = String. Column name containing values
  #  y_col      = String. Column name containing grouping factor
  #  ci.method  = Confidence interval method from auRoc. Must be one of 
  #               "newcombe", "pepe", "delong", "jackknife", "bootstrapP",
  #               "bootstrapBCa"
  #  alpha      = Alpha level
  #  nboot      = Number of bootstrap repetitions
  # Value:
  #   e.size    = Point estimate and lower and upper bounds of the confidence
  #               interval
  
  groups <- split(data.fr, data.fr[[y_col]])
  e.size <- auc.nonpara.mw(groups[[1]][[x_col]], groups[[2]][[x_col]], 
    method = ci.method, conf.level = 1 - alpha, nboot = nboot)
  return(e.size)
}

# example.data <- na.omit(genus_data[, c('log.area.zscore', 'logvol.zscore', 'log.hatching.size.zscore')])
# Function to fit a kernel density approximation to data, and then sample from
# the approximation
make.three.random.samples <- function(data.fr, n) {
  
  # Function to extract N observations from kernel density approximations
  # This function assumes a data structure with three variables of interest
  
  # Empirical covariance matrix of predictor variables
  cov.emp <- cor(data.fr, method = 'pearson')
  
  # Sample from data with replacements, and smooth out with a random 
  # Gaussian error with default bandwidth 'bw'
  den1 <- density(data.fr[ , 1])
  den2 <- density(data.fr[ , 2])
  den3 <- density(data.fr[ , 3])
  
  # Generate multivariate Gaussian distribution with desired correlation structure
  mult.norm <- mvrnorm(n, c(0, 0, 0), cov.emp)
  # Transform to multivariate Uniform distribution with corrs
  mult.unif <- pnorm(mult.norm, mean = 0, sd = 1)
  
  # Transform from uniform to empirically observed distributions using the 
  # quantile function of the fitted kernel density estimates. These should retain
  # the covariance structure of the empirical data
  X1 <- quantile(den1, probs = mult.unif[ , 1])
  X2 <- quantile(den2, probs = mult.unif[ , 2])
  X3 <- quantile(den3, probs = mult.unif[ , 3])
  return(list(X1, X2, X3))
}

# Empirical coefficient estimates
# b0 <- regr_all$coefficients[[1]]
# b1 <- regr_all$coefficients[[2]]
# b2 <- regr_all$coefficients[[3]]
# b3 <- regr_all$coefficients[[4]]

fit.logistic.regression <- function(coeffs, samples, alpha = 0.05) {
  # For the logistic regression, I have three explanatory variables:
  # X1 == Geographic range sizes
  # X2 == Hatching sizes
  # X3 == Body sizes
  # These have a certain correlative structure and are related to the response
  # variable Y via a logit-link function. The correlation between the three 
  # variables and Y is the effect size of interest (B1, B2, B3)
  
  # Get out the random samples
  X1 <- samples[[1]]
  X2 <- samples[[2]]
  X3 <- samples[[3]]
  
  # Get out the coefficients
  b0 <- coeffs[1]
  b1 <- coeffs[2]
  b2 <- coeffs[3]
  b3 <- coeffs[4]
  
  # Structure of logistic regression model. p is the vector of probabilities of 
  # the response variable Y, which in this case is extinction/survival
  p <- plogis(b0 + 
              b1 * X1 +
              b2 * X2 +
              b3 * X3)
  Y <- rbinom(n = length(X1), size = 1, prob = p)
  # Function to fit logistic regression to simulated data X times. Extract power.
  sim_data <- data.frame(y = Y, x1 = X1, x2 = X2, x3 = X3)
  regr_sim <- glm(y ~ x1 + x2 + x3, data = sim_data, family = binomial)
  
  # Check if coefficients are significant or not
  outcome <- c(NA, NA, NA)
  outcome[1] <- summary(regr_sim)$coefficients[2,4] < alpha
  outcome[2] <- summary(regr_sim)$coefficients[3,4] < alpha
  outcome[3] <- summary(regr_sim)$coefficients[4,4] < alpha
  
  return(outcome)
}

power.logistic.regression <- function(coeff.range, nsteps, data.fr, n, 
  nrepl = 100, alpha = 0.05) {
  # Function to perform power analysis for logistic regression with three
  # variables of interest
  # Variables:
  #  coeff.range  = Range of hypothetical coefficients to explore
  #  nsteps       = Number of intervals between hypothetical cofficients
  #  data.fr      = Empirical data to parameterise kernel density probabilities
  #  n            = Number of data points to simulate
  #  nrepl        = Number of replicates per set of coefficient values
  #  alpha        = Alpha level to test for
  # Value
  #  power        = Power at each set of coefficient values
  
  # define coefficient space
  
  # iterate over coefficient space. Perform nrepl replicates per coefficient set
  samples <- make.three.random.samples(data.fr = data.fr, n = n)
  outcome <- fit.logistic.regression(coeffs, samples, alpha = alpha)
}



##################### PLOTTING FUNCTIONS ##################### 

# For notation in plotting.
# source: https://www.aj2duncan.com/blog/latex-notation-for-scientific-notation-in-ggplot2/
scientific_10 <- function(x) {
  scales::scientific_format()(x) %>% 
    str_remove_all("e\\+00") %>% # strip e+00 as we don't need it
    str_replace_all("e\\+", " %*% 10^") %>%
    parse(text = .)
}