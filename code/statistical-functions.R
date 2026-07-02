# Helper functions for statistical analyses and figure plotting

library(scales)
library(auRoc)

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

wilmanwhit.effect.size <- function(data, x_col, y_col = 'survival',
  ci.method = 'DL.corr', alpha = 0.05, nboot = 1000) {
  # Function to calculate the common language effect size with its associated
  # confidence interval.
  # Variables:
  #  data       = Data frame containing values in x_col and grouping factor in 
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
  
  groups <- split(data, data[[y_col]])
  e.size <- auc.nonpara.mw(groups[[1]][[x_col]], groups[[2]][[x_col]], 
    method = ci.method, conf.level = 1 - alpha, nboot = nboot)
  return(e.size)
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