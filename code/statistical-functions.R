# Helper functions for statistical analyses and figure plotting

library(scales)

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

##################### PLOTTING FUNCTIONS ##################### 

# For notation in plotting.
# source: https://www.aj2duncan.com/blog/latex-notation-for-scientific-notation-in-ggplot2/
scientific_10 <- function(x) {
  scales::scientific_format()(x) %>% 
    str_remove_all("e\\+00") %>% # strip e+00 as we don't need it
    str_replace_all("e\\+", " %*% 10^") %>%
    parse(text = .)
}