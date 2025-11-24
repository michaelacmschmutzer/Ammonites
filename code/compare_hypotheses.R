# Compare the explanatory power of different extinction hypotheses. 
setwd("~/Documents/Projects/Ammonites/code/")

library(dplyr)

# Command line arguments if run from bash
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0){
  cephas <- 'nautilids'  # 'nautilids' 'ammonoids'
} else {
  cephas <- args[1]
}
