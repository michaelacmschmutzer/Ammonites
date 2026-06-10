# Compare strength of evidence for different hypotheses for all datasets
setwd("~/Documents/Projects/Ammonites/code/")

library(dplyr)
library(auRoc)

# Specify plotting aesthetics
darkblue  <- '#12255A'
orange    <- '#E38C59'
errbarcol <- '#CA282C'
font <- 'Palatino'
fsize <- 9

################### Compile hypothesis testing data set ################### 

hatch_nauti <- read.csv('../data/nautilids_embryonic_shell_size.csv')
hatch_ammon <- read.csv('../data/ammonoids_embryonic_shell_size.csv')
surv_nauti_genus <- read.csv('../data/nautilids_extinction_genus.csv')
surv_ammon_genus <- read.csv('../data/ammonoids_extinction_genus.csv')

# Simplify data, take median of species to get at genus level
hatch_nauti_genus <- hatch_nauti %>%
  group_by(genus) %>%
  summarise(med.hatching.size = median(hatching.size..mm.))
hatch_ammon_genus <- hatch_ammon %>%
  group_by(genus) %>%
  summarise(med.hatching.size = median(hatching.size..mm.))

surv_nauti_genus <- dplyr::select(surv_nauti_genus, c('genus', 'survival'))
surv_ammon_genus <- dplyr::select(surv_ammon_genus, c('genus', 'survival'))

# Merge with survival data
hatch_nauti_surv <- na.omit(merge(hatch_nauti_genus, surv_nauti_genus))
hatch_ammon_surv <- na.omit(merge(hatch_ammon_genus, surv_ammon_genus))

################### Permutation testing ################### 

# Use coin 'approximate' distribution This either uses exact tests or Monte Carlo
# subsampling (default 10000 random samples from all possible permutations, or
# all permutations if that is less)

# Get p-value
wt_hatch <- wilcox.test(med.hatching.size ~ factor(survival),
  data = hatch_ammon_surv, exact = TRUE)

surv <- hatch_ammon_surv[hatch_ammon_surv$survival == TRUE, 'med.hatching.size']
exti <- hatch_ammon_surv[hatch_ammon_surv$survival == FALSE, 'med.hatching.size']

surv <- hatch_nauti_surv[hatch_nauti_surv$survival == TRUE, 'med.hatching.size']
exti <- hatch_nauti_surv[hatch_nauti_surv$survival == FALSE, 'med.hatching.size']

# Get common language effect size (U / (n1 * n2)) with confidence interval
auc.nonpara.mw(exti, surv, method = 'DL.corr')

wilmanwhit.effect.size <- function(data, x_col, y_col = 'survival',
  ci.method = 'DL.corr', alpha = 0.05, nboot = 1000) {
  # Function to calculate the common language effect size with its associated
  # confidence interval
  groups <- split(data, data[[y_col]])
  e.size <- auc.nonpara.mw(groups[[1]][[x_col]], groups[[2]][[x_col]], method = ci.method, 
    conf.level = 1 - alpha, nboot = nboot)
  return(e.size)
}
