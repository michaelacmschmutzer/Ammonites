# Compare strength of evidence for different hypotheses for all datasets
setwd("~/Documents/Projects/Ammonites/code/")

library(dplyr)
library(auRoc)

source('statistical-functions.R')

# Specify plotting aesthetics
darkblue  <- '#12255A'
orange    <- '#E38C59'
errbarcol <- '#CA282C'
font <- 'Palatino'
fsize <- 9

################### Compile hypothesis testing data set ################### 
body_sizes <- read.csv(
  '../data/doi_10_5061_dryad_zpc866t5b__v20210211/genus.sizes.ranges.rev.csv')
geo_nauti <- read.csv(
  '../results/genus/geographic_distributions/nautilids_distributions_genus.csv')
geo_ammon <- read.csv(
  '../results/genus/geographic_distributions/ammonoids_distributions_genus.csv')
geo_nauti_boot <- read.csv(
  '../results/genus/subsampling_distributions/nautilids/bootstrap.csv')
geo_ammon_boot <- read.csv(
  '../results/genus/subsampling_distributions/ammonoids/bootstrap.csv')
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
body_sizes <- dplyr::select(body_sizes, c('genus', 'logvol'))
surv_nauti_genus <- dplyr::select(surv_nauti_genus, c('genus', 'survival'))
surv_ammon_genus <- dplyr::select(surv_ammon_genus, c('genus', 'survival'))

# Merge with survival data
hatch_nauti_surv <- na.omit(merge(hatch_nauti_genus, surv_nauti_genus))
hatch_ammon_surv <- na.omit(merge(hatch_ammon_genus, surv_ammon_genus))

################### Calculating effect sizes ################### 

# Get p-value
wt_hatch <- wilcox.test(med.hatching.size ~ factor(survival),
  data = hatch_ammon_surv, exact = TRUE)

surv <- hatch_ammon_surv[hatch_ammon_surv$survival == TRUE, 'med.hatching.size']
exti <- hatch_ammon_surv[hatch_ammon_surv$survival == FALSE, 'med.hatching.size']

surv <- hatch_nauti_surv[hatch_nauti_surv$survival == TRUE, 'med.hatching.size']
exti <- hatch_nauti_surv[hatch_nauti_surv$survival == FALSE, 'med.hatching.size']

# Get common language effect size (U / (n1 * n2)) with confidence interval
variables <- c(
  'Hatching size (genus)',
  'Body size (genus)',
  'Geographic range (genus)',
  'Geographic range (species)',
  'Geographic range (bootstrap, genus)',
  'Geographic range (boostrap, species)',
  'Geographic range (jackknife, genus)',
  'Geographic range (jackknife, species)'
)

effect.sizes.ammon <- data.frame(
  variable = variables, eff.size = 0, ci.lower = 0, ci.upper = 0)
effect.sizes.nauti <- data.frame(
  variable = variables, eff.size = 0, ci.lower = 0, ci.upper = 0)


effect.sizes.ammon[1, 2:4] <-
  wilmanwhit.effect.size(hatch_ammon_surv, 'med.hatching.size')
effect.sizes.ammon[2, 2:4] <- 
