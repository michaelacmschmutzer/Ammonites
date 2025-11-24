# Compare the explanatory power of different extinction hypotheses. 
setwd("~/Documents/Projects/Ammonites/code/")

library(dplyr)

################### TODO Compile hypothesis testing data frame 

# Read in all required data
cepha <- read.csv('../data/cephalopods.csv')
body_sizes <- read.csv(
  '../data/doi_10_5061_dryad_zpc866t5b__v20210211/genus.sizes.ranges.rev.csv')
geo_nauti <- read.csv(
  '../results/geographic_distributions/nautilids_distributions_genus.csv')
geo_ammon <- read.csv(
  '../results/geographic_distributions/ammonoids_distributions_genus.csv')
abun_nauti <- read.csv(
  '../results/anundance_survivalship/nautilids_abundance_raw.csv')
abun_ammon <- read.csv(
  '../results/anundance_survivalship/ammonoids_abundance_raw.csv')
hatch_nauti <- read.csv('../data/nautilids_embryonic_shell_size.csv')
hatch_ammon <- read.csv('../data/ammonoids_embryonic_shell_size.csv')
surv_nauti <- read.csv('../data/nautilids_extinction_genus.csv')
surv_ammon <- read.csv('../data/ammonoids_extinction_genus.csv')

# Assemble the following data frame:
# Order Subord Supfam Fam Genus Geo.range Body.vol Hatch.size Abundance Survival

genus_data <- cepha %>%
# TODO Resolve taxonomy
#  select(c('order', 'family', 'genus')) %>%
  select(genus) %>%
  distinct()

################### TODO Check multivariate normality

################### TODO Check equality covariance matrix

################### TODO Hotellinger's t-test

################### TODO Estimating effect size