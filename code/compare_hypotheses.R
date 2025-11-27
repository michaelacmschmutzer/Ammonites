# Compare the explanatory power of different extinction hypotheses. 
setwd("~/Documents/Projects/Ammonites/code/")

library(dplyr)
library(DataExplorer)
library(corrplot)
library(mnt)

################### Compile hypothesis testing data set ################### 

# Read in all required data
cepha <- read.csv('../data/cephalopods.csv')
body_sizes <- read.csv(
  '../data/doi_10_5061_dryad_zpc866t5b__v20210211/genus.sizes.ranges.rev.csv')
geo_nauti <- read.csv(
  '../results/geographic_distributions/nautilids_distributions_genus.csv')
geo_ammon <- read.csv(
  '../results/geographic_distributions/ammonoids_distributions_genus.csv')
abun_nauti <- read.csv(
  '../results/abundance_survivalship/nautilids_abundance_raw.csv')
abun_ammon <- read.csv(
  '../results/abundance_survivalship/ammonoids_abundance_raw.csv')
hatch_nauti <- read.csv('../data/nautilids_embryonic_shell_size.csv')
hatch_ammon <- read.csv('../data/ammonoids_embryonic_shell_size.csv')
surv_nauti <- read.csv('../data/nautilids_extinction_genus.csv')
surv_ammon <- read.csv('../data/ammonoids_extinction_genus.csv')

# Simplify data, take median of species to get at genus level
hatch_nauti_genus <- hatch_nauti %>%
  select(genus, species, hatching.size..mm.) %>%
  group_by(genus) %>%
  summarise(med.hatching.size = median(hatching.size..mm.))
hatch_ammon_genus <- hatch_ammon %>%
  select(genus, species, hatching.size..mm.) %>%
  group_by(genus) %>%
  summarise(med.hatching.size = median(hatching.size..mm.))
body_sizes <- select(body_sizes, c('genus', 'logvol'))
surv_nauti <- select(surv_nauti, c('genus', 'survival'))
surv_ammon <- select(surv_ammon, c('genus', 'survival'))

# Combine nautilid and ammonoid data
geo <- rbind(geo_ammon, geo_nauti)
hatch <- rbind(hatch_ammon_genus, hatch_nauti_genus)
abun <- rbind(abun_ammon, abun_nauti)
surv <- rbind(surv_ammon, surv_nauti)

# Assemble the following data frame:
# Order Subord Supfam Fam Genus Geo.range Body.vol Hatch.size Abundance Survival

genus_data <- cepha %>%
# TODO Resolve taxonomy
#  select(c('order', 'family', 'genus')) %>%
  select(genus, is.nautilid) %>%
  distinct()

genus_data <- left_join(genus_data, geo, by = 'genus')
genus_data <- left_join(genus_data, body_sizes, by = 'genus')
genus_data <- left_join(genus_data, hatch, by = 'genus')
genus_data <- left_join(genus_data, abun, by = 'genus')
genus_data <- left_join(genus_data, surv, by = 'genus')

# Try logging the median hatching sizes and abundances
genus_data$log.hatching.size <- log10(genus_data$med.hatching.size)
genus_data$log.abun <- log10(genus_data$n)

################### Get to know dataset

# DataExplorer
introduce(genus_data)
plot_intro(genus_data)
plot_missing(genus_data)
plot_bar(genus_data)
plot_boxplot(genus_data, by = 'is.nautilid')
plot_boxplot(genus_data, by = 'survival')
plot_correlation(na.omit(genus_data), type = 'all')

pca_df <- na.omit(genus_data[ , c('PALEOMAP.area.km2', 'logvol', 
  'med.hatching.size', 'n', 'survival', 'is.nautilid')])

plot_prcomp(pca_df, variance_cap = 0.9, nrow = 2L, ncol = 2L)

# QQ plots
qq_data <- na.omit(genus_data[ ,  c('PALEOMAP.area.km2', 'logvol', 
  'med.hatching.size', 'n', 'log.hatching.size', 'log.abun')])
plot_qq(qq_data)

# Correlation plot
corrM <- cor(na.omit(genus_data[ , c('PALEOMAP.area.km2', 'logvol',
  'med.hatching.size', 'n', 'survival')]))
corrplot(corrM, order = 'AOE', type = 'lower')


################### Check multivariate normality

genus_data_cont <- as.matrix(
  na.omit(genus_data[ , c('PALEOMAP.area.km2', 'logvol', 'med.hatching.size',
  'n')]))

# Cox and Small (1978) multivariate normality test
test.CS(genus_data_cont)

################### TODO Check equality covariance matrix

################### TODO Hotellinger's t-test

################### TODO Estimating effect size