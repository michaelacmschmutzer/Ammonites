# Compare the explanatory power of different extinction hypotheses. 
setwd("~/Documents/Projects/Ammonites/code/")

library(dplyr)
library(DataExplorer)
library(corrplot)
library(QuantPsyc)
library(MVTests)
library(ggplot2)
library(gridExtra)
library(latex2exp)

# Specify plotting aesthetics
darkblue  <- '#12255A'
orange    <- '#E38C59'
errbarcol <- '#CA282C'
font <- 'Palatino'
fsize <- 9


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
  group_by(genus) %>%
  summarise(med.hatching.size = median(hatching.size..mm.))
hatch_ammon_genus <- hatch_ammon %>%
  group_by(genus) %>%
  summarise(med.hatching.size = median(hatching.size..mm.))
body_sizes <- dplyr::select(body_sizes, c(genus, logvol))
surv_nauti <- dplyr::select(surv_nauti, c('genus', 'survival'))
surv_ammon <- dplyr::select(surv_ammon, c('genus', 'survival'))

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
  dplyr::select(genus, is.nautilid) %>%
  distinct()

genus_data <- left_join(genus_data, geo, by = 'genus')
genus_data <- left_join(genus_data, body_sizes, by = 'genus')
genus_data <- left_join(genus_data, hatch, by = 'genus')
genus_data <- left_join(genus_data, abun, by = 'genus')
genus_data <- left_join(genus_data, surv, by = 'genus')

# Try logging the median hatching sizes and abundances
genus_data$log.hatching.size <- log10(genus_data$med.hatching.size)
genus_data$log.abun <- log10(genus_data$n)
genus_data$log.area <- log10(genus_data$PALEOMAP.area.km2)

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
qq_data <- na.omit(genus_data[genus_data$is.nautilid == TRUE,
  c('PALEOMAP.area.km2', 'logvol', 'med.hatching.size', 'n')])
qq_data_log <- na.omit(genus_data[genus_data$is.nautilid == TRUE,
  c('log.area', 'logvol', 'log.hatching.size', 'log.hatching.size',
  'log.abun')])
plot_qq(qq_data)
plot_qq(qq_data_log)

qq_data <- na.omit(genus_data[genus_data$is.nautilid == FALSE,
  c('PALEOMAP.area.km2', 'logvol', 'med.hatching.size', 'n')])
qq_data_log <- na.omit(genus_data[genus_data$is.nautilid == FALSE,
  c('log.area', 'logvol', 'log.hatching.size','log.abun')])
plot_qq(qq_data)
plot_qq(qq_data_log)

# Correlation plot
corrM <- cor(na.omit(genus_data[genus_data$is.nautilid == FALSE,
  c('log.area', 'logvol', 'log.hatching.size', 'log.abun', 
  'survival')]))
amcorr <- corrplot(corrM, order = 'original', type = 'lower', diag = FALSE)

corrM <- cor(na.omit(genus_data[genus_data$is.nautilid == TRUE,
  c('log.area', 'logvol', 'log.hatching.size', 'log.abun',
  'survival')]))
naucorr <- corrplot(corrM, order = 'original', type = 'lower', diag = FALSE)

nautilid_data <- genus_data[genus_data$is.nautilid == TRUE, ]
ammonoid_data <- genus_data[genus_data$is.nautilid == FALSE, ]

ammaa <- ggplot(data = ammonoid_data, aes(x = log.area, y = log.abun)) +
  geom_point(color = darkblue) + 
  labs(
    x = TeX('Geographic range ($log_{10}$ $km^2$)'),
    y = TeX('Abundance ($log_{10}$)')) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.y = element_text(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font),
  )
nauaa <- ggplot(data = nautilid_data, aes(x = log.area, y = log.abun)) +
  geom_point(color = orange) + 
  labs(
    x = TeX('Geographic range ($log_{10}$ $km^2$)'),
    y = TeX('Abundance ($log_{10}$)')) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.y = element_text(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font),
  )
p <- grid.arrange(ammaa, nauaa, ncol = 2)
ggsave('../results/comparing_hypotheses/Area_abundance_genus.png',
       width = 12, height = 6, units = 'cm', dpi = 600, plot = p)

cor(nautilid_data$log.area, nautilid_data$log.abun, method = 'spearman')
cor(ammonoid_data$log.area, ammonoid_data$log.abun, method = 'spearman')

plot(genus_data[genus_data$is.nautilid == FALSE, 'log.abun'],
     genus_data[genus_data$is.nautilid == FALSE, 'log.area'])

################### Check multivariate normality ###################

genus_data_cont <- na.omit(
  genus_data[genus_data$is.nautilid == FALSE, c('log.area', 'logvol',
  'log.hatching.size', 'log.abun')])

# Mardia's test. Both alternative hypotheses must be rejected
mult.norm(genus_data_cont)$mult.test

################### Check homogeneity covariance matrix ################### 

# Box M test. Test for homogeneity of the covariance matrices of extinct vs 
# surviving genera

survs <- na.omit(
  genus_data[genus_data$is.nautilid == FALSE, c('log.area', 'logvol',
  'log.hatching.size', 'log.abun', 'survival')])$survival

BoxM(genus_data_cont, survs)

################### Hotelling's t-test ASSUMPTIONS NOT MET

################### TODO Estimating effect size