# Compare the explanatory power of different extinction hypotheses. 
setwd("~/Documents/Projects/Ammonites/code/")

library(dplyr)
library(DataExplorer)
library(corrplot)
library(ggplot2)
library(gridExtra)
library(latex2exp)
library(lme4)


# Specify plotting aesthetics
darkblue  <- '#12255A'
orange    <- '#E38C59'
errbarcol <- '#CA282C'
font <- 'Palatino'
fsize <- 9

################### Convenience functions ################### 

zscore <- function(x, na.rm = TRUE) {
  # Calculate the Z-score normalized values for x
  z <- (x - mean(x, na.rm = na.rm) ) / sd(x, na.rm = na.rm)
  return(z)
}

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
genus_data$sqrt.abun <- sqrt(genus_data$n)
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
  'sqrt.abun')])
plot_qq(qq_data)
plot_qq(qq_data_log)

qq_data <- na.omit(genus_data[genus_data$is.nautilid == FALSE,
  c('PALEOMAP.area.km2', 'logvol', 'med.hatching.size', 'n')])
qq_data_log <- na.omit(genus_data[genus_data$is.nautilid == FALSE,
  c('log.area', 'logvol', 'log.hatching.size','sqrt.abun')])
plot_qq(qq_data)
plot_qq(qq_data_log)

# Correlation plot
corrM <- cor(na.omit(genus_data[genus_data$is.nautilid == FALSE,
  c('log.area', 'logvol', 'log.hatching.size', 'sqrt.abun', 
  'survival')]), method = 'spearman')
amcorr <- corrplot(corrM, order = 'original', type = 'lower', diag = FALSE)

corrM <- cor(na.omit(genus_data[genus_data$is.nautilid == TRUE,
  c('log.area', 'logvol', 'log.hatching.size', 'sqrt.abun',
  'survival')]), method = 'spearman')
naucorr <- corrplot(corrM, order = 'original', type = 'lower', diag = FALSE)

################### Correlation abundance and area ###################

nautilid_data <- genus_data[genus_data$is.nautilid == TRUE, ]
ammonoid_data <- genus_data[genus_data$is.nautilid == FALSE, ]

ammaa <- ggplot(data = ammonoid_data, aes(x = log.area, y = sqrt.abun)) +
  geom_point(color = darkblue) + 
  labs(
    x = TeX('Geographic range ($log_{10}$ $km^2$)'),
    y = TeX('Abundance (square-root)')) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.y = element_text(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font),
  )
nauaa <- ggplot(data = nautilid_data, aes(x = log.area, y = sqrt.abun)) +
  geom_point(color = orange) + 
  labs(
    x = TeX('Geographic range ($log_{10}$ $km^2$)'),
    y = TeX('Abundance (square-root)')) +
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

cor(nautilid_data$log.area, nautilid_data$sqrt.abun, method = 'spearman')
cor(ammonoid_data$log.area, ammonoid_data$sqrt.abun, method = 'spearman')

################### Standardize predictor variables ###################

# From all predictors, subtract mean and divide by standard deviation
# (Z-score normalization)

genus_data$log.area.zscore <- zscore(genus_data$log.area)
genus_data$logvol.zscore <- zscore(genus_data$logvol)
genus_data$sqrt.abun.zscore <- zscore(genus_data$sqrt.abun)
genus_data$log.hatching.size.zscore <- zscore(genus_data$log.hatching.size)

# Convert boolean factors to numeric
genus_data$survival <- as.numeric(genus_data$survival)
genus_data$is.nautilid <- as.numeric(genus_data$is.nautilid)

# Try out a model without random effect and include area
# Hatching size has too many missing values
glm.test.area <- glm(
  survival ~ log.area.zscore + logvol.zscore + sqrt.abun.zscore,
  data = genus_data,
  family = binomial())

summary(glm.test.area)

# Try out a model without random effect to see if a random effect is necessary
glm.test <- glm(
  survival ~ logvol.zscore + sqrt.abun.zscore,
  data = genus_data,
  family = binomial())

summary(glm.test)

# Plot the residuals vs nautilid/ammonoid to see how much residual variation 
# it might explain
glm.test.resid <- rstandard(glm.test)

plot(
  glm.test.resid ~ as.factor(genus_data[!is.na(genus_data$logvol), ]$is.nautilid), 
    xlab = 'Nautilid',
    ylab = 'Standardized residuals')
abline(0, 0, lty = 2)

# is.nautilid does explain quite some residual variation. Distinguishing between
# nautilids and ammonoids makes sense.  

# Model fitting
# Having area in the model causes singularity, except if only area is included
# in the model, in which case the AIC is at 56.7, while otherwise it is at 30. 
genus_model <- glmer(
  survival ~ logvol.zscore + sqrt.abun.zscore + (1|is.nautilid), 
  data = genus_data,
  family = binomial())

summary(genus_model)
