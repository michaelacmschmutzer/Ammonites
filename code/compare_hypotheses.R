# Compare strength of evidence for different hypotheses for all datasets
setwd("~/Documents/Projects/Ammonites/code/")

library(dplyr)
library(ggplot2)
library(forcats)

source('statistical-functions.R')

# Specify plotting aesthetics
darkblue  <- '#12255A'
orange    <- '#E38C59'
errbarcol <- '#CA282C'
font <- 'Palatino'
fsize <- 9

################### Compile hypothesis testing data set ###################
cepha <- read.csv('../data/cephalopods.csv')

# Genus-level data
body_sizes <- read.csv(
  '../data/doi_10_5061_dryad_zpc866t5b__v20210211/genus.sizes.ranges.rev.csv')
geo_nauti_genus <- read.csv(
  '../results/genus/geographic_distributions/nautilids_distributions_genus.csv')
geo_ammon_genus <- read.csv(
  '../results/genus/geographic_distributions/ammonoids_distributions_genus.csv')
geo_nauti_boot_genus <- read.csv(
  '../results/genus/subsampling_distributions/nautilids/bootstrap.csv')
geo_ammon_boot_genus <- read.csv(
  '../results/genus/subsampling_distributions/ammonoids/bootstrap.csv')
geo_nauti_jack_genus <- read.csv(
  '../results/genus/subsampling_distributions/nautilids/jackknife.csv')
geo_ammon_jack_genus <- read.csv(
  '../results/genus/subsampling_distributions/ammonoids/jackknife.csv')

hatch_nauti <- read.csv('../data/nautilids_embryonic_shell_size.csv')
hatch_ammon <- read.csv('../data/ammonoids_embryonic_shell_size.csv')
surv_nauti_genus <- read.csv('../data/nautilids_extinction_genus.csv')
surv_ammon_genus <- read.csv('../data/ammonoids_extinction_genus.csv')

# Ammonoid genera
ammon_genus <- unique(cepha[cepha$is.nautilid == FALSE, 'genus'])
# Nautilid genera
nauti_genus <- unique(cepha[cepha$is.nautilid == TRUE, 'genus'])

# Species-level data
geo_nauti_species <- read.csv(
  '../results/species/geographic_distributions/nautilids_distributions_species.csv')
geo_ammon_species <- read.csv(
  '../results/species/geographic_distributions/ammonoids_distributions_species.csv')
geo_nauti_boot_species <- read.csv(
  '../results/species/subsampling_distributions/nautilids/bootstrap.csv')
geo_ammon_boot_species <- read.csv(
  '../results/species/subsampling_distributions/ammonoids/bootstrap.csv')
geo_nauti_jack_species <- read.csv(
  '../results/species/subsampling_distributions/nautilids/jackknife.csv')
geo_ammon_jack_species <- read.csv(
  '../results/species/subsampling_distributions/ammonoids/jackknife.csv')

surv_nauti_species <- read.csv('../data/nautilids_extinction_species.csv')
surv_ammon_species <- read.csv('../data/ammonoids_extinction_species.csv')

# Place genus and species name in single column
surv_nauti_species <- surv_nauti_species %>%
  mutate(species = paste(genus, species))
surv_ammon_species <- surv_ammon_species %>%
  mutate(species = paste(genus, species))

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
surv_nauti_species <- dplyr::select(surv_nauti_species, c('species', 'survival'))
surv_ammon_species <- dplyr::select(surv_ammon_species, c('species', 'survival'))

# Simplify sub-sampling results
# Use PALEOMAP for now
boot_nauti_genus <- geo_nauti_boot_genus %>%
  filter(model == 'PALEOMAP') %>%
  group_by(genus) %>%
  summarise(boot.median.area = median(area.km2))
boot_ammon_genus <- geo_ammon_boot_genus %>%
  filter(model == 'PALEOMAP') %>%
  group_by(genus) %>%
  summarise(boot.median.area = median(area.km2))
jack_nauti_genus <- geo_nauti_jack_genus %>%
  filter(model == 'PALEOMAP') %>%
  group_by(genus) %>%
  summarise(jack.median.area = median(area.km2))
jack_ammon_genus <- geo_ammon_jack_genus %>%
  filter(model == 'PALEOMAP') %>%
  group_by(genus) %>%
  summarise(jack.median.area = median(area.km2))

boot_nauti_species <- geo_nauti_boot_species %>%
  filter(model == 'PALEOMAP') %>%
  group_by(species) %>%
  summarise(boot.median.area = median(area.km2))
boot_ammon_species <- geo_ammon_boot_species %>%
  filter(model == 'PALEOMAP') %>%
  group_by(species) %>%
  summarise(boot.median.area = median(area.km2))
jack_nauti_species <- geo_nauti_jack_species %>%
  filter(model == 'PALEOMAP') %>%
  group_by(species) %>%
  summarise(jack.median.area = median(area.km2))
jack_ammon_species <- geo_ammon_jack_species %>%
  filter(model == 'PALEOMAP') %>%
  group_by(species) %>%
  summarise(jack.median.area = median(area.km2))

# Handle body sizes
bodyvol_nauti <- body_sizes %>% filter(genus %in% nauti_genus)
bodyvol_ammon <- body_sizes %>% filter(genus %in% ammon_genus)

# Merge with survival data
hatch_nauti <- na.omit(merge(hatch_nauti_genus, surv_nauti_genus))
hatch_ammon <- na.omit(merge(hatch_ammon_genus, surv_ammon_genus))

bodyvol_nauti <- merge(bodyvol_nauti, surv_nauti_genus)
bodyvol_ammon <- merge(bodyvol_ammon, surv_ammon_genus)

geo_nauti_genus <- merge(geo_nauti_genus, surv_nauti_genus)
geo_ammon_genus <- merge(geo_ammon_genus, surv_ammon_genus)

boot_nauti_genus <- merge(boot_nauti_genus, surv_nauti_genus)
boot_ammon_genus <- merge(boot_ammon_genus, surv_ammon_genus)

jack_nauti_genus <- merge(jack_nauti_genus, surv_nauti_genus)
jack_ammon_genus <- merge(jack_ammon_genus, surv_ammon_genus)

geo_nauti_species <- merge(geo_nauti_species, surv_nauti_species)
geo_ammon_species <- merge(geo_ammon_species, surv_ammon_species)

boot_nauti_species <- merge(boot_nauti_species, surv_nauti_species)
boot_ammon_species <- merge(boot_ammon_species, surv_ammon_species)

jack_nauti_species <- merge(jack_nauti_species, surv_nauti_species)
jack_ammon_species <- merge(jack_ammon_species, surv_ammon_species)

################### Calculating effect sizes ################### 

# Get p-value
wt_hatch <- wilcox.test(med.hatching.size ~ factor(survival),
  data = hatch_ammon, exact = TRUE)

surv <- hatch_ammon[hatch_ammon$survival == TRUE, 'med.hatching.size']
exti <- hatch_ammon[hatch_ammon$survival == FALSE, 'med.hatching.size']

surv <- hatch_nauti[hatch_nauti$survival == TRUE, 'med.hatching.size']
exti <- hatch_nauti[hatch_nauti$survival == FALSE, 'med.hatching.size']

# Get common language effect size (U / (n1 * n2)) with confidence interval
variables <- c(
  'Hatching size (genus)',
  'Body size (genus)',
  'Geographic range (genus)',
  'Geographic range (species)',
  'Geographic range (bootstrap, genus)',
  'Geographic range (bootstrap, species)',
  'Geographic range (jackknife, genus)',
  'Geographic range (jackknife, species)'
)

# Create dataframe for effect sizes
effect.sizes.ammon <- data.frame(
  variable = variables, eff.size = 0, ci.lower = 0, ci.upper = 0, 
  num.ext = 0, num.surv = 0, var.ext = 0, var.surv = 0)
effect.sizes.nauti <- data.frame(
  variable = variables, eff.size = 0, ci.lower = 0, ci.upper = 0,
  num.ext = 0, num.surv = 0, var.ext = 0, var.surv = 0)

# Fill in effect sizes ... ammonoids
effect.sizes.ammon[effect.sizes.ammon$variable == 'Hatching size (genus)', 2:4] <-
  wilmanwhit.effect.size(hatch_ammon, 'med.hatching.size')
effect.sizes.ammon[effect.sizes.ammon$variable == 'Body size (genus)', 2:4] <- 
  wilmanwhit.effect.size(bodyvol_ammon, 'logvol')
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (genus)', 2:4] <- 
  wilmanwhit.effect.size(geo_ammon_genus, 'PALEOMAP.area.km2')
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (bootstrap, genus)', 2:4] <- 
  wilmanwhit.effect.size(boot_ammon_genus, 'boot.median.area')
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (jackknife, genus)', 2:4] <- 
  wilmanwhit.effect.size(jack_ammon_genus, 'jack.median.area')
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (species)', 2:4] <- 
  wilmanwhit.effect.size(geo_ammon_species, 'PALEOMAP.area.km2')
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (bootstrap, species)', 2:4] <- 
  wilmanwhit.effect.size(boot_ammon_species, 'boot.median.area')
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (jackknife, species)', 2:4] <- 
  wilmanwhit.effect.size(jack_ammon_species, 'jack.median.area')

# Fill in number of extinct taxa ... ammonoids
effect.sizes.ammon[effect.sizes.ammon$variable == 'Hatching size (genus)', 5] <-
  sum(hatch_ammon$survival == FALSE)
effect.sizes.ammon[effect.sizes.ammon$variable == 'Body size (genus)', 5] <- 
  sum(bodyvol_ammon$survival == FALSE)
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (genus)', 5] <- 
  sum(geo_ammon_genus$survival == FALSE)
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (bootstrap, genus)', 5] <- 
  sum(boot_ammon_genus$survival == FALSE)
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (jackknife, genus)', 5] <- 
  sum(jack_ammon_genus$survival == FALSE)
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (species)', 5] <- 
  sum(geo_ammon_species$survival == FALSE)
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (bootstrap, species)', 5] <- 
  sum(boot_ammon_species$survival == FALSE)
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (jackknife, species)', 5] <- 
  sum(jack_ammon_species$survival == FALSE)

# Fill in number of surviving taxa ... ammonoids
effect.sizes.ammon[effect.sizes.ammon$variable == 'Hatching size (genus)', 6] <-
  sum(hatch_ammon$survival == TRUE)
effect.sizes.ammon[effect.sizes.ammon$variable == 'Body size (genus)', 6] <- 
  sum(bodyvol_ammon$survival == TRUE)
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (genus)', 6] <- 
  sum(geo_ammon_genus$survival == TRUE)
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (bootstrap, genus)', 6] <- 
  sum(boot_ammon_genus$survival == TRUE)
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (jackknife, genus)', 6] <- 
  sum(jack_ammon_genus$survival == TRUE)
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (species)', 6] <- 
  sum(geo_ammon_species$survival == TRUE)
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (bootstrap, species)', 6] <- 
  sum(boot_ammon_species$survival == TRUE)
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (jackknife, species)', 6] <- 
  sum(jack_ammon_species$survival == TRUE)

# Fill in variance for extinct taxa ... ammonoids
effect.sizes.ammon[effect.sizes.ammon$variable == 'Hatching size (genus)', 7] <-
  var(hatch_ammon[hatch_ammon$survival == FALSE, 'med.hatching.size'])
effect.sizes.ammon[effect.sizes.ammon$variable == 'Body size (genus)', 7] <- 
  var(bodyvol_ammon[bodyvol_ammon$survival == FALSE, 'logvol'])
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (genus)', 7] <- 
  var(geo_ammon_genus[geo_ammon_genus$survival == FALSE, 'PALEOMAP.area.km2'])
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (bootstrap, genus)', 7] <- 
  var(boot_ammon_genus[boot_ammon_genus$survival == FALSE, 'boot.median.area'])
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (jackknife, genus)', 7] <- 
  var(jack_ammon_genus[jack_ammon_genus$survival == FALSE, 'jack.median.area'])
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (species)', 7] <- 
  var(geo_ammon_species[geo_ammon_species$survival == FALSE, 'PALEOMAP.area.km2'])
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (bootstrap, species)', 7] <- 
  var(boot_ammon_species[boot_ammon_species$survival == FALSE, 'boot.median.area'])
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (jackknife, species)', 7] <- 
  var(jack_ammon_species[jack_ammon_species$survival == FALSE, 'jack.median.area'])

# Fill in variance for surviving taxa ... ammonoids
effect.sizes.ammon[effect.sizes.ammon$variable == 'Hatching size (genus)', 8] <-
  var(hatch_ammon[hatch_ammon$survival == TRUE, 'med.hatching.size'])
effect.sizes.ammon[effect.sizes.ammon$variable == 'Body size (genus)', 8] <- 
  var(bodyvol_ammon[bodyvol_ammon$survival == TRUE, 'logvol'])
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (genus)', 8] <- 
  var(geo_ammon_genus[geo_ammon_genus$survival == TRUE, 'PALEOMAP.area.km2'])
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (bootstrap, genus)', 8] <- 
  var(boot_ammon_genus[boot_ammon_genus$survival == TRUE, 'boot.median.area'])
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (jackknife, genus)', 8] <- 
  var(jack_ammon_genus[jack_ammon_genus$survival == TRUE, 'jack.median.area'])
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (species)', 8] <- 
  var(geo_ammon_species[geo_ammon_species$survival == TRUE, 'PALEOMAP.area.km2'])
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (bootstrap, species)', 8] <- 
  var(boot_ammon_species[boot_ammon_species$survival == TRUE, 'boot.median.area'])
effect.sizes.ammon[effect.sizes.ammon$variable == 'Geographic range (jackknife, species)', 8] <- 
  var(jack_ammon_species[jack_ammon_species$survival == TRUE, 'jack.median.area'])

# Fill in effect sizes ... nautilids
effect.sizes.nauti[effect.sizes.nauti$variable == 'Hatching size (genus)', 2:4] <-
  wilmanwhit.effect.size(hatch_nauti, 'med.hatching.size')
effect.sizes.nauti[effect.sizes.nauti$variable == 'Body size (genus)', 2:4] <- 
  wilmanwhit.effect.size(bodyvol_nauti, 'logvol')
effect.sizes.nauti[effect.sizes.nauti$variable == 'Geographic range (genus)', 2:4] <- 
  wilmanwhit.effect.size(geo_nauti_genus, 'PALEOMAP.area.km2')
effect.sizes.nauti[effect.sizes.nauti$variable == 'Geographic range (bootstrap, genus)', 2:4] <- 
  wilmanwhit.effect.size(boot_nauti_genus, 'boot.median.area')
effect.sizes.nauti[effect.sizes.nauti$variable == 'Geographic range (jackknife, genus)', 2:4] <- 
  wilmanwhit.effect.size(jack_nauti_genus, 'jack.median.area')
# effect.sizes.nauti[effect.sizes.nauti$variable == 'Geographic range (species)', 2:4] <- 
#   wilmanwhit.effect.size(geo_nauti_species, 'PALEOMAP.area.km2')
# effect.sizes.nauti[effect.sizes.nauti$variable == 'Geographic range (bootstrap, species)', 2:4] <- 
#   wilmanwhit.effect.size(boot_nauti_species, 'boot.median.area')
# effect.sizes.nauti[effect.sizes.nauti$variable == 'Geographic range (jackknife, species)', 2:4] <- 
#   wilmanwhit.effect.size(jack_nauti_species, 'jack.median.area')

# Fill in number of extinct taxa ... nautilids
effect.sizes.nauti[effect.sizes.nauti$variable == 'Hatching size (genus)', 5] <-
  sum(hatch_nauti$survival == FALSE)
effect.sizes.nauti[effect.sizes.nauti$variable == 'Body size (genus)', 5] <- 
  sum(bodyvol_nauti$survival == FALSE)
effect.sizes.nauti[effect.sizes.nauti$variable == 'Geographic range (genus)', 5] <- 
  sum(geo_nauti_genus$survival == FALSE)
effect.sizes.nauti[effect.sizes.nauti$variable == 'Geographic range (bootstrap, genus)', 5] <- 
  sum(boot_nauti_genus$survival == FALSE)
effect.sizes.nauti[effect.sizes.nauti$variable == 'Geographic range (jackknife, genus)', 5] <- 
  sum(jack_nauti_genus$survival == FALSE)

# Fill in number of surviving taxa ... nautilids
effect.sizes.nauti[effect.sizes.nauti$variable == 'Hatching size (genus)', 6] <-
  sum(hatch_nauti$survival == TRUE)
effect.sizes.nauti[effect.sizes.nauti$variable == 'Body size (genus)', 6] <- 
  sum(bodyvol_nauti$survival == TRUE)
effect.sizes.nauti[effect.sizes.nauti$variable == 'Geographic range (genus)', 6] <- 
  sum(geo_nauti_genus$survival == TRUE)
effect.sizes.nauti[effect.sizes.nauti$variable == 'Geographic range (bootstrap, genus)', 6] <- 
  sum(boot_nauti_genus$survival == TRUE)
effect.sizes.nauti[effect.sizes.nauti$variable == 'Geographic range (jackknife, genus)', 6] <- 
  sum(jack_nauti_genus$survival == TRUE)

# Fill in variance for extinct taxa ... nautilids
effect.sizes.nauti[effect.sizes.nauti$variable == 'Hatching size (genus)', 7] <-
  var(hatch_nauti[hatch_nauti$survival == FALSE, 'med.hatching.size'])
effect.sizes.nauti[effect.sizes.nauti$variable == 'Body size (genus)', 7] <- 
  var(bodyvol_nauti[bodyvol_nauti$survival == FALSE, 'logvol'])
effect.sizes.nauti[effect.sizes.nauti$variable == 'Geographic range (genus)', 7] <- 
  var(geo_nauti_genus[geo_nauti_genus$survival == FALSE, 'PALEOMAP.area.km2'])
effect.sizes.nauti[effect.sizes.nauti$variable == 'Geographic range (bootstrap, genus)', 7] <- 
  var(boot_nauti_genus[boot_nauti_genus$survival == FALSE, 'boot.median.area'])
effect.sizes.nauti[effect.sizes.nauti$variable == 'Geographic range (jackknife, genus)', 7] <- 
  var(jack_nauti_genus[jack_nauti_genus$survival == FALSE, 'jack.median.area'])

# Fill in variance for surviving taxa ... nautilids
effect.sizes.nauti[effect.sizes.nauti$variable == 'Hatching size (genus)', 8] <-
  var(hatch_nauti[hatch_nauti$survival == TRUE, 'med.hatching.size'])
effect.sizes.nauti[effect.sizes.nauti$variable == 'Body size (genus)', 8] <- 
  var(bodyvol_nauti[bodyvol_nauti$survival == TRUE, 'logvol'])
effect.sizes.nauti[effect.sizes.nauti$variable == 'Geographic range (genus)', 8] <- 
  var(geo_nauti_genus[geo_nauti_genus$survival == TRUE, 'PALEOMAP.area.km2'])
effect.sizes.nauti[effect.sizes.nauti$variable == 'Geographic range (bootstrap, genus)', 8] <- 
  var(boot_nauti_genus[boot_nauti_genus$survival == TRUE, 'boot.median.area'])
effect.sizes.nauti[effect.sizes.nauti$variable == 'Geographic range (jackknife, genus)', 8] <- 
  var(jack_nauti_genus[jack_nauti_genus$survival == TRUE, 'jack.median.area'])

# Remove missing rows from nautilid effect sizes
effect.sizes.nauti <- effect.sizes.nauti %>%
  filter(eff.size != 0)

################### Save effect sizes ################### 

write.csv(effect.sizes.ammon, 
  '../results/comparing_hypotheses/ammonoids_effect_sizes.csv', 
  row.names = FALSE)

write.csv(effect.sizes.nauti, 
  '../results/comparing_hypotheses/nautilids_effect_sizes.csv', 
  row.names = FALSE)

################### Plotting effect sizes ################### 

# Make sure variables are ordered as specified in data frame row order
effect.sizes.ammon$variable <- fct_rev(factor(effect.sizes.ammon$variable, 
  levels = effect.sizes.ammon$variable))
effect.sizes.nauti$variable <- fct_rev(factor(effect.sizes.nauti$variable,
  levels = effect.sizes.nauti$variable))

# Plotting ... ammonoids
efam <- ggplot(effect.sizes.ammon, aes(x = eff.size, y = variable)) +
  geom_hline(yintercept = nrow(effect.sizes.ammon) + 0.6, linewidth = 1) +
  geom_vline(xintercept = 0.5, linetype = 'longdash', linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci.lower, xmax = ci.upper, width = 0.3),
    color = 'black') +
  geom_point(size = 3, color = darkblue) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1), 
    expand = c(0,0)) +
  labs(x = 'Effect size', y = '') +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.x = element_text(size = fsize, colour = 'black', family = font),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(t = 2, r = 10, b = 0, l = 0, unit = 'pt')
  )
ggsave(
  '../results/comparing_hypotheses/Effect_sizes_ammonoids.png',
  width = 12, height = 6, units = 'cm', dpi = 600, plot = efam)

# And nautilids
efna <- ggplot(effect.sizes.nauti, aes(x = eff.size, y = variable)) +
  geom_hline(yintercept = nrow(effect.sizes.nauti) + 0.6, linewidth = 1) +
  geom_vline(xintercept = 0.5, linetype = 'longdash', linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci.lower, xmax = ci.upper, width = 0.3),
                color = 'black') +
  geom_point(size = 3, color = orange) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1), 
                     expand = c(0,0)) +
  labs(x = 'Effect size', y = '') +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.x = element_text(size = fsize, colour = 'black', family = font),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(t = 2, r = 10, b = 0, l = 0, unit = 'pt')
  )
ggsave(
  '../results/comparing_hypotheses/Effect_sizes_nautilids.png',
  width = 12, height = 6, units = 'cm', dpi = 600, plot = efna)
