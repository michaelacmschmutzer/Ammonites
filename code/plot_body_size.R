# Plot body size vs survival
setwd('~/Documents/Projects/Ammonites/code/')

library(dplyr)
library(ggplot2)
library(gridExtra)
library(latex2exp)

# Specify plotting aesthetics
darkblue  <- '#12255A'
orange    <- '#E38C59'
errbarcol <- '#CA282C'
font <- 'Palatino'
fsize <- 9

cepha <- read.csv('../data/cephalopods.csv')

body_sizes <- read.csv(
  '../data/doi_10_5061_dryad_zpc866t5b__v20210211/genus.sizes.ranges.rev.csv')

# Read in survivorship data 
surv_nauti <- read.csv('../data/nautilids_extinction_genus.csv') 
surv_ammon <- read.csv('../data/ammonoids_extinction_genus.csv')
# Simplify
surv_nauti <- select(surv_nauti, c('genus', 'survival'))
surv_ammon <- select(surv_ammon, c('genus', 'survival'))

# Which genera are nautilids?
is_nautilid <- cepha %>%
  select(genus, order) %>%
  distinct(genus, .keep_all = TRUE) %>%
  filter(order == 'Nautilida')

# Get body sizes
body_sizes_ammon <- body_sizes %>%
  filter(genus %in% cepha$genus) %>%
  filter(!genus %in% is_nautilid$genus)
body_sizes_nauti <- body_sizes %>%
  filter(genus %in% cepha$genus) %>%
  filter(genus %in% is_nautilid$genus)

# Merge with survivorship
body_ammon_surv <- merge(body_sizes_ammon, surv_ammon)
body_nauti_surv <- merge(body_sizes_nauti, surv_nauti)

# Plot
pbam <- ggplot(data = body_ammon_surv, aes(x = survival, y = logvol)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, color = darkblue) +
  xlab('Ammonoids') +
  ylab(TeX('Body volume ($log_{10}$ $mm^3$)')) +
  scale_x_discrete(labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  ylim(0, 7) +
  stat_summary(
    fun = median, geom = 'point', shape = 18, size = 3.5, color = errbarcol) +
  stat_summary(
    fun.max = function(z) { quantile(z, 0.75) },
    fun.min = function(z) { quantile(z, 0.25) },
    geom = 'errorbar', color = errbarcol, width = 0.1) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.y = element_text(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font),
  )
pbnau <- ggplot(data = body_nauti_surv, aes(x = survival, y = logvol)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, color = orange) +
  xlab('Nautilids') +
  ylab(TeX('Body volume ($log_{10}$ $mm^3$)')) +
  scale_x_discrete(labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  ylim(0, 7) +
  stat_summary(
    fun = median, geom = 'point', shape = 18, size = 3.5, color = errbarcol) +
  stat_summary(
    fun.max = function(z) { quantile(z, 0.75) },
    fun.min = function(z) { quantile(z, 0.25) },
    geom = 'errorbar', color = errbarcol, width = 0.1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.y = element_text(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font),
  )
p <- grid.arrange(pbam, pbnau, ncol = 2)
ggsave('../results/body_size/Genus_body_size_survival.png', device = pdf,
      width = 12, height = 6, units = 'cm', dpi = 600, plot = p)

wilcox.test(logvol ~ survival, data = body_ammon_surv)
