# Analysis of abundance datasets
setwd("~/Documents/Projects/Ammonites/code/")

library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(scales)
library(stringr)
library(systemfonts)

source('statistical-functions.R')
source('geospatial-functions.R')

# Specify plotting aesthetics
darkblue  <- '#12255A'
orange    <- '#E38C59'
errbarcol <- '#CA282C'
font <- 'Palatino'
fsize <- 9

###################  Read in data ################### 

cepha_rota <- read.csv('../data/cephalopods_palaeorotated.csv')
maast <- read.csv('../data/cephalopods_maastrichtian.csv')
campa <- read.csv('../data/campanian_cephalopods.csv')

# Read in survivorship data 
surv_nauti <- read.csv('../data/nautilids_extinction_genus.csv') 
surv_ammon <- read.csv('../data/ammonoids_extinction_genus.csv')
# Simplify
surv_nauti <- select(surv_nauti, c('genus', 'survival'))
surv_ammon <- select(surv_ammon, c('genus', 'survival'))

################### Maastrichtian analyses ################### 

# Count raw abundance in dataset vs survival
tally_genus_nauti <- cepha_rota %>% 
  filter(order == 'Nautilida') %>%
  count(genus)
tally_genus_ammon <- cepha_rota %>% 
  filter(!order == 'Nautilida') %>%
  count(genus)

# Save the raw abundance data
write.csv(tally_genus_nauti, 
  '../results/abundance_survivalship/nautilids_abundance_raw.csv', 
  row.names = FALSE)
write.csv(tally_genus_ammon, 
  '../results/abundance_survivalship/ammonoids_abundance_raw.csv',
  row.names = FALSE)

# Add survivalship info
tally_nauti_surv <- merge(tally_genus_nauti, surv_nauti)
tally_ammon_surv <- merge(tally_genus_ammon, surv_ammon)

# Plotting
set.seed(36)
pammon <- ggplot(data = tally_ammon_surv, aes(x = survival, y = n)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, color = darkblue) +
  geom_text_repel(data = subset(tally_ammon_surv, survival == TRUE),
    aes(label = genus),
    nudge_x = - 0.1,
    family = font,
    size = 2.5,
    fontface = 'italic',
    box.padding   = unit(0.5, 'lines'),
    point.padding = unit(0.6, 'lines')) +
  labs(x = 'Ammonoids', y = 'Number of occurrences per genus') +
  scale_x_discrete(labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.y = element_text(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font),
  )
pnauti <- ggplot(data = tally_nauti_surv, aes(x = survival, y = n)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, colour = orange) +
  geom_text_repel(data = subset(tally_nauti_surv, n > 100),
    aes(label = genus),
    family = font,
    size = 2.5,
    fontface = 'italic',
    box.padding   = unit(0.35, 'lines'),
    point.padding = unit(0.3, 'lines')) +
  labs(x = 'Nautilids', y = 'Number of occurrences per genus') +
  scale_x_discrete(labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.y = element_text(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font),
  )
p <- grid.arrange(pammon, pnauti, ncol = 2)
ggsave('../results/abundance_survivalship/Genus_abundance_survival.png',
      width = 12, height = 6, units = 'cm', dpi = 600, plot = p)


wilcox.test(n ~ survival, data = tally_ammon_surv)
wilcox.test(n ~ survival, data = tally_nauti_surv)

# Plot the relative abundance of each genus in locations with multiple
# genera. Express relative abundance as a percentage of total number
# of occurrences at that location. Ignore all locations with occurrences from
# a single genus
# Probably want to put this through an equal-area filter too to synonymise
# highly similar localities
abd_ammon <- cepha_rota %>% 
  filter(!order == 'Nautilida') %>%
  group_by(lat, lng, genus) %>%
  count(genus)
abd_ammon$coords <- paste(abd_ammon$lat, abd_ammon$lng) # easier to handle
# Get the number of genera per fossil locality
ngenus_ammon <- abd_ammon %>%
  group_by(coords) %>%
  summarise(ngenus = n_distinct(genus))
# Which localities have more than one genus?
multiloc_ammon <- ngenus_ammon %>%
  filter(ngenus > 1) %>%
  select(coords)
# Restrict analysis to those localities that have more than one genus
abd_ammon <- abd_ammon %>%
  filter(coords %in% multiloc_ammon$coords)
# Calculate relative abundance of genus as a percentage
# Total number of occurrences at locus
loc_total_ammon <- abd_ammon %>%
  group_by(coords) %>%
  summarise(total.occ = sum(n))
abd_ammon <- merge(abd_ammon, loc_total_ammon)
# Calculate the percentages
abd_ammon <- abd_ammon %>%
  mutate(percent.occ = 100*n/total.occ)
# Get the median percentage occupancy per genus
abd_ammon_genus <- abd_ammon %>%
  group_by(genus) %>%
  summarise(med.percent.occ = median(percent.occ))
# Add survivorship
abd_ammon_genus <- merge(abd_ammon_genus, surv_ammon)

# Same for nautilids
abd_nauti <- cepha_rota %>% 
  filter(order == 'Nautilida') %>%
  group_by(lat, lng, genus) %>%
  count(genus)
abd_nauti$coords <- paste(abd_nauti$lat, abd_nauti$lng) # easier to handle
# Get the number of genera per fossil locality
ngenus_nauti <- abd_nauti %>%
  group_by(coords) %>%
  summarise(ngenus = n_distinct(genus))
# Which localities have more than one genus?
multiloc_nauti <- ngenus_nauti %>%
  filter(ngenus > 1) %>%
  select(coords)
# Restrict analysis to those localities that have more than one genus
abd_nauti <- abd_nauti %>%
  filter(coords %in% multiloc_nauti$coords)
# Calculate relative abundance of genus as a percentage
# Total number of occurrences at locus
loc_total_nauti <- abd_nauti %>%
  group_by(coords) %>%
  summarise(total.occ = sum(n))
abd_nauti <- merge(abd_nauti, loc_total_nauti)
abd_nauti <- abd_nauti %>%
  mutate(percent.occ = 100*n/total.occ)
# Get the median percentage occupancy per genus
abd_nauti_genus <- abd_nauti %>%
  group_by(genus) %>%
  summarise(med.percent.occ = median(percent.occ))
# Add survivorship
abd_nauti_genus <- merge(abd_nauti_genus, surv_nauti)

# Using this method, the relative abundance of dominating genera at localities
# with lots of different genera will be smaller than at localities where few
# genera have been sampled. This makes a direct comparison difficult.
# It is easier to be a big fish in a small pond. 
pabdam <- ggplot(data = abd_ammon_genus, 
  aes(x = survival, y = med.percent.occ)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, color = darkblue) +
  labs(x = 'Ammonoids', y = 'Median abundance \n at fossil site (percent)') +
  scale_x_discrete(labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  ylim(0, 100) +
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
pabdna <- ggplot(data = abd_nauti_genus, 
  aes(x = survival, y = med.percent.occ)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, color = orange) +
  labs(x = 'Nautilids', y = 'Median abundance \n at fossil site (percent)') +
  scale_x_discrete(labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  ylim(0, 100) +
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
p <- grid.arrange(pabdam, pabdna, ncol = 2)
ggsave('../results/abundance_survivalship/Genus_locality_survival.png',
       width = 12, height = 6, units = 'cm', dpi = 600, plot = p)

  
median(abd_ammon_genus[abd_ammon_genus$survival == FALSE, 'med.percent.occ'])
median(abd_ammon_genus[abd_ammon_genus$survival ==  TRUE, 'med.percent.occ'])
median(abd_nauti_genus[abd_nauti_genus$survival == FALSE, 'med.percent.occ'])
median(abd_nauti_genus[abd_nauti_genus$survival ==  TRUE, 'med.percent.occ'])

wilcox.test(med.percent.occ ~ survival, data = abd_ammon_genus)
wilcox.test(med.percent.occ ~ survival, data = abd_nauti_genus)

# Look at total number of fossil localities per genus and survival
tally_loc_ammon <- cepha_rota %>% 
  filter(!order == 'Nautilida') %>%
  group_by(genus) %>%
  summarise(n.loc = n_distinct(lat, lng))
loc_ammon <- merge(tally_loc_ammon, surv_ammon)

tally_loc_nauti <- cepha_rota %>% 
  filter(order == 'Nautilida') %>%
  group_by(genus) %>%
  summarise(n.loc = n_distinct(lat, lng))
loc_nauti <- merge(tally_loc_nauti, surv_nauti)

plocam <- ggplot(data = loc_ammon, aes(x = survival, y = n.loc)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, color = darkblue) +
  labs(x = 'Ammonoids', y = 'Number of unique sites') +
  scale_x_discrete(labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  ylim(0, max(loc_ammon$n.loc) + 10) +
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
plocna <- ggplot(data = loc_nauti, aes(x = survival, y = n.loc)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, color = orange) +
  labs(x = 'Nautilids', y = 'Number of unique sites') +
  scale_x_discrete(labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  ylim(0, max(loc_ammon$n.loc) + 10) +
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
p <- grid.arrange(plocam, plocna, ncol = 2)
ggsave('../results/abundance_survivalship/Genus_numloc_survival.png',
       width = 12, height = 6, units = 'cm', dpi = 600, plot = p)

# TODO Get the number of unique grid cells per genus
# Use PALEOMAP for now
numcell <- distinct(cepha_rota, genus, order)
numcell$n.cell <- 0

for (g in unique(cepha_rota$genus)) {
  genus.df <- cepha_rota %>%
    filter(genus == g)
  cell.df <- unique.within.grid(genus.df, 'p_lat_PALEOMAP', 'p_lng_PALEOMAP')
  n <- nrow(cell.df)
  numcell[numcell$genus == g, 'n.cell'] <- n
}

# Split into nautilids and ammonoids
ncell_ammon <- numcell[numcell$order != 'Nautilida', ]
ncell_nauti <- numcell[numcell$order == 'Nautilida', ]

# Add survival information
ncell_ammon_surv <- merge(ncell_ammon, surv_ammon)
ncell_nauti_surv <- merge(ncell_nauti, surv_nauti)

pcellam <- ggplot(data = ncell_ammon_surv, aes(x = survival, y = n.cell)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, color = darkblue) +
  labs(x = 'Ammonoids', y = 'Number of grid cells') +
  scale_x_discrete(labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  ylim(0, max(ncell_ammon_surv$n.cell) + 10) +
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
pcellna <- ggplot(data = ncell_nauti_surv, aes(x = survival, y = n.cell)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, color = orange) +
  labs(x = 'Nautilids', y = 'Number of grid cells') +
  scale_x_discrete(labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  ylim(0, max(ncell_ammon_surv$n.cell) + 10) +
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
p <- grid.arrange(pcellam, pcellna, ncol = 2)
ggsave('../results/abundance_survivalship/Genus_numcell_survival.png',
       width = 12, height = 6, units = 'cm', dpi = 600, plot = p)

################### Comparison with Campanian ################### 

# Count raw abundance in Maastrichtian (all cephalopods, also those ammonoids
# not extant at the end of the Maastrichtian)
tally_genus_ammon_maas <- maast %>% 
  filter(!order == 'Nautilida') %>%
  count(genus)
tally_genus_nauti_maas <- maast %>% 
  filter(order == 'Nautilida') %>%
  count(genus)

# Count raw abundance in Campanian
tally_genus_ammon_camp <- campa %>% 
  filter(!order == 'Nautilida') %>%
  count(genus)
tally_genus_nauti_camp <- campa %>% 
  filter(order == 'Nautilida') %>%
  count(genus)

# Compare raw abundances in the Campanian and the Maastrichtian
comp_abun_ammon <- merge(tally_genus_ammon_maas, tally_genus_ammon_camp, 
  by = 'genus')
names(comp_abun_ammon) <- c('genus', 'n.maas', 'n.camp')
comp_abun_ammon[is.na(comp_abun_ammon$n.camp), 'n.camp'] <- 0

comp_abun_nauti <- merge(tally_genus_nauti_maas, tally_genus_nauti_camp, 
  by = 'genus')
names(comp_abun_nauti) <- c('genus', 'n.maas', 'n.camp')
comp_abun_nauti[is.na(comp_abun_nauti$n.camp), 'n.camp'] <- 0

# Add survivorship data for ammonoids
survivors <- surv_ammon[surv_ammon$survival == TRUE, 'genus']
comp_abun_ammon <- comp_abun_ammon %>%
  mutate(survival = ifelse(genus %in% survivors, TRUE, FALSE))

# Add information on which ammonoids were present at the end of the 
# Maastrichtian
end_maas <- unique(cepha_rota[cepha_rota$order != 'Nautilida', 'genus'])
comp_abun_ammon <- comp_abun_ammon %>%
  mutate(end.maas = ifelse(genus %in% end_maas, TRUE, FALSE))

# Plot correlations for both ammonoids and nautilids
pamcor <- ggplot(data = comp_abun_ammon, aes(x = n.camp, y = n.maas, 
    alpha = end.maas)) +
  geom_point(cex = 3, color = darkblue) +
  geom_abline(slope = 1, intercept = 0, linetype = 'twodash',
    color = errbarcol) +
  geom_text_repel(
    data = subset(comp_abun_ammon, survival == TRUE),
    aes(label = genus),
    family = font,
    size = 2.5,
    fontface = 'italic',
    force = 3,
    box.padding = unit(0.5, 'lines'),
    point.padding = unit(0.2, 'lines')
  ) +
  scale_alpha_manual(values = c(0.6, 1.0)) +
  labs(
    x = 'Number of Campanian occurrences',
    y = 'Number of Maastrichtian occurrences'
  ) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.y = element_text(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font),
    legend.position = 'none'
  )
pnacor <- ggplot(data = comp_abun_nauti, aes(x = n.camp, y = n.maas)) +
  geom_point(cex = 3, color = orange) +
  geom_abline(slope = 1, intercept = 0, linetype = 'twodash',
    color = errbarcol) +
  geom_text_repel(
    data = subset(comp_abun_nauti, n.maas >= 20),
    aes(label = genus),
    family = font,
    size = 2.5,
    fontface = 'italic',
    box.padding = unit(0.5, 'lines'),
    point.padding = unit(0.2, 'lines')
  ) +
  labs(
    x = 'Number of Campanian occurrences',
    y = 'Number of Maastrichtian occurrences'
  ) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  expand_limits(y = 1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.y = element_text(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font)
  )
p <- grid.arrange(pamcor, pnacor, ncol = 2)
ggsave('../results/abundance_survivalship/Genus_abundance_comparison.png',
       width = 12, height = 6, units = 'cm', dpi = 600, plot = p)

cor.test(comp_abun_ammon$n.camp, comp_abun_ammon$n.maas, method = 'spearman')
cor.test(comp_abun_nauti$n.camp, comp_abun_nauti$n.maas, method = 'spearman')

comp_abun_ammon <- comp_abun_ammon %>%
  mutate(change.n = n.maas - n.camp)

comp_end_maas <- comp_abun_ammon %>%
  filter(end.maas == TRUE)

wilcox.test(change.n ~ survival, data = comp_abun_ammon)
wilcox.test(change.n ~ survival, data = comp_end_maas)
