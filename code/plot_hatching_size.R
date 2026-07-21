# Plot hatching size vs survival
setwd('~/Documents/Projects/Ammonites/code/')

library(dplyr)
library(ggplot2)
library(ggtext)
library(forcats)
library(gridExtra)

# Specify plotting aesthetics
darkblue  <- '#12255A'
orange    <- '#E38C59'
errbarcol <- '#CA282C'
font <- 'Palatino'
fsize <- 9

###################  Read in data  ################### 

# Read in data on hatching size, survival, and taxonomy
ammon_hatch <- read.csv('../data/ammonoids_embryonic_shell_size.csv')
ammon_surv_genus <- read.csv('../data/ammonoids_extinction_genus.csv')
ammon_surv_spp <- read.csv('../data/ammonoids_extinction_species.csv')
cepha <- read.csv('../data/cephalopods.csv')

nauti_hatch <- read.csv('../data/nautilids_embryonic_shell_size.csv')
nauti_surv_genus <- read.csv('../data/nautilids_extinction_genus.csv')
nauti_surv_spp <- read.csv('../data/nautilids_extinction_species.csv')

# Read in data on species presence/absence end Maastrichtian
ammon_extant_spp <- read.csv('../data/ammonoids_end_maastrichtian_species.csv')
ammon_extant_spp <- ammon_extant_spp %>% mutate(species = paste(genus, species))
ammon_extant_spp <- select(ammon_extant_spp, c('species', 'extant'))
# Remove NA entries
ammon_extant_spp <- na.omit(ammon_extant_spp)

# Paste together binomial name for easier processing
ammon_surv_spp <- ammon_surv_spp %>% mutate(species = paste(genus, species))
nauti_surv_spp <- nauti_surv_spp %>% mutate(species = paste(genus, species))
ammon_hatch <- ammon_hatch %>% mutate(species = paste(genus, species))

# Simplify for further processing
ammon_surv_genus <- select(ammon_surv_genus, c('genus', 'survival'))
nauti_surv_genus <- select(nauti_surv_genus, c('genus', 'survival'))
ammon_surv_spp <- select(ammon_surv_spp, c('species', 'survival'))
nauti_surv_spp <- select(nauti_surv_spp, c('species', 'survival'))

###################  Genus-specific preparations  ################### 
# Take median to look at genus-level
ammon_hatch_median <- ammon_hatch %>%
  group_by(genus) %>%
  summarise(median.size = median(hatching.size..mm.))
nauti_hatch_median <- nauti_hatch %>%
  group_by(genus) %>%
  summarise(median.size = median(hatching.size..mm.))

# Merge with survival
ammon_hatch_surv_genus <- merge(ammon_hatch_median, ammon_surv_genus) %>%
  filter(!is.na(survival))
nauti_hatch_surv_genus <- merge(nauti_hatch_median, nauti_surv_genus) %>%
  filter(!is.na(survival))

###################  Species-specific preparations  ################### 

# Restrict species to those present at end of Maastrichtian
extant <- ammon_extant_spp[ammon_extant_spp$extant == TRUE, 'species']
ammon_hatch_spp <- ammon_hatch %>%
  filter(species %in% extant)

# Merge with survival
ammon_hatch_surv_spp <- merge(ammon_hatch_spp, ammon_surv_spp) 

###################  Sort ammonoids by suborder  ################### 

unique_genus <- unique(ammon_hatch_surv_genus$genus)
ammon_taxon <- cepha %>%
  filter(genus %in% unique_genus) %>%
  filter(is.na(suborder) == FALSE) %>%
  filter(is.na(family) == FALSE) %>%
  group_by(genus) %>%
  summarise(suborder = unique(suborder),
            family = unique(family))

ammon_hatch <- merge(ammon_hatch, ammon_taxon)
# add in genus-level survival
ammon_hatch <- merge(ammon_hatch, ammon_surv_genus)

###################  Plotting and analysis genus  ################### 

pamhatch <- ggplot(data = ammon_hatch_surv_genus, 
  aes(x = survival, y = median.size)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, color = darkblue) +
  labs(x = 'Ammonoids', y = 'Median hatching size (mm)') +
  scale_x_discrete(labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
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
pnahatch <- ggplot(data = nauti_hatch_surv_genus,
  aes(x = survival, y = median.size)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, color = orange) +
  labs(x = 'Nautilids', y = 'Median hatching size (mm)') +
  scale_x_discrete(labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
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
p <- grid.arrange(pamhatch, pnahatch, ncol = 2)
ggsave('../results/hatching_size/Genus_hatching_survival.png',
       width = 12, height = 6, units = 'cm', dpi = 600, plot = p)

wilcox.test(median.size ~ survival, data = ammon_hatch_surv_genus)
wilcox.test(median.size ~ survival, data = nauti_hatch_surv_genus)

###################  Plotting individual hatching sizes  ################### 

# Need different genus name rotation for ammonoids/nautiloids
if (cephas == 'ammonoids') {
  rot <- 45
} else {
  rot <- 15
}

ammon_hatch <- ammon_hatch %>%
  arrange(survival, suborder, family) %>%
  mutate(genus = factor(genus, levels = unique(genus))) 

  ggplot(data = ammon_hatch, aes(x = hatching.size..mm., y = genus, 
    color = survival)) +
  geom_jitter(position = position_jitter(0.05), cex = 3, shape = 1) +
  labs(x = 'Hatching size (mm)', y = 'Genus') 
  
  
  scale_x_discrete(labels = ~ if_else(
    .x %in% extinct, paste0("<span style='color: red3'>", .x, "</span>"),
    if_else(
      .x %in% already, paste0("<span style='color: blue4'>", .x, "</span>"), .x)
  )) +
  theme_classic() +
  theme(
    axis.text.x = element_markdown(size = 8, face = 'italic',  colour = 'black',
                                   angle = rot, hjust = 0.9),
    axis.text.y = element_text(size = 8, colour = 'black'),
    axis.title = element_text(size = 8)
    )

# plotfile <- paste(cephas, '_hatching_sizes.pdf', sep = '')
# ggsave(paste('../results/hatching_size/', plotfile, sep = ''),
#        width = 12, height = 8, units = 'cm', plot = p)
# 
# 
# wilcox.test(hatching.size..mm. ~ survival, data = hatch)
