# Plot hatching size vs survival
setwd('~/Documents/Projects/Ammonites/code/')

library(dplyr)
library(ggplot2)
library(ggtext)
library(forcats)

# Specify plotting aesthetics
darkblue  <- '#12255A'
orange    <- '#E38C59'
errbarcol <- '#CA282C'
font <- 'Palatino'
fsize <- 9

# Read in data on hatching size and survival
ammon_hatch <- read.csv('../data/ammonoids_embryonic_shell_size.csv')
ammon_surv <- read.csv('../data/ammonoids_extinction_genus.csv')

nauti_hatch <- read.csv('../data/nautilids_embryonic_shell_size.csv')
nauti_surv <- read.csv('../data/nautilids_extinction_genus.csv')

# Simplify for further processing
ammon_surv <- select(ammon_surv, c('genus', 'survival'))
nauti_surv <- select(nauti_surv, c('genus', 'survival'))

# Merge with survival
ammon_hatch_median <- ammon_hatch %>%
  group_by(genus) %>%
  summarise(median.size = median(hatching.size..mm.))
ammon_hatch_surv <- merge(ammon_hatch_median, ammon_surv) %>%
  filter(!is.na(survival))

nauti_hatch_median <- nauti_hatch %>%
  group_by(genus) %>%
  summarise(median.size = median(hatching.size..mm.))
nauti_hatch_surv <- merge(nauti_hatch_median, nauti_surv) %>%
  filter(!is.na(survival))

pamhatch <- ggplot(data = ammon_hatch_surv, aes(x = survival, y = median.size)) +
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
pnahatch <- ggplot(data = nauti_hatch_surv, aes(x = survival, y = median.size)) +
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

wilcox.test(median.size ~ survival, data = ammon_hatch_surv)
wilcox.test(median.size ~ survival, data = nauti_hatch_surv)
# 
# 
# # Extinct at end-Cretaceous
# extinct <- survival[survival$survival==FALSE, 'genus'] 
# # Extinct before end of Maastrichtian
# already <- survival[is.na(survival$survival), 'genus']
# 
# # Use numerical for plotting
# survival[is.na(survival$survival), 'survival'] <- 3
# survival[survival$survival == TRUE, 'survival'] <- 2
# survival[survival$survival == FALSE, 'survival'] <- 1
# hatch <- merge(hatch, survival, by = 'genus')
# 
# # Need different genus name rotation for ammonoids/nautiloids
# if (cephas == 'ammonoids') {
#   rot <- 45
# } else {
#   rot <- 15
# }
# 
# p <- hatch %>%
#   mutate(genus = fct_reorder(genus, survival)) %>%
#   ggplot(aes(x = genus, y = hatching.size..mm.)) +
#   geom_jitter(position = position_jitter(0.05), cex = 3, shape = 1) +
#   labs(x = 'Genus', y = 'Hatching size (mm)') +
#   scale_x_discrete(labels = ~ if_else(
#     .x %in% extinct, paste0("<span style='color: red3'>", .x, "</span>"), 
#     if_else(
#       .x %in% already, paste0("<span style='color: blue4'>", .x, "</span>"), .x)
#   )) +
#   theme_classic() +
#   theme(
#     axis.text.x = element_markdown(size = 8, face = 'italic',  colour = 'black',
#                                    angle = rot, hjust = 0.9),
#     axis.text.y = element_text(size = 8, colour = 'black'),
#     axis.title = element_text(size = 8)
#     )
# 
# plotfile <- paste(cephas, '_hatching_sizes.pdf', sep = '')
# ggsave(paste('../results/hatching_size/', plotfile, sep = ''),
#        width = 12, height = 8, units = 'cm', plot = p)
# 
# 
# wilcox.test(hatching.size..mm. ~ survival, data = hatch)
