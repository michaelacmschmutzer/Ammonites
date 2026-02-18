# Filter paleorotated data so that only those occurrences are retained that are
# classified to the species level
setwd("~/Documents/Projects/Ammonites/code/")

library(dplyr)

cepha_rota <- read.csv('../data/cephalopods_palaeorotated.csv')

# Which ammonoid species existed at the end of the Maastrichtian? 
end_maas_ammon <- read.csv('../data/ammonoids_end_maastrichtian_species.csv')

# Only retain those occurrences that can be allocated to a species
cepha_rota_clean <- cepha_rota %>%
  filter(!is.na(species)) %>%
  # For convenience, to have a unique identifier for every species
  mutate(species = paste(genus, species)) %>%
  select(-c(genus))

# As on the genus level, assume all nautilid species are extant at the end
# of the Maastrichtian
ext_ammon <- end_maas_ammon %>%
  mutate(species = paste(genus, species)) %>%
  filter(extant == TRUE) %>%
  select(species)
# add in nautilids (these were extant)
ext_nauti <-
  unique(cepha_rota_clean[cepha_rota_clean$order == 'Nautilida', 'species'])
extant <- c(ext_ammon$species, ext_nauti)
# Restrict data to those present at the end of the Maastrichtian
cepha_rota_end <- filter(cepha_rota_clean, species %in% extant)

write.csv(cepha_rota_end, '../data/cephalopods_palaeorotated_species.csv')