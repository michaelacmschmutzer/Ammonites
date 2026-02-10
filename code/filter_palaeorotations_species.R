# Filter paleorotated data so that only those occurrences are retained that are
# classified to the species level
setwd("~/Documents/Projects/Ammonites/code/")

library(dplyr)

cepha_rota <- read.csv('../data/cephalopods_palaeorotated.csv')
cepha_rota_clean <- cepha_rota %>%
  filter(!is.na(species)) %>%
  mutate(species = paste(genus, species)) %>%
  select(-c(genus))

write.csv(cepha_rota_clean, '../data/cephalopods_palaeorotated_species.csv')