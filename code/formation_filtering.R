# Quick script to make formation data human-modifiable
setwd("~/Documents/Projects/Ammonites/code/")

library(dplyr)
library(stringr)

# Read in raw PBDB file
pbdb_data <-read.csv('../data/formation_info/pbdb_data.csv', skip = 18)
flsu_data <- read.csv('../data/ESM_LateCretAmmon_2024/data/occs.csv')

# Reduce to ammonites from Maastrichtian only
pbdb_data_cepha <- pbdb_data %>% filter(class == 'Cephalopoda')
pbdb_data_cepha <- pbdb_data %>% filter(early_interval != 'Danian')
cornu_orders = c('Ammonitida', 'Ammonoidea', 'Phylloceratida')
cornu_new <- pbdb_data_cepha %>% filter(order %in% cornu_orders)

# Some basic cleaning
# use entries with "accepted_rank" either genus or species
cornu_filter <- cornu_new %>% filter(accepted_rank %in% c('genus', 'species'))
# Some entries have subgenera/no longer accepted genera in brackets. Remove
# the info between brackets (e.g. 'Pachydiscus (Pachydiscus) gollevillensis')
cornu_bino <- cornu_filter %>%
  mutate(accepted_name = str_remove(accepted_name, '\\([:alpha:]*\\)\\s'))
# split species and genus names. Mark missing species names "NA"
cornu_split <- cornu_bino %>% mutate(
  genus = if_else(accepted_rank=='genus',
                  true = accepted_name,
                  false = sapply(str_split(accepted_name, " "), '[', 1)),
  species = if_else(accepted_rank=='species',
                    true = sapply(str_split(accepted_name, " "), '[', 2),
                    false = NA))

# N.B. identified_name often has "." or "?" as marker of uncertainty. 
# Use only genus name (option: leave out entirely). 
cornu <- cornu_split %>% mutate(
  species = replace(species, 
                    str_detect(identified_name, pattern = '[[:punct:]]'), NA))

# Simplify a bit
cornu <- select(
  cornu, 
  -c(record_type, reid_no, flags, identified_name, identified_rank, 
  identified_no, difference, accepted_name, accepted_rank, accepted_no,
  reference_no
))

# Filer out all Flannery-Sutherland occurrences
# add "PB" to all values in cornu occurrence_no & collection_no to maintain 
# consistency with flsu
cornu <- cornu %>% mutate_at(c('occurrence_no', 'collection_no'), as.character)
cornu <- cornu %>% mutate(occurrence_no = paste('PB', occurrence_no, sep=''),
                          collection_no = paste('PB', collection_no, sep=''))

# Remove cornu entries already in Flannery-Sutherland et al
cornu <- cornu %>% filter(!occurrence_no %in% flsu_data$occurrence_no)

# Save
write.csv(cornu, '../data/formation_info/cephalopods_formation_info.csv',
  row.names = FALSE)
