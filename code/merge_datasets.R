# Merge Flannery-Sutherland et al. (2024) and PBDB nautiloid and ammonoid
# datasets
setwd("~/Documents/Projects/Ammonites/code/")

# WARNING: Observed Nautilus campbelli in dataset. It appears to be a synonym of
# Eutrephoceras campbelli. Are there any true members of the genus Nautilus in 
# the Cretaceous? (Ward et al (2015) seem to think so)

library(dplyr)
library(stringr)
library(fossilbrush)

# Read in fossil occurrence data
ammon <- read.csv('../data/ESM_LateCretAmmon_2024/data/occs.csv')
#nauti <- read.csv('../data/pbdb_data_nautiloids.csv', skip=17)
pbdb_data <- read.csv('../data/pbdb_campanian_danian_all.csv', skip = 17)

# Read in data on which ammonoid genera were extant at the end of the
# Maastrichtian
end_maas <- read.csv('../data/ammonoids_end_maastrichtian.csv')

################### Data filtering, correct taxonomy ###################

# Filter out all entries that are not ammonoid or nautilid
pbdb_data_cepha <- pbdb_data %>% filter(class == 'Cephalopoda')
# Only accept these orders
cornu_orders = c('Ammonitida', 'Nautilida', 'Ammonoidea', 'Phylloceratida')
cornu_new <- pbdb_data_cepha %>% filter(order %in% cornu_orders)

# Update ages using fossilbrush
cornu_chrono <- chrono_scale(cornu_new, tscale = 'GTS_2020', srt = 'early_interval', 
  end = 'late_interval', max_ma = 'max_ma', min_ma = 'min_ma', verbose = FALSE)
cornu_new$max_ma <- cornu_chrono$newFAD
cornu_new$min_ma <- cornu_chrono$newLAD

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

# Make new columns for classification (nauti only)
#nauti <- nauti_clean %>% mutate(phylum      = 'Mollusca',
#                                class       = 'Cephalopoda',
#                                order       = 'Nautilida',
#                                suborder    =  NA,
#                                superfamily =  NA,
#                                family      = 'Nautilidae')

################### Merging main dataset ###################

# Restrict datasets to columns actually needed
# use same columns as Flannery-Sutherland et al. (2024). Remove region and 
# globe. I don't need them for this analysis. Also remove palaeo-positions. 
# I will rotate again & can use those to check

ammon <- select(ammon, 
  -c(region, globe, p_lng, p_lat))
cornu <- select(cornu, 
  -c(record_type, reid_no, flags, identified_name, identified_rank, 
    identified_no, difference, accepted_name, accepted_rank, accepted_no,
    early_interval, late_interval, reference_no, primary_reference
    ))

# add "PB" to all values in cornu occurrence_no & collection_no to maintain 
# consistency with ammon
cornu <- cornu %>% mutate_at(c('occurrence_no', 'collection_no'), as.character)
cornu <- cornu %>% mutate(occurrence_no = paste('PB', occurrence_no, sep=''),
                          collection_no = paste('PB', collection_no, sep=''))

# Remove cornu entries already in Flannery-Sutherland et al
cornu <- cornu %>% filter(!occurrence_no %in% ammon$occurrence_no)

# TO DO: Ammonite dataset contains occurrences with order "NA" inquire/remove
# Why are there specimen with unknown order?
ammon <- ammon %>% filter(!is.na(order)) 

# Merge datasets
cepha <- bind_rows(ammon, cornu)

################### Clean and disambiguate taxonomy ###################

genus_taxon <- cepha %>%
  select(c('order', 'suborder', 'superfamily', 'family', 'genus')) %>%
  distinct() %>%
  arrange(genus)

genus_taxon_ambiguous <- genus_taxon %>%
  count(genus) %>%
  filter(n > 1)

# For many genera, the problem is not different taxonomy but missing values
# Automatically correct for those by filling in missing data 

revise_taxonomy <- c()
taxon_levels <- c('order', 'suborder', 'superfamily', 'family', 'genus')

for (genus in genus_taxon_ambiguous$genus) {
  continue <- TRUE # Continue with this genus until a true taxonomic ambiguity 
  # is found, then move on to next genus
  for (taxlev in taxon_levels) {
    if (continue == TRUE) {
      clade <- unique(cepha[cepha$genus == genus, taxlev])
      # What if there is no ambiguity at this taxonomic level? Skip
      if (length(clade) > 1) {
        # In case there is ambiguity at this taxonomic level, check if it is 
        # missing values or true ambiguity
        if (sum(is.na(clade) == FALSE) > 1) {
          # In this case, there is more than one name for the taxonomic level
          # and the taxonomy has to be resolved by hand
          revise_taxonomy <- c(revise_taxonomy, genus)
          continue <- FALSE
        } else {
          # Replace all NAs with the already known taxonomic classification
          cepha[cepha$genus == genus, taxlev] <- clade[!is.na(clade)]
        }
      }
    }
  }
}

# TODO To avoid being dependent on conflicting taxonomy, which I don't want to 
# resolve right now, make a column is.nautilid TRUE/FALSE to bypass the issue 
cepha <- cepha %>%
  mutate(is.nautilid = ifelse(order == 'Nautilida', TRUE, FALSE))

# Remove the genus Conchorhynchus. These are nautiloid jaws.
# Also remove the genus Nautilus. Its presence is controversial. 
cepha <- filter(cepha, !genus %in% c('Conchorhynchus', 'Nautilus'))

# Reassign all Jeletzkytes occurrences to Hoploscaphites. These are 
# synonyms (Landman et al., 2010, B Am Mus Nat Hist)
cepha <- cepha %>%
  mutate(genus = ifelse(genus == 'Jeletzkytes', 'Hoploscaphites', genus))

################### Create Campanian dataset & save ################### 

# Get all occurrences from the Campanian

campa <- filter(cepha, max_ma <= 84 & min_ma >= 72)

write.csv(campa, file = '../data/campanian_cephalopods.csv', row.names = FALSE)

################### Create Maastrichtian data set & save ###################

# Remove all species/occurrences that do not lie around 66 mya
# (date of impact). Exclude all occurrences that are not in the Maastrichtian
cepha_maas <- filter(cepha, max_ma <= 72.5 & min_ma >= 65)

# Remove all ammonoid genera not extant at the end of the Maastrichtian
ext_ammon <- end_maas %>%
  filter(extant == TRUE) %>%
  select(genus)
# add in nautilids (these were extant)
ext_nauti <- unique(cepha_maas[cepha_maas$order == 'Nautilida', 'genus'])
extant <- c(ext_ammon$genus, ext_nauti)
# Restrict data to those present at the end of the Maastrichtian
cepha_end <- filter(cepha_maas, genus %in% extant)

# Save the end-Maastrichtian cephalopods (main analysis)
write.csv(cepha_end, file = '../data/cephalopods.csv', row.names = FALSE)

# Save the all-Maastrichtian cephalopod data
write.csv(cepha_maas, file = '../data/cephalopods_maastrichtian.csv', 
          row.names = FALSE)
