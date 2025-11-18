# Get references for PBDB data to ward off adding duplicate info

setwd("~/Documents/Projects/Ammonites/code/")

library(dplyr)

cepha <- read.csv("../data/cephalopods.csv")
pbdb_maas <- read.csv("../data/pbdb_maastrichtian_danian_all.csv", , skip=17)

# add "PB" to all values in pbdb occurrence_no & collection_no to maintain 
# consistency with cepha
pbdb_maas <- pbdb_maas %>% mutate_at(c('occurrence_no', 'collection_no'), as.character)
pbdb_maas <- pbdb_maas %>% mutate(occurrence_no = paste('PB', occurrence_no, sep=''),
                                  collection_no = paste('PB', collection_no, sep=''))

# Get all the unique references where data in cepha comes from
pbdb_ref <- pbdb_maas %>% filter(occurrence_no %in% cepha$occurrence_no)
uniq_ref <- unique(pbdb_ref$primary_reference)

fileConn <- file("../data/pbdb_reflist.txt")
writeLines(uniq_ref, fileConn)
close(fileConn)

