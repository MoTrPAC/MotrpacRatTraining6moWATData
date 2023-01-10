library(MotrpacRatTraining6moWAT)
library(dplyr)

# Get Gene Ontology terms from MSigDB
MSIGDB_PATHWAYS <- msigdbr2(species = "Rattus norvegicus",
                            genes = "entrez_gene",
                            gs_subcat = c("GO:BP", "GO:MF", "GO:CC"),
                            capitalize = TRUE) %>%
  mutate(set_size = lengths(entrez_gene)) %>%
  filter(set_size >= 15, set_size <= 300)

table(MSIGDB_PATHWAYS$gs_subcat)
# GO:BP GO:CC GO:MF
#  3913   437   756

# Save
usethis::use_data(MSIGDB_PATHWAYS, internal = FALSE,
                  overwrite = TRUE, version = 3,
                  compress = "bzip2")

