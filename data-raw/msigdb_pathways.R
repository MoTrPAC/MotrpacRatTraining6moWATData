library(motrpacWAT)
library(dplyr)

# Get Gene Ontology terms from MSigDB
MSIGDB_PATHWAYS <- msigdbr2(species = "Rattus norvegicus",
                            genes = "entrez_gene",
                            gs_subcat = c("GO:BP", "GO:MF", "GO:CC"),
                            capitalize = TRUE) %>%
  mutate(set_size = lengths(entrez_gene)) %>%
  filter(set_size >= 15, set_size <= 300)
# nrow(MSIGDB_PATHWAYS) # 5693

# During FGSEA, sets are further filtered based on a membership ratio
# and are limited to sets with no more than 300 genes

# Save
usethis::use_data(MSIGDB_PATHWAYS, internal = FALSE,
                  overwrite = TRUE, version = 3,
                  compress = "bzip2")

