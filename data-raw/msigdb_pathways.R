library(dplyr)

# Get Gene Ontology pathways from MSigDB
MSIGDB_PATHWAYS <- msigdbr2(species = "Rattus norvegicus",
                            genes = "entrez_gene",
                            gs_subcat = c("GO:BP", "GO:MF", "GO:CC"),
                            capitalize = TRUE) %>%
  mutate(set_size = lengths(entrez_gene)) %>%
  filter(set_size >= 15, set_size <= 500)
# nrow(MSIGDB_PATHWAYS) # 6552

# Save
usethis::use_data(MSIGDB_PATHWAYS, internal = FALSE,
                  overwrite = TRUE, version = 3,
                  compress = "bzip2")

