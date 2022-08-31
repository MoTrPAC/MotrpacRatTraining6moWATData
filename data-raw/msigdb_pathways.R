library(dplyr)

# Get pathways from MSigDB
msigdb_pathways <- msigdbr2(species = "rat",
                            genes = "entrez_gene",
                            gs_subcat = c("GO:BP", "GO:MF", "GO:CC",
                                          "CP:REACTOME", "CP:KEGG"),
                            capitalize = TRUE) %>%
  mutate(num_genes = lengths(entrez_gene)) %>%
  filter(num_genes >= 15, num_genes <= 500)
# nrow(msigdb_pathways) # 6552

# Save
usethis::use_data(msigdb_pathways, internal = FALSE,
                  overwrite = TRUE, version = 3,
                  compress = "bzip2")
