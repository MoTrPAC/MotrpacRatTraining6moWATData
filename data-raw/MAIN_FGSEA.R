library(motrpacWATData)
library(motrpacWAT) # FGSEA functions
library(tidyverse)


## Transcriptomics ----------------------------------------------

table(grepl(";", TRNSCRPT_DA$MvF_SED$entrez_gene))
# 155 transcripts have multiple genes. Separate to multiple rows

TRNSCRPT_DA_SEP <- TRNSCRPT_DA %>%
  map(.f = ~ separate_rows(.x, entrez_gene, gene_symbol, sep = ";") %>%
        filter(entrez_gene != "NA"))

# Entrez to gene symbol conversion vector for leading edge
TRNSCRPT_entrez_to_symbol <- TRNSCRPT_DA_SEP$MvF_SED %>%
  select(entrez_gene, gene_symbol) %>%
  distinct() %>%
  deframe()

# Select gene sets that are largely unchanged when filtering
# to only those Entrez IDs present in the results.
TRNSCRPT_MSIGDB <- MSIGDB_PATHWAYS %>%
  mutate(entrez_gene = map(.x = entrez_gene, .f = intersect,
                           y = names(TRNSCRPT_entrez_to_symbol)),
         set_size_post = lengths(entrez_gene),
         ratio = set_size_post / set_size) %>%
  filter(ratio >= 0.85, # at least 85% of the original set remains
         set_size_post >= 15)

table(TRNSCRPT_MSIGDB$gs_subcat) # how many gene sets remain?
# GO:BP GO:CC GO:MF
#  2737   322   487

# FGSEA
TRNSCRPT_FGSEA <- map(TRNSCRPT_DA_SEP, function(res_i) {
  fgsea2(pathways = TRNSCRPT_MSIGDB,
         stats = get_ranking(res_i, genes = "entrez_gene"),
         seed = 0, nPermSimple = 10000,
         adjust.globally = TRUE, nproc = 1) %>%
    # Map Entrez IDs in leading edge subset to gene symbols
    mutate(leadingEdge_genes = map(
      .x = leadingEdge,
      .f = ~ na.omit(TRNSCRPT_entrez_to_symbol[as.character(.x)])
    )) %>%
    # Reorder columns
    select(pathway, gs_subcat, gs_description, everything()) %>%
    relocate(contrast, .before = leadingEdge)
})

# Save
usethis::use_data(TRNSCRPT_FGSEA, internal = FALSE,
                  overwrite = TRUE, version = 3,
                  compress = "bzip2")


## Proteomics --------------------------------------------------------

# Entrez to gene symbol conversion vector for leading edge
PROT_entrez_to_symbol <- fData(PROT_MSNSET) %>%
  select(entrez_gene, gene_symbol) %>%
  distinct() %>%
  deframe()

# Select gene sets that are largely unchanged when filtering
# to only those Entrez IDs present in the results.
PROT_MSIGDB <- MSIGDB_PATHWAYS %>%
  mutate(entrez_gene = map(.x = entrez_gene, .f = intersect,
                           y = names(PROT_entrez_to_symbol)),
         set_size_post = lengths(entrez_gene),
         ratio = set_size_post / set_size) %>%
  filter(ratio >= 0.85, # at least 85% of the original set remains
         set_size_post >= 15)

table(PROT_MSIGDB$gs_subcat) # how many gene sets remain?
# GO:BP GO:CC GO:MF
#   227    99    74

# FGSEA
PROT_FGSEA <- map(PROT_DA, function(res_i) {
  fgsea2(pathways = PROT_MSIGDB,
         stats = get_ranking(res_i, genes = "entrez_gene"),
         seed = 0, nPermSimple = 10000,
         adjust.globally = TRUE, nproc = 1) %>%
    # Map Entrez IDs in leading edge subset to gene symbols
    mutate(leadingEdge_genes = map(
      .x = leadingEdge,
      .f = ~ na.omit(PROT_entrez_to_symbol[as.character(.x)])
    )) %>%
    # Reorder columns
    select(pathway, gs_subcat, gs_description, everything()) %>%
    relocate(contrast, .before = leadingEdge)
})

# Save
usethis::use_data(PROT_FGSEA, internal = FALSE,
                  overwrite = TRUE, version = 3,
                  compress = "bzip2")


## Metabolomics -----------------------------------------------------

# Unlike with proteomics and transcriptomics, we are not
# limited to testing terms that are largely unchanged after
# filtering. This is because the RefMet chemical subclasses
# are homogenous groups (e.g., subsetting the acyl carnitines
# will still result in a group of just acyl carnitines).

# Reformat fData for use with fgsea2
sub_classes <- fData(METAB_MSNSET) %>%
  group_by(refmet_sub_class) %>%
  summarise(feature = list(feature_ID)) %>%
  mutate(gs_subcat = "refmet_subclass",
         gs_exact_source = refmet_sub_class,
         gs_description = refmet_sub_class,
         set_size = lengths(feature)) %>%
  filter(set_size >= 10)

nrow(sub_classes) # 19

# FGSEA
METAB_FGSEA <- map(METAB_DA, function(res_i) {
  fgsea2(pathways = sub_classes,
         gene_column = "feature",
         stats = get_ranking(res_i, genes = "feature"),
         seed = 0, nPermSimple = 10000,
         adjust.globally = TRUE, nproc = 1) %>%
    # Reorder columns
    dplyr::rename(refmet_sub_class = pathway) %>%
    select(-starts_with("gs_")) %>%
    relocate(contrast, .before = leadingEdge)
})

# Save
usethis::use_data(METAB_FGSEA, internal = FALSE,
                  overwrite = TRUE, version = 3,
                  compress = "bzip2")


## Save xlsx worksheets to inst/supp-tables ----
all_FGSEA <- list(TRNSCRPT_GO    = TRNSCRPT_FGSEA,
                  PROT_GO        = PROT_FGSEA,
                  PROT_MITOCARTA = PROT_MITOCARTA_FGSEA, # MitoCarta FGSEA
                  PHOSPHO_KSEA   = PHOSPHO_KSEA,         # KSEA
                  METAB_REFMET   = METAB_FGSEA) %>%
  list_transpose() %>% # nice function
  # Flatten leadingEdge columns to work with write_xlsx
  map_depth(.depth = 2, mutate, across(contains("leadingEdge"),
                                       ~ map_chr(.x, paste, collapse = ", ")))

paths <- file.path("inst", "supp-tables",
                   sprintf("FGSEA_%s.xlsx", names(all_FGSEA)))

map2(.x = all_FGSEA, .y = paths, .f = writexl::write_xlsx, .progress = TRUE)

