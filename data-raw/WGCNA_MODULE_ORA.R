library(MotrpacRatTraining6moWATData)
library(MotrpacRatTraining6moWAT)
library(tidyverse)

## Proteomics ----------------------------------------------------------------
# Map Entrez IDs to gene symbols
PROT_entrez_to_symbol <- pluck(PROT_WGCNA, "modules") %>%
  select(entrez_gene, gene_symbol) %>%
  distinct() %>%
  deframe()

# Genes by module
PROT_MOD_LIST <- pluck(PROT_WGCNA, "modules") %>%
  filter(moduleColor != "grey") %>%
  group_by(moduleID) %>%
  summarise(feature = list(unique(entrez_gene))) %>%
  deframe()

# Number of unique genes in each module
lengths(PROT_MOD_LIST)
#   P1   P2   P3   P4   P5   P6   P7   P8   P9  P10  P11
# 3860 1403 1391  728  667  434  414  233  227  191  159

# all genes
PROT_UNIVERSE <- unique(unlist(PROT_MOD_LIST))

## Select gene sets that are largely unchanged when filtering
## to only those Entrez IDs present in the results.
PROT_MSIGDB <- MSIGDB_PATHWAYS %>%
  mutate(entrez_gene = map(entrez_gene, intersect, y = PROT_UNIVERSE),
         set_size_post = lengths(entrez_gene),
         ratio = set_size_post / set_size) %>%
  filter(ratio >= 0.85, # at least 85% of the original set remains
         set_size_post >= 15)

table(PROT_MSIGDB$gs_subcat)
# GO:BP GO:CC GO:MF
#   230    97    77

## ORA ----
PROT_MODULE_ORA <- fora2(pathways = PROT_MSIGDB,
                         genes = PROT_MOD_LIST,
                         universe = PROT_UNIVERSE,
                         adjust.method = "scale",
                         adjust.globally = TRUE) %>%
  mutate(overlapGenes_symbol = map(
    .x = overlapGenes,
    .f = ~ na.omit(PROT_entrez_to_symbol[as.character(.x)])
  )) %>%
  relocate(module, gs_subcat, pathway, gs_description) %>%
  relocate(overlap, maxOverlap, overlapRatio, .after = size)

# Save
usethis::use_data(PROT_MODULE_ORA, internal = FALSE, overwrite = TRUE,
                  version = 3, compress = "bzip2")

## Transcriptomics -----------------------------------------------------------
TRNSCRPT_MODULES <- pluck(TRNSCRPT_WGCNA, "modules") %>%
  separate_rows(entrez_gene, gene_symbol, sep = ";") %>%
  filter(entrez_gene != "NA",
         moduleColor != "grey")

# Entrez to gene symbol conversion vector
TRNSCRPT_entrez_to_symbol <- TRNSCRPT_MODULES %>%
  select(entrez_gene, gene_symbol) %>%
  distinct() %>%
  deframe()

# Genes by module
TRNSCRPT_MOD_LIST <- TRNSCRPT_MODULES %>%
  group_by(moduleID) %>%
  summarise(feature = list(unique(entrez_gene))) %>%
  deframe()

# Number of unique genes in each module
lengths(TRNSCRPT_MOD_LIST)
#   T1   T2   T3   T4   T5   T6   T7   T8   T9  T10  T11  T12  T13  T14
# 4425 3061 2014 1435  528  418  314  221  202  153  112   97   49   32

# all genes
TRNSCRPT_UNIVERSE <- unique(unlist(TRNSCRPT_MOD_LIST))

## Select gene sets that are largely unchanged when filtering
## to only those Entrez IDs present in the results.
TRNSCRPT_MSIGDB <- MSIGDB_PATHWAYS %>%
  mutate(entrez_gene = map(entrez_gene, intersect, y = TRNSCRPT_UNIVERSE),
         set_size_post = lengths(entrez_gene),
         ratio = set_size_post / set_size) %>%
  filter(ratio >= 0.85, # at least 85% of the original set remains
         set_size_post >= 15)

table(TRNSCRPT_MSIGDB$gs_subcat)
# GO:BP GO:CC GO:MF
#   927   110   158

## ORA ----
TRNSCRPT_MODULE_ORA <- fora2(pathways = TRNSCRPT_MSIGDB,
                             genes = TRNSCRPT_MOD_LIST,
                             universe = TRNSCRPT_UNIVERSE,
                             adjust.method = "scale",
                             adjust.globally = TRUE) %>%
  mutate(overlapGenes_symbol = map(
    .x = overlapGenes,
    .f = ~ na.omit(TRNSCRPT_entrez_to_symbol[as.character(.x)])
  )) %>%
  relocate(module, gs_subcat, pathway, gs_description) %>%
  relocate(overlap, maxOverlap, overlapRatio, .after = size)

# Save
usethis::use_data(TRNSCRPT_MODULE_ORA, internal = FALSE, overwrite = TRUE,
                  version = 3, compress = "bzip2")

## Metabolomics ----------------------------------------------
# List of features by module
METAB_MOD_LIST <- pluck(METAB_WGCNA, "modules") %>%
  filter(moduleColor != "grey") %>%
  group_by(moduleID) %>%
  summarise(feature = list(unique(feature_ID))) %>%
  deframe()

lengths(METAB_MOD_LIST)
#  M1  M2  M3  M4  M5  M6  M7
# 415 221 137  99  86  69  30

# all features (excluding those 6 from M0)
METAB_UNIVERSE <- unique(unlist(METAB_MOD_LIST))

# Unlike with proteomics and transcriptomics, we are not
# limited to testing terms that are largely unchanged after
# filtering. This is because the RefMet chemical subclasses
# are homogenous groups (e.g., subsetting the acyl carnitines
# will still result in a group of just acyl carnitines).

# Reformat fData for use with fora2
REFMET_SUBCLASSES <- pluck(METAB_WGCNA, "modules") %>%
  filter(moduleColor != "grey") %>%
  group_by(refmet_sub_class) %>%
  summarise(feature = list(feature_ID)) %>%
  mutate(gs_subcat = "refmet_sub_class",
         gs_exact_source = refmet_sub_class,
         gs_description = refmet_sub_class,
         set_size = lengths(feature)) %>%
  filter(set_size >= 10)

length(REFMET_SUBCLASSES) # 6

## ORA ----
METAB_MODULE_ORA <- fora2(pathways = REFMET_SUBCLASSES,
                          genes = METAB_MOD_LIST,
                          gene_column = "feature",
                          universe = METAB_UNIVERSE,
                          adjust.method = "scale",
                          adjust.globally = TRUE) %>%
  relocate(module, gs_subcat, pathway, gs_description) %>%
  relocate(overlap, maxOverlap, overlapRatio, .after = size) %>%
  dplyr::rename(refmet_sub_class = pathway,
                overlapMetabolites = overlapGenes) %>%
  select(-starts_with("gs_"))

# Save
usethis::use_data(METAB_MODULE_ORA, internal = FALSE, overwrite = TRUE,
                  version = 3, compress = "bzip2")

