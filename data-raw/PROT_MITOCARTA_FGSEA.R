library(MotrpacRatTraining6moData)
library(MotrpacRatTraining6moWATData)
library(MotrpacRatTraining6moWAT)
library(tidyverse)
library(fgsea)

# Entrez to gene symbol conversion vector for leading edge
entrez_to_symbol <- pluck(PROT_DA, "MvF_SED") %>%
  filter(!is.na(entrez_gene)) %>%
  dplyr::select(entrez_gene, gene_symbol) %>%
  deframe()

# Human to rat gene conversion
human_to_rat <- RAT_TO_HUMAN_GENE %>%
  dplyr::select(HUMAN_ORTHOLOG_SYMBOL, RAT_NCBI_GENE_ID) %>%
  distinct() %>%
  deframe()
# two genes map to two entrez genes each - leave as-is

# Convert human gene symbols to rat Entrez IDs
PROT_MITOCARTA <- MITOCARTA_HS %>%
  mutate(rat_entrez = map(human_genes,
                          ~ as.character(na.omit(human_to_rat[.x]))),
         set_size = lengths(rat_entrez),
         rat_entrez = map(rat_entrez, intersect, names(entrez_to_symbol)),
         set_size_post = lengths(rat_entrez),
         ratio = set_size_post / set_size) %>%
  filter(ratio >= 0.85,
         set_size <= 300, # Same as MSIGDB_PATHWAYS filter
         set_size_post >= 5)
# 4 more pathways than before:
# [1] "Amino acid metabolism"   "CI assembly factors"
# [3] "OXPHOS assembly factors" "mt-rRNA modifications"

# These pathways previously had ratios just below 0.85 because set_size
# was incorrectly based on human genes

# List of 68 pathways to test
MITOCARTA_PATHWAYS <- PROT_MITOCARTA %>%
  select(pathway, rat_entrez) %>%
  tibble::deframe()

## FGSEA
PROT_MITOCARTA_FGSEA <- map(PROT_DA, function(res_i) {
  stats <- rank_genes(res_i, genes = "entrez_gene")

  map(names(stats), function(contrast_i) {
    message(contrast_i)
    set.seed(0)
    fgseaMultilevel(pathways = MITOCARTA_PATHWAYS,
                    stats = stats[[contrast_i]],
                    nPermSimple = 10000, nproc = 1) %>%
      mutate(contrast = contrast_i)
  }) %>%
    data.table::rbindlist() %>%
    mutate(contrast = factor(contrast, levels = unique(contrast)),
           padj = p.adjust(pval, method = "BH"),
           leadingEdge_genes = map(leadingEdge,
                                   ~ na.omit(entrez_to_symbol[.x])),
           leadingEdge = map(leadingEdge, as.numeric)) %>% # consistency
    left_join(dplyr::select(PROT_MITOCARTA, pathway, hierarchy),
              by = "pathway") %>%
    relocate(hierarchy, .after = pathway) %>%
    relocate(contrast, .before = leadingEdge)
})

# Save
usethis::use_data(PROT_MITOCARTA_FGSEA, internal = FALSE,
                  overwrite = TRUE, version = 3, compress = "bzip2")

