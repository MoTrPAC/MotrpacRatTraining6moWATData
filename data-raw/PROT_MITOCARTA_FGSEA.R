library(MotrpacRatTraining6moData)
library(MotrpacRatTraining6moWATData)
library(tidyverse)
library(fgsea)

# Entrez to gene symbol conversion vector for leading edge
entrez_to_symbol <- pluck(PROT_DA, "MvF_SED") %>%
  filter(!is.na(entrez_gene)) %>%
  dplyr::select(entrez_gene, gene_symbol) %>%
  deframe()

# Human to rat gene conversion
rat2human <- RAT_TO_HUMAN_GENE %>%
  dplyr::select(human_ortholog = HUMAN_ORTHOLOG_SYMBOL,
                entrez_gene = RAT_NCBI_GENE_ID) %>%
  distinct()

# Use human data from MitoCarta3.0
PROT_MITOCARTA <- file.path("data-raw", "Human.MitoCarta3.0.xls") %>%
  readxl::read_xls(sheet = "C MitoPathways") %>%
  .[, 2:4] %>%
  setNames(c("pathway", "hierarchy", "human_ortholog")) %>%
  separate_rows(human_ortholog, sep = ", ") %>%
  left_join(rat2human, by = "human_ortholog") %>%
  group_by(pathway) %>%
  mutate(num_genes = n()) %>%
  filter(entrez_gene %in% as.numeric(names(entrez_to_symbol)),
         num_genes <= 300) %>% # Same as MSIGDB_PATHWAYS size filter
  group_by(pathway) %>%
  mutate(num_genes_post = n(),
         ratio = num_genes_post / num_genes) %>%
  filter(ratio >= 0.85,
         num_genes_post >= 5) %>%
  group_by(pathway, hierarchy, num_genes,
           num_genes_post, ratio) %>%
  summarise(entrez_gene = list(entrez_gene)) %>%
  ungroup()

# List of 65 pathways to test
MITOCARTA_PATHWAYS <- PROT_MITOCARTA %>%
  select(pathway, entrez_gene) %>%
  tibble::deframe()

## FGSEA
PROT_MITOCARTA_FGSEA <- map(PROT_DA, function(res_i) {
  stats <- get_ranking(res_i, genes = "entrez_gene")

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

