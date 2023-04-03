library(MotrpacRatTraining6moWATData)
library(MotrpacRatTraining6moWAT)
library(tidyverse)
library(data.table)
library(fgsea)

# Reformat DEA results ----
human_res <- PHOSPHO_DA %>%
  map(function(res_i) {
    filter(res_i, !is.na(human_uniprot)) %>%
      mutate(num_sites = str_count(human_site, ";") + 1) %>%
      separate_rows(human_site) %>% # single-site-level data
      mutate(human_feature = paste0(human_uniprot, "_", human_site)) %>%
      select(contrast, human_feature, logFC, P.Value, num_sites) %>%
      distinct()
  })

# List of substrate sites by kinase (379 kinases before filtering)
KS_sets <- PSP_KINASE_SUBSTRATE %>%
  transmute(kinase = GENE,
            substrate = paste0(SUB_ACC_ID, "_", SUB_MOD_RSD)) %>%
  # Filter to what is in the DEA results
  filter(substrate %in% human_res$MvF_SED$human_feature) %>%
  group_by(kinase) %>%
  summarise(substrate = list(substrate)) %>%
  deframe()

# How many substrates are in PSP?
table(unique(human_res$MvF_SED$human_feature) %in% unlist(KS_sets))
# FALSE  TRUE
# 18019  1118

# Only about 5.8% of the substrates are in PSP

# How many will pass the size filter?
KS_sets <- KS_sets[lengths(KS_sets) >= 3]
length(KS_sets) # 121 kinases

# Removing small kinase sets drops some substrate sites
table(unique(human_res$MvF_SED$human_feature) %in% unlist(KS_sets))
# FALSE  TRUE
# 18063  1074

## Conversion vectors for KSEA leadingEdge
# Uniprot to gene symbol sites
human_uniprot_to_symbol <- PSP_KINASE_SUBSTRATE %>%
  transmute(across(c(SUB_ACC_ID, SUB_GENE),
                   ~ paste0(.x, "_", SUB_MOD_RSD))) %>%
  distinct() %>%
  deframe()

# Map human to rat sites (single-site level)
human_to_rat <- fData(PHOSPHO_EXP) %>%
  filter(!is.na(human_uniprot)) %>%
  separate_rows(site, human_site, sep = ";") %>%
  transmute(human = paste0(human_uniprot, "_", human_site),
            rat = paste0(gene_symbol, "_", site)) %>%
  distinct() %>%
  deframe()

# KSEA
PHOSPHO_KSEA <- map(human_res, function(res_i)
{
  # Ranking metric is weighted more heavily by singly-phosphorylated peptides
  # rank_list <- rank_genes(res_i, genes = "human_feature",
  #                          metric = "-log10(P.Value)*sign(logFC)/num_sites")

  # Ranking metric is weighted equally (default metric)
  rank_list <- rank_genes(res_i, genes = "human_feature")

  map(names(rank_list), function(contr_i) {
    set.seed(0)
    fgseaMultilevel(pathways = KS_sets,
                    stats = rank_list[[contr_i]],
                    minSize = 3,
                    nproc = 1, nPermSimple = 10000) %>%
      mutate(contrast = contr_i)
  }) %>%
    rbindlist() %>%
    mutate(padj = p.adjust(pval, method = "BH"),
           contrast = factor(contrast, levels = unique(contrast)),
           gs_subcat = "kinase",
           leadingEdge_symbol = map(.x = leadingEdge,
                                    .f = ~ human_uniprot_to_symbol[.x]),
           leadingEdge_rno = map(.x = leadingEdge,
                                 .f = ~ human_to_rat[.x])) %>%
    dplyr::rename(kinase = pathway) %>%
    relocate(contrast, .before = leadingEdge) %>%
    select(-gs_subcat)
})

# Save
usethis::use_data(PHOSPHO_KSEA, internal = FALSE,
                  overwrite = TRUE, version = 3, compress = "bzip2")

