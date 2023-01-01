library(motrpacWATData)
library(motrpacWAT) # limma_full
library(tidyverse)
library(MSnbase)


## Contrasts to test -----------------------------------------------------------
# Sex-specific training differences
contr_train <- sprintf("%s_%s - %s_SED",
                       rep(c("F", "M"), each = 4),
                       rep(paste0(2^(0:3), "W"), times = 2),
                       rep(c("F", "M"), each = 4))

# Training-induced sexual dimorphism (sex by timepoint interaction)
contr_diff <- sprintf("(%s) - (%s)",
                      contr_train[5:8],
                      contr_train[1:4])

# List of contrast groups
contr_list <- list("trained_vs_SED" = contr_train,
                   "MvF_SED" = "M_SED - F_SED",
                   "MvF_exercise_response" = contr_diff)


## Proteomics ------------------------------------------------------------------
PROT_DA <- map(contr_list, function(contrasts) {
  limma_full(eset = PROT_MSNSET,
             model.str = "~ 0 + exp_group",
             coef.str = "exp_group",
             contrasts = contrasts,
             var.group = "viallabel") %>%
    arrange(contrast, feature) %>%
    select(-B)
}, .progress = TRUE)

usethis::use_data(PROT_DA, internal = FALSE, overwrite = TRUE,
                  version = 3, compress = "bzip2")


## Phosphoproteomics -----------------------------------------------------------
PHOSPHO_DA <- map(contr_list, function(contrasts) {
  limma_full(eset = PHOSPHO_MSNSET,
             model.str = "~ 0 + exp_group",
             coef.str = "exp_group",
             contrasts = contrasts,
             var.group = "vialLabel") %>%
    arrange(contrast, feature) %>%
    select(-B)
}, .progress = TRUE)

usethis::use_data(PHOSPHO_DA, internal = FALSE, overwrite = TRUE,
                  version = 3, compress = "bzip2")


## Transcriptomics -------------------------------------------------------------
covariates <- "rin + pct_globin + pct_umi_dup + median_5_3_bias"

TRNSCRPT_DA <- map(contr_list, function(contrasts) {
  limma_full(eset = TRNSCRPT_MSNSET,
             model.str = sprintf("~ 0 + exp_group + %s", covariates),
             coef.str = "exp_group",
             contrasts = contrasts,
             var.group = "viallabel") %>%
    arrange(contrast, feature) %>%
    select(-B)
    # entrez_gene is of type character because of the one-to-many
    # transcript to gene mapping. Keep this in mind.
})

usethis::use_data(TRNSCRPT_DA, internal = FALSE, overwrite = TRUE,
                  version = 3, compress = "bzip2")


## Metabolomics ----------------------------------------------------------------

## Platforms for DEA
# Does not make sense to run eBayes on everything, so
# we subset to a specific platform and then stack results.
assays <- unique(fData(METAB_MSNSET)[["dataset"]])

## Differential analysis results list
METAB_DA <- map(contr_list, function(contrasts) {
  map(assays, function(assay) {
    message(assay)
    # subset to features in group to model separate mean-variance trends
    METAB_MSNSET[fData(METAB_MSNSET)[["dataset"]] == assay, ] %>%
      limma_full(model.str = "~ 0 + exp_group",
                 coef.str = "exp_group",
                 contrasts = contrasts,
                 var.group = "vialLabel") %>%
      arrange(contrast, feature) %>%
      select(-B)
  }) %>%
    data.table::rbindlist() %>%
    mutate(contrast = factor(contrast, levels = unique(contrast)),
           adj.P.Val = p.adjust(P.Value, method = "BH")) %>%
    arrange(contrast, feature)
})

usethis::use_data(METAB_DA, internal = FALSE, overwrite = TRUE,
                  version = 3, compress = "bzip2")


## Save xlsx worksheets to inst/supp-tables ----
all_DA <- list(TRNSCRPT = TRNSCRPT_DA,
               PROT     = PROT_DA,
               PHOSPHO  = PHOSPHO_DA,
               METAB    = METAB_DA) %>%
  list_transpose() # nice function

paths <- file.path("inst", "supp-tables",
                   sprintf("DA_%s.xlsx", names(all_DA)))

map2(.x = all_DA, .y = paths, .f = writexl::write_xlsx, .progress = TRUE)

