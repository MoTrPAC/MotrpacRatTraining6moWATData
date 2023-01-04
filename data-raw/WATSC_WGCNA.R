library(motrpacWAT)
library(motrpacWATData)
library(edgeR)
library(tidyverse)


## Metabolomics (quick) --------------------------------------------------------
METAB_WGCNA <- run_WGCNA(object = METAB_MSNSET,
                         power = 12,
                         module_prefix = "M",
                         merge_modules = FALSE)

table(METAB_WGCNA$modules$moduleID)
# M0  M1  M2  M3  M4  M5  M6  M7
#  6 415 221 137  99  86  69  30

usethis::use_data(METAB_WGCNA, internal = FALSE, overwrite = TRUE,
                  version = 3, compress = "bzip2")


## Proteomics (slow) -----------------------------------------------------------
# Proportion of missing values?
prop.table(table(is.na(exprs(PROT_MSNSET)))) # ~5.7%

PROT_WGCNA <- run_WGCNA(object = PROT_MSNSET,
                        power = 12,
                        module_prefix = "P",
                        merge_modules = TRUE)

table(PROT_WGCNA$modules$moduleID)
#   P1   P2   P3   P4   P5   P6   P7   P8   P9  P10  P11
# 3984 1444 1412  734  696  440  435  235  227  192  165

usethis::use_data(PROT_WGCNA, internal = FALSE, overwrite = TRUE,
                  version = 3, compress = "bzip2")


## Transcriptomics (incredibly slow ~1 hr) -------------------------------------------

# Convert filtered counts to normalized log2 counts-per-million reads
dge <- DGEList(counts = exprs(TRNSCRPT_MSNSET),
               samples = pData(TRNSCRPT_MSNSET),
               group = TRNSCRPT_MSNSET$exp_group)
dge <- calcNormFactors(dge, method = "TMM")
exprs(TRNSCRPT_MSNSET) <- cpm(dge, log = TRUE)

TRNSCRPT_WGCNA <- run_WGCNA(object = TRNSCRPT_MSNSET,
                            power = 25, # use power = 20:30 to see plots
                            module_prefix = "T",
                            merge_modules = TRUE)

table(TRNSCRPT_WGCNA$modules$moduleID)
#   T0   T1   T2   T3   T4   T5   T6   T7   T8   T9  T10  T11  T12  T13
# 2599 4681 3267 2025 1677  540  461  326  225  220  178  113   98   33

#   T0   T1   T2   T3   T4   T5   T6   T7   T8   T9  T10  T11  T12  T13  T14
# 2587 4683 3251 2134 1517  541  448  325  226  210  186  114   98   51   33

usethis::use_data(TRNSCRPT_WGCNA, internal = FALSE, overwrite = TRUE,
                  version = 3, compress = "bzip2")


## Save xlsx worksheets to inst/supp-tables ----
all_WGCNA <- list(TRNSCRPT = TRNSCRPT_WGCNA,
                  PROT     = PROT_WGCNA,
                  METAB    = METAB_WGCNA) %>%
  list_transpose() # nice function

paths <- file.path("inst", "supp-tables",
                   sprintf("WGCNA_%s.xlsx",
                           c("modules", "module_eigenfeatures")))

map2(.x = all_WGCNA, .y = paths, .f = writexl::write_xlsx)

