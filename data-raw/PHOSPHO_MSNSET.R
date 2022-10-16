library(MotrpacRatTraining6moData)
library(dplyr)
library(MSnbase)


# Phenodata
p_data <- PHENO %>%
  filter(tissue == "WAT-SC") %>%
  select(pid:viallabel, sex, timepoint = group) %>%
  mutate(sex = factor(stringr::str_to_title(sex),
                      levels = c("Female", "Male")),
         timepoint = ifelse(timepoint == "control", "SED",
                            toupper(timepoint)),
         timepoint = factor(timepoint,
                            levels = c("SED", paste0(2^(0:3), "W"))),
         exp_group = interaction(sex, timepoint, sep = "_")) %>%
  arrange(sex, timepoint) %>%
  mutate(exp_group = factor(exp_group, levels = unique(exp_group)))

# Normalized proteomics data
expr_mat <- PHOSPHO_WATSC_NORM_DATA %>%
  select(feature_ID, where(is.numeric)) %>%
  tibble::column_to_rownames("feature_ID") %>%
  as.matrix() %>%
  round(digits = 4)
# dim(expr_mat) # 30304  60

p_data <- p_data[colnames(expr_mat), ]

# Rat-to-human phosphosite mapping
rat2human <- RAT_TO_HUMAN_PHOSPHO %>%
  `colnames<-`(c("feature_ID", "human_feature_ID")) %>%
  filter(!is.na(feature_ID)) %>%
  mutate(site = sub(".*_", "", feature_ID),
         human_uniprot = sub("_.*", "", human_feature_ID),
         human_site = sub(".*_", "", human_feature_ID),
         across(contains("site"),
                ~ sub(";$", "", gsub("[sty]", ";", .x)))) %>%
  select(feature_ID, site, everything())

# Feature data
f_data <- FEATURE_TO_GENE %>%
  filter(feature_ID %in% rownames(expr_mat)) %>%
  select(feature_ID, gene_symbol, entrez_gene) %>%
  distinct() %>%
  # For each phosphosite, remove rows with missing gene symbols
  # unless all genes are missing
  group_by(feature_ID) %>%
  filter(!(is.na(gene_symbol) & any(!is.na(gene_symbol)))) %>%
  as.data.frame() %>%
  left_join(rat2human, by = "feature_ID") %>%
  `rownames<-`(.[["feature_ID"]]) %>%
  .[rownames(expr_mat), ] # reorder features

# Create MSnset
PHOSPHO_MSNSET <- MSnSet(exprs = expr_mat, fData = f_data, pData = p_data)

# Save
usethis::use_data(PHOSPHO_MSNSET, internal = FALSE, overwrite = TRUE,
                  version = 3, compress = "bzip2")

