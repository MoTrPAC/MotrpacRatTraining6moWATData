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
expr_mat <- PROT_WATSC_NORM_DATA %>%
  select(feature_ID, where(is.numeric)) %>%
  tibble::column_to_rownames("feature_ID") %>%
  as.matrix()
# dim(expr_mat) # 9964  60

# Subset and reorder samples
p_data <- p_data[colnames(expr_mat), ]

# Feature data
f_data <- FEATURE_TO_GENE %>%
  filter(feature_ID %in% rownames(expr_mat)) %>%
  select(feature_ID, gene_symbol, entrez_gene) %>%
  distinct() %>%
  # For each protein, remove rows with a missing gene symbol unless
  # all gene symbols are missing
  group_by(feature_ID) %>%
  filter(!(is.na(gene_symbol) & any(!is.na(gene_symbol)))) %>%
  as.data.frame() %>%
  `rownames<-`(.[["feature_ID"]]) %>%
  .[rownames(expr_mat), ] # reorder features

# Create MSnset
PROT_MSNSET <- MSnSet(exprs = expr_mat, fData = f_data, pData = p_data)

# Save
usethis::use_data(PROT_MSNSET, internal = FALSE, overwrite = TRUE,
                  version = 3, compress = "bzip2")

