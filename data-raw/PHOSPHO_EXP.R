library(MotrpacRatTraining6moData)
library(dplyr)
library(Biobase)

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
         exp_group = interaction(substr(sex, 1, 1), timepoint, sep = "_")) %>%
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
  tidyr::separate(human_feature_ID,
                  into = c("human_uniprot", "human_site"),
                  sep = "_", remove = FALSE) %>%
  mutate(site = sub(".*_", "", feature_ID),
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

## Create ExpressionSet object
data_dict <- read.delim(file.path("data-raw", "data_dictionary.txt"),
                        row.names = 1) %>%
  tibble::deframe()

phenoData <- AnnotatedDataFrame(data = p_data)
varMetadata(phenoData)[["labelDescription"]] <- data_dict[colnames(p_data)]

featureData <- AnnotatedDataFrame(data = f_data)
varMetadata(featureData)[["labelDescription"]] <-
  c("character; Reference Sequence (RefSeq) protein identifier followed by an underscore and position(s) of phosphorylation.",
    data_dict[colnames(f_data)][-1])

PHOSPHO_EXP <- ExpressionSet(assayData = expr_mat,
                             phenoData = phenoData,
                             featureData = featureData)

# Save
usethis::use_data(PHOSPHO_EXP, internal = FALSE, overwrite = TRUE,
                  version = 3, compress = "bzip2")

