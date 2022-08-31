# library(dplyr)
# library(MSnbase)
# library(readxl)
# library(tibble)
#
# Phenotype data
# pheno_data <- read_xlsx("../../data/PASS1B_T70_adipose_metadata.xlsx") %>%
#   as.data.frame() %>%
#   mutate(sex_training = paste0(sex, "_", Training),
#          sex_training = factor(sex_training))
# rownames(pheno_data) <- make.names(pheno_data$vialLabel)
#
# # Expression data
# exprs_data <- read.delim("../../data/motrpac_pass1b-06_t70-white-adipose_prot-pr_med-mad-normalized-logratio.txt") %>%
#   .[3:nrow(.), ]
#
# rownames(exprs_data) <- NULL
#
# exprs_data <- exprs_data %>%
#   column_to_rownames("viallabel") %>%
#   as.matrix()
#
# # Feature data
# feature_data <- read.delim("../../../master_feature_to_gene_small.txt") %>%
#   filter(feature %in% rownames(exprs_data))
#
# # Create MSnSet
# m1 <- MSnSet(exprs = exprs_data, pData = pheno_data)
#
#
# # Add feature data
# fData(m1) <- feature_data %>%
#   group_by(feature) %>%
#   summarise(entrez_id = paste(entrez_id, collapse = ", "),
#             gene_symbol = paste(gene_symbol, collapse = ", ")) %>%
#   column_to_rownames("feature") %>%
#   .[featureNames(m1), ]
#
#
# # Save MSnSet
# usethis::use_data(m_prot, internal = FALSE, overwrite = TRUE,
#                   version = 3, compress = "bzip2")

