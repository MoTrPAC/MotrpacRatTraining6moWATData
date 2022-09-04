# library(edgeR)
#
# load("../../data/t70-white-adipose_prepped_RNA-seq_MSnSet.RData")
#
# # All RNA-Seq feature data
# feature_data <- read.delim("../../../master_feature_to_gene_small.txt")
# feature_data <- feature_data[feature_data$feature %in%
#                                featureNames(m1), ]
#
# library(dplyr)
#
# fData(m1) <- feature_data %>%
#   {rownames(.) <- NULL; .} %>%
#   group_by(feature) %>%
#   mutate(n = n(), all_LOC = all(grepl("^LOC", gene_symbol))) %>%
#   filter(!(n > 1 & grepl("^LOC", gene_symbol) & !all_LOC)) %>%
#   mutate(n = n()) %>%
#   summarise(across(c(entrez_id, gene_symbol), .fns = list),
#             n = unique(n)) %>%
#   tibble::column_to_rownames("feature") %>%
#   .[featureNames(m1), ]
#
# # Clinical traits for WGCNA
# clin <- readxl::read_xlsx("../../WGCNA/Zhenxin_work/WGCNA/WGCNA/WGCNA_clinical_traits.xlsx")
#
# sampleNames(m1) <- make.names(sampleNames(m1))
#
# pData(m1) <- pData(m1) %>%
#   mutate(bid = substr(viallabel, 1, 5)) %>%
#   # left_join(clin, by = "pid") %>%
#   dplyr::rename(exp_group = sex_training,
#                 timepoint = training_group) %>%
#   mutate(sex = toupper(substr(sex, 1, 1))) %>%
#   {rownames(.) <- make.names(pull(., viallabel)); .}
#
# library(edgeR)
#
# dge <- DGEList(counts = exprs(m1), group = m1$exp_group)
# keep <- filterByExpr(dge)
# dge <- dge[keep, , keep.lib.sizes = FALSE]
#
# m1 <- m1[featureNames(m1) %in% rownames(dge$counts), ]
# m1$covariates <- NULL
# m_trnscrpt <- m1
#
# fData(m_trnscrpt) <- fData(m_trnscrpt) %>%
#   mutate(across(where(is.list),
#                 ~ unlist(lapply(.x, function(xi) {
#                   paste(ifelse(xi == "NA", NA_character_, xi),
#                         collapse = ", ")
#                 }))
#   ))
#
# usethis::use_data(m_trnscrpt, internal = FALSE, overwrite = TRUE,
#                   version = 3, compress = "bzip2")

