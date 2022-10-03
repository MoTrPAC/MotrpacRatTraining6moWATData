# ## Phosphoproteomics MSnSet
#
# library(dplyr)
#
# m_phospho <- m1
# featureNames(m_phospho) <- gsub("[sty]", ";", featureNames(m_phospho))
# featureNames(m_phospho) <- sub(";$", "", featureNames(m_phospho))
#
# fData(m_phospho) <- data.frame(feature = featureNames(m_phospho)) %>%
#   mutate(refseq = sub("_[STY].*", "", feature),
#          site = sub(".*_", "", feature)) %>%
#   {rownames(.) <- pull(., feature); select(., -feature)}
#
# feature_data <- read.delim("../../../master_feature_to_gene_small.txt") %>% filter(feature %in% fData(m_phospho)[["refseq"]])
#
# feature_data <- filter(feature_data, entrez_id != "100910990")
#
# tmp <- merge(fData(m_phospho), feature_data, by = c("refseq"), sort = FALSE, all.x = T, all.y = FALSE)
# tmp2 <- merge(tmp, human_sites, by = c("refseq", "site"),
#               sort = FALSE, all.x = T, all.y = FALSE)
# tmp2 <- tmp2[, c(1:4, 6, 5)]
# rownames(tmp2) <- paste0(tmp2$refseq, "_", tmp2$site)
#
# fData(m_phospho) <- tmp2[featureNames(m_phospho), ]
#
# save(m_phospho, file = "../m_phospho.RData", compress = TRUE)
#
#
# pData(m_phospho) <- pData(m_phospho) %>%
#   dplyr::rename(timepoint = Training, exp_group = sex_training) %>%
#   mutate(across(c(pid, vialLabel, bid, plex), as.character),
#          sex = factor(sex, level = c("F", "M")),
#          timepoint = factor(timepoint,
#                             levels = c("SED", paste0(2^(0:3), "W"))))
#
# # Save MSnSet
# usethis::use_data(m_phospho, internal = FALSE, overwrite = TRUE,
#                   version = 3, compress = "bzip2")


