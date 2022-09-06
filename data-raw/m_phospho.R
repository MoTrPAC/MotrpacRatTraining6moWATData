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


