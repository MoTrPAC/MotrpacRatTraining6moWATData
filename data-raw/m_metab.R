# library(dplyr)
# library(MSnbase)
#
#
# # Sample data
# p_data <- read.delim("./data/Pass1B_T70_WhiteAdipose_Class_Labels.txt") %>%
#   transmute(viallabel = as.character(vialLabel),
#             sex = substr(class, nchar(class), nchar(class)),
#             timepoint = ifelse(grepl("^C", class),
#                                "SED", timepoint)) %>%
#   mutate(sex = factor(sex, levels = c("F", "M"),
#                       labels = c("Female", "Male")),
#          timepoint = factor(
#            timepoint, levels = c("SED", paste0(2^(0:3), "W"))
#          )
#   ) %>%
#   arrange(sex, timepoint) %>%
#   mutate(exp_group = paste0(substr(sex, 1, 1),
#                             "_", timepoint),
#          exp_group = factor(exp_group,
#                             levels = unique(exp_group))) %>%
#   {rownames(.) <- .[["viallabel"]]; .}
#
# # Abundance data
# load("./data/Met-01-dataInput.RData")
# x <- t(datMet)[, rownames(p_data)]
#
# # Add platform info
# platform_df <- read.delim("./data/WAT_Met_NonRedun_data.txt") %>%
#   select(feature_ID, dataset)
#
# # Feature data
# f_data <- read.csv("./data/WAT_Dictionary_new.csv") %>%
#   left_join(platform_df, by = "feature_ID") %>%
#   filter(feature_ID %in% rownames(x)) %>%
#   {rownames(.) <- .[["feature_ID"]]; .[rownames(x), ]}
#
# # Create MSnSet
# m_metab <- MSnSet(exprs = x, fData = f_data, pData = p_data)

# m_metab$bid <- m_metab$viallabel
# m_metab$sex <- substr(m_metab$sex, 1, 1)

# Save MSnSet
# usethis::use_data(m_metab, internal = FALSE, overwrite = TRUE,
#                   version = 3, compress = "bzip2")

