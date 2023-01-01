library(motrpacWATData)
library(tidyverse)
library(MSnbase)


# Triacylglyceride concentrations
conc_data <- file.path("data-raw", "WAT_TAG_concentration.csv") %>%
  read.csv(check.names = FALSE) %>%
  dplyr::rename(bid = vialLabel) %>%
  dplyr::select(bid, any_of(featureNames(METAB_MSNSET))) %>%
  column_to_rownames("bid") %>%
  t()

# Create MSnSet
METAB_TG_CONC <- METAB_MSNSET[rownames(conc_data), colnames(conc_data)]

exprs(METAB_TG_CONC) <- conc_data

pData(METAB_TG_CONC) <- pData(METAB_TG_CONC) %>%
  mutate(total_TG = colSums(conc_data)) %>%
  group_by(sex, timepoint, exp_group) %>%
  transmute(median_total_TG = median(total_TG)) %>%
  ungroup() %>%
  as.data.frame() %>%
  `rownames<-`(sampleNames(METAB_TG_CONC))

METAB_TG_CONC <- METAB_TG_CONC[, order(METAB_TG_CONC$exp_group)]

# Median concentration by metabolite (used to select features for heatmap later)
fData(METAB_TG_CONC) <- fData(METAB_TG_CONC) %>%
  mutate(median_conc = apply(exprs(METAB_TG_CONC), 1, median),
         rank = rank(-median_conc))

# Extract median values and standardize for heatmap
METAB_TG_CONC <- METAB_TG_CONC %>%
  t() %>%
  combineFeatures(groupBy = METAB_TG_CONC$exp_group,
                  method = median,
                  cv = FALSE) %>%
  scale() %>%
  t() %>%
  .[, levels(METAB_MSNSET$exp_group)]

# Convert characters back to factors
pData(METAB_TG_CONC) <- pData(METAB_TG_CONC) %>%
  mutate(across(where(is.character), ~ factor(.x, levels = unique(.x))))

# Save
usethis::use_data(METAB_TG_CONC, internal = FALSE, overwrite = TRUE,
                  version = 3, compress = "bzip2")


