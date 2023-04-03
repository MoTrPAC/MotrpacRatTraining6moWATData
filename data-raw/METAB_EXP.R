library(MotrpacRatTraining6moData)
library(tidyverse)
library(Biobase)
library(ggplot2)

## Remove redundant metabolites ----

# Better names (is there an easier way to do this?)
new_refmet <- c(
  "Car(16:1)_feature1" = "CAR 16:1 (1)",
  "Car(16:1)" = "CAR 16:1 (2)",
  "C18:1 carnitine" = "CAR 18:1 (1)",
  "Car(18:1)_feature1" = "CAR 18:1 (2)",
  "C18:2 carnitine" = "CAR 18:2 (1)",
  "Car(18:2)_feature1" = "CAR 18:2 (2)",
  "C4 carnitine" = "CAR 4:0 (1)",
  "Car(4:0)" = "CAR 4:0 (2)",
  "C5 carnitine" = "CAR 5:0 (1)",
  "Car(5:0)" = "CAR 5:0 (2)",
  "C3-DC-CH3 carnitine" = "CAR 3:0;2Me (1)",
  "Car(3:0, 2-CH3)" = "CAR 3:0;2Me (2)",

  "TG(36:1)>TG(8:0_10:0_18:1)_and_TG(10:0_10:0_16:1)_feature3" = "TG 36:1 iso1",
  "TG(36:1)>TG(4:0_16:0_16:1)_and_TG(2:0_16:0_18:1)_feature4" = "TG 36:1 iso2",
  "TG(36:1)>TG(4:0_16:0_16:1)_and_TG(2:0_16:0_18:1)_M+NH3_feature1" = "TG 36:1 iso3",
  "TG(38:0)>TG(10:0_12:0_16:0)_and_TG(10:0_10:0_18:0)_and_TG(8:0_12:0_18:0)_M+NH3_feature3" = "TG 38:0 iso1",
  "TG(38:0)>TG(10:0_12:0_16:0)_and_TG(10:0_10:0_18:0)_and_TG(8:0_12:0_18:0)_M+NH3_feature2" = "TG 38:0 iso2",
  "TG(38:0)>TG(10:0_12:0_16:0)_and_TG(10:0_10:0_18:0)_and_TG(8:0_12:0_18:0)_M+NH3_feature1" = "TG 38:0 iso3",
  "TG(38:1)>TG(10:0_10:0_18:1)_and_TG(10:0_12:0_16:1)_and_TG(8:0_12:0_18:1)_feature2" = "TG 38:1 iso1",
  "TG(38:1)>TG(4:0_16:0_18:1)_feature3" = "TG 38:1 iso2",
  "TG(38:1)>TG(10:0_10:0_18:1)_and_TG(10:0_12:0_16:1)_and_TG(8:0_12:0_18:1)_M+NH3_feature1" = "TG 38:1 iso3",
  "TG(38:2)>TG(10:0_10:0_18:2)_and_TG(8:0_12:0_18:2)_feature2" = "TG 38:2 iso1",
  "TG(38:2)>TG(4:0_16:1_18:1)_and_TG(4:0_16:0_18:2)_feature3" = "TG 38:2 iso2",
  "TG(38:2)>TG(4:0_16:1_18:1)_and_TG(4:0_16:0_18:2)_M+NH3_feature1" = "TG 38:2 iso3"
)

# Pre-filtered data with additional columns
f_data <- file.path("data-raw", "nonredundant_METAB_fData.txt") %>%
  read.delim() %>%
  mutate(feature = feature_ID) %>%
  # isomers: 50:5, 56:6
  mutate(name_in_figures = ifelse(CURRENT_REFMET_NAME %in% names(new_refmet),
                                  new_refmet[CURRENT_REFMET_NAME],
                                  CURRENT_REFMET_NAME)) %>%
  select(feature_ID, dataset, name_in_figures,
         contains("class", ignore.case = FALSE),
         chain_length, double_bond) %>%
  dplyr::rename(refmet_super_class = super_class,
                refmet_main_class = main_class,
                refmet_sub_class = sub_class) %>%
  left_join(METAB_FEATURE_ID_MAP %>%
              filter(tissue == "WAT-SC") %>%
              select(dataset, feature_ID = metabolite_name,
                     rt, mz, neutral_mass, formula),
            by = c("feature_ID", "dataset")) %>%
  # Remove spaces in formulas
  mutate(formula = gsub(" ", "", formula)) %>%
  relocate(contains("refmet"), .after = name_in_figures) %>%
  `rownames<-`(.[["feature_ID"]])

f_data[f_data == ""] <- NA

# Phenodata
p_data <- PHENO %>%
  filter(tissue == "WAT-SC",
         grepl("georgia tech", biclabeldata.shiptositeid)) %>%
  select(pid, bid, labelid, sex, timepoint = group) %>%
  mutate(sex = factor(sex,
                      levels = c("female", "male"),
                      labels = c("Female", "Male")),
         timepoint = ifelse(timepoint == "control", "SED",
                            toupper(timepoint)),
         timepoint = factor(timepoint,
                            levels = c("SED", paste0(2^(0:3), "W"))),
         exp_group = paste(substr(sex, 1, 1), timepoint, sep = "_")) %>%
  arrange(sex, timepoint) %>%
  mutate(exp_group = factor(exp_group, levels = unique(exp_group))) %>%
  `rownames<-`(.[["labelid"]])

# Viallabel to labelid conversion table for expression data (next step)
vial_to_label <- select(PHENO, viallabel, labelid) %>%
  mutate(across(.fns = as.character))

# Normalized expression data
expr_mat <- METAB_NORM_DATA_NESTED %>%
  list_transpose() %>%
  pluck("WAT-SC") %>%
  map(rownames_to_column, var = "feature_ID") %>%
  enframe(name = "dataset") %>%
  unnest(value) %>%
  # subset to nonredundant combinations
  inner_join(select(f_data, feature_ID, dataset),
             by = c("feature_ID", "dataset")) %>%
  pivot_longer(cols = -c(feature_ID, dataset),
               names_to = "viallabel") %>%
  left_join(vial_to_label, by = "viallabel") %>%
  select(-c(viallabel, dataset)) %>%
  filter(!is.na(value)) %>%
  pivot_wider(id_cols = feature_ID,
              values_from = value, names_from = labelid) %>%
  column_to_rownames("feature_ID") %>%
  as.matrix() %>%
  .[rownames(f_data), rownames(p_data)]

## Current version is more flexible. This one is clearer
# expr_mat <- METAB_NORM_DATA_NESTED %>%
#   lapply(function(xi) {
#     xi <- xi[["WAT-SC"]]
#     # If we don't switch from viallabels to bid,
#     # none of the datasets will have the same column names,
#     # and we will end up with 200+ columns
#     colnames(xi) <- vial_to_bid[colnames(xi)]
#
#     rownames_to_column(xi, "feature_ID")
#   }) %>%
#   enframe(name = "dataset") %>%
#   unnest(value) %>%
#   inner_join(select(f_data, feature_ID, dataset),
#              by = c("feature_ID", "dataset")) %>%
#   column_to_rownames("feature_ID") %>%
#   select(-dataset) %>%
#   as.matrix() %>%
#   .[f_data$feature_ID, rownames(p_data)]

## Create ExpressionSet object
data_dict <- read.delim(file.path("data-raw", "data_dictionary.txt"),
                        row.names = 1) %>%
  tibble::deframe()

phenoData <- AnnotatedDataFrame(data = p_data)
varMetadata(phenoData)[["labelDescription"]] <- data_dict[colnames(p_data)]

featureData <- AnnotatedDataFrame(data = f_data)
varMetadata(featureData)[["labelDescription"]] <-
  c("character; metabolite RefMet identifier.",
    data_dict[colnames(f_data)][-1])

METAB_EXP <- ExpressionSet(assayData = expr_mat,
                           phenoData = phenoData,
                           featureData = featureData)

# Save
usethis::use_data(METAB_EXP, internal = FALSE, overwrite = TRUE,
                  version = 3, compress = "bzip2")

