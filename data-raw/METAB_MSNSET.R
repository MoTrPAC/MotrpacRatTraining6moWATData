library(MotrpacRatTraining6moData)
library(dplyr)
library(tibble)
library(tidyr)
library(MSnbase)


## Expression data
expr_mat <- METAB_NORM_DATA_NESTED %>%
  lapply(function(.) {
    rownames_to_column(.[["WAT-SC"]], "feature_ID")
    }) %>%
  enframe(name = "dataset") %>%
  unnest(value)

## Remove redundant metabolites

new_refmet <- c("Car(16:1)_feature1" = "CAR 16:1 (1)",
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
                "Car(3:0, 2-CH3)" = "CAR 3:0;2Me (2)")

tmp2 <- read.delim("./data-raw/metab_tmp_data.txt")

tmp <- read.delim("./data-raw/metab_tmp_data.txt") %>%
  mutate(feature = feature_ID) %>%
  # isomers: 50:5, 56:6
  mutate(REFMET_NAME_new = ifelse(feature == "TG(36:1)>TG(8:0_10:0_18:1)_and_TG(10:0_10:0_16:1)_feature3", "TG 36:1 iso1",
                                  ifelse(feature == "TG(36:1)>TG(4:0_16:0_16:1)_and_TG(2:0_16:0_18:1)_feature4", "TG 36:1 iso2",
                                         ifelse(feature == "TG(36:1)>TG(4:0_16:0_16:1)_and_TG(2:0_16:0_18:1)_M+NH3_feature1", "TG 36:1 iso3",
                                                ifelse(feature == "TG(38:0)>TG(10:0_12:0_16:0)_and_TG(10:0_10:0_18:0)_and_TG(8:0_12:0_18:0)_M+NH3_feature3", "TG 38:0 iso1",
                                                       ifelse(feature == "TG(38:0)>TG(10:0_12:0_16:0)_and_TG(10:0_10:0_18:0)_and_TG(8:0_12:0_18:0)_M+NH3_feature2", "TG 38:0 iso2",
                                                              ifelse(feature == "TG(38:0)>TG(10:0_12:0_16:0)_and_TG(10:0_10:0_18:0)_and_TG(8:0_12:0_18:0)_M+NH3_feature1", "TG 38:0 iso3",
                                                                     ifelse(feature == "TG(38:1)>TG(10:0_10:0_18:1)_and_TG(10:0_12:0_16:1)_and_TG(8:0_12:0_18:1)_feature2", "TG 38:1 iso1",
                                                                            ifelse(feature == "TG(38:1)>TG(4:0_16:0_18:1)_feature3", "TG 38:1 iso2",
                                                                                   ifelse(feature == "TG(38:1)>TG(10:0_10:0_18:1)_and_TG(10:0_12:0_16:1)_and_TG(8:0_12:0_18:1)_M+NH3_feature1", "TG 38:1 iso3",
                                                                                          ifelse(feature == "TG(38:2)>TG(10:0_10:0_18:2)_and_TG(8:0_12:0_18:2)_feature2", "TG 38:2 iso1",
                                                                                                 ifelse(feature == "TG(38:2)>TG(4:0_16:1_18:1)_and_TG(4:0_16:0_18:2)_feature3", "TG 38:2 iso2",
                                                                                                        ifelse(feature == "TG(38:2)>TG(4:0_16:1_18:1)_and_TG(4:0_16:0_18:2)_M+NH3_feature1", "TG 38:2 iso3",
                                                                                                               CURRENT_REFMET_NAME))))))
                                                              ))))))) %>%
  select(feature_ID, dataset, name_in_figures = REFMET_NAME_new,
         lipid_class:sub_class) %>%
  dplyr::rename(refmet_super_class = super_class,
                refmet_main_class = main_class,
                refmet_sub_class = sub_class) %>%
  mutate(name_in_figures = ifelse(
    feature_ID %in% names(new_refmet),
    new_refmet[feature_ID],
    name_in_figures
  )) %>%
  `rownames<-`(NULL)

tmp[tmp == ""] <- NA

write.table(tmp, "./data-raw/nonredundant_metabolites.txt",
            quote = F, row.names = F, sep = "\t")

filtered_metabolites <- read.delim("./data-raw/nonredundant_metabolites.txt")





# ## Filtered metabolites and additional refmet data
# filtered_metabolites <- read.delim("./data-raw/nonredundant_metabolites.txt")
#
# ## Metab feature data
# f_data <- METAB_FEATURE_ID_MAP %>%
#   filter(tissue == "WAT-SC") %>%
#   select(feature_ID = feature_ID_sample_data,
#          dataset, metabolite_refmet, rt:formula)
# f_data <- inner_join(filtered_metabolites, f_data,
#                      by = c("dataset", "feature_ID")) %>%
#   `rownames<-`(.[["f"]])
#
