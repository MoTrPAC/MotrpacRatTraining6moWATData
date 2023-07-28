library(readxl)
library(dplyr)
library(tidyr)
# devtools::install_github("MoTrPAC/MotrpacRatTraining6moData")
library(MotrpacRatTraining6moData)

path <- file.path("data-raw",
                  "WAT_histology_adipocyte_area.xlsx")

# Phenotype information (need bid to join with histology data)
pheno <- MotrpacRatTraining6moData::PHENO %>%
  dplyr::select(pid, bid, sex, timepoint = group) %>%
  distinct() %>%
  mutate(across(.cols = c(pid, bid), as.character),
         timepoint = ifelse(timepoint == "control", "SED", toupper(timepoint)),
         timepoint = factor(timepoint,
                            levels = c("SED", paste0(2 ^ (0:3), "W"))),
         sex = factor(sex, levels = c("female", "male"),
                      labels = c("Female", "Male")))

ADIPOCYTE_SIZE <- read_xlsx(path, sheet = 2) %>%
  pivot_longer(cols = everything(), names_to = "bid", values_to = "area") %>%
  filter(!is.na(area)) %>%
  group_by(bid) %>%
  mutate(n_adipocytes = n()) %>%
  ungroup() %>%
  mutate(diameter = 2 * sqrt(area / pi), # assume circular cross-sections
         volume = 4 / 3 * pi * (diameter / 2) ^ 3, # assume spherical adipocytes
         # Bin adipocytes in 5 micron intervals by diameter
         diameter_bin = cut(diameter, dig.lab = 4,
                            breaks = c(min(diameter),
                                       seq(20, 60, 5),
                                       max(diameter)),
                            include.lowest = TRUE, right = FALSE,
                            ordered_result = FALSE)) %>%
  relocate(diameter, diameter_bin, area, volume, .before = n_adipocytes) %>%
  left_join(pheno, by = "bid") %>%
  dplyr::select(-bid) %>%
  relocate(pid, sex, timepoint, .before = everything()) %>%
  arrange(sex, timepoint, pid, diameter)

usethis::use_data(ADIPOCYTE_SIZE, internal = FALSE,
                  overwrite = TRUE, version = 3, compress = "bzip2")
