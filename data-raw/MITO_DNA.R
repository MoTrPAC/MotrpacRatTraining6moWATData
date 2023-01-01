library(dplyr)
library(ggplot2)
library(MotrpacRatTraining6moData)


# Sample phenotype data
p_data <- PHENO %>%
  filter(tissue == "WAT-SC") %>%
  dplyr::select(bid, sex, timepoint = group) %>%
  distinct() %>%
  mutate(sex = factor(sex, levels = c("female", "male"),
                      labels = c("Female", "Male")),
         timepoint = ifelse(timepoint == "control", "SED",
                            toupper(timepoint)),
         timepoint = factor(timepoint,
                            levels = c("SED", paste0(2^(0:3), "W"))))

# qPCR data with CT values
MITO_DNA <- file.path("data-raw", "WAT_mtDNA.xlsx") %>%
  readxl::read_xlsx(skip = 2) %>%
  dplyr::select(2:4) %>%
  setNames(c("name", "bid", "CT")) %>%
  mutate(well = rep(1:(n() / 2), each = 2), # pair the values
         name = sub(" ", "_", tolower(name))) %>%
  group_by(well, bid) %>%
  summarise(delta_CT = CT[name == "mito_dloop"] -
              CT[name == "beta_actin"]) %>%
  group_by(bid) %>%
  # Summarize samples over 2 replicates
  summarise(mean_delta_CT = mean(delta_CT),
            SE_delta_CT = sd(delta_CT) / sqrt(2)) %>%
  left_join(p_data, by = "bid") %>%
  ungroup() %>%
  # ddCT is dCT normalized to control group. We use SED female mean as control
  mutate(
    delta_delta_CT = mean_delta_CT -
      mean(mean_delta_CT[timepoint == "SED" & sex == "Female"]),
    relative_expr = 2 ^ (-delta_delta_CT)
  ) %>%
  arrange(sex, timepoint) %>%
  mutate(exp_group = paste(sex, timepoint, sep = "_"),
         exp_group = factor(exp_group, levels = unique(exp_group))) %>%
  select(bid, sex, timepoint, exp_group, everything())

# Save
usethis::use_data(MITO_DNA, internal = FALSE,
                  overwrite = TRUE, version = 3, compress = "bzip2")

