library(MotrpacRatTraining6moData)
library(MotrpacRatTraining6moWATData)
library(MotrpacRatTrainingPhysiologyData)
library(dplyr)


# Extract necessary physiology data
PHENO_WAT <- MotrpacRatTrainingPhysiologyData::PHYSIO %>%
  filter(age == "6M") %>%
  select(sex:group, matches(".*_.*_weight"),
         matches("nmr_.*_fat.*"), matches("nmr_.*_lean"),
         matches(".*vo2max_ml_kg_min$"),
         -age, -starts_with("term")) %>%
  relocate(sex, .after = omics_analysis) %>%
  rename(timepoint = group) %>%
  rename_with(~ sub("^vo2_|^nmr_", "", .x),
              .cols = matches("vo2max|nmr")) %>%
  filter(!timepoint %in% c("1W", "2W")) %>%
  droplevels.data.frame() %>%
  mutate(pre_vo2max_ml_kg_lean_min = pre_vo2max_ml_kg_min *
           vo2_pre_weight / pre_lean,
         post_vo2max_ml_kg_lean_min = post_vo2max_ml_kg_min *
           vo2_post_weight / post_lean)

# Save
usethis::use_data(PHENO_WAT, overwrite = TRUE,
                  compress = "bzip2", version = 3)

