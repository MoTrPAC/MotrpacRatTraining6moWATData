library(tidyverse)

# Reformat, add new columns
PLASMA_ANALYTES <-
  file.path("data-raw",
            "20220830_PASS1B-06_clinical_analytes_updated.txt") %>%
  read.delim() %>%
  rename_with(tolower) %>%
  dplyr::rename(total_ketones = total.ketones) %>%
  mutate(plate = floor(platepos / 100),
         omics_subset = omics_subset == "x",
         sex = factor(str_to_title(sex),
                      levels = c("Female", "Male")),
         timepoint = factor(timepoint,
                            levels = c("SED", paste0(2^(0:3), "W"))),
         ins_gcg_ratio = (insulin / 5.804) / glucagon, # molar ratio
         insulin_iu = insulin * 0.023, # 1 pg/mL = 0.023 uIU/mL
         homa_ir = insulin_iu * glucose / 405) %>% # HOMA-IR
  relocate(plate, .before = runseq) %>%
  relocate(bid, viallabel, .after = pid) %>%
  relocate(insulin_iu, .after = insulin) %>%
  arrange(runseq)

# Number of missing values per column
map_int(PLASMA_ANALYTES, ~ sum(length(which(is.na(.x)))))

# Save
usethis::use_data(PLASMA_ANALYTES, internal = FALSE,
                  overwrite = TRUE, version = 3, compress = "bzip2")

