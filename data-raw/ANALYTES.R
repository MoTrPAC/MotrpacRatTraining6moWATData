library(dplyr)

# Reformat, add new columns
ANALYTES <-
  file.path("data-raw",
            "20220830_PASS1B-06_clinical_analytes_updated.txt") %>%
  read.delim() %>%
  rename_with(tolower) %>%
  dplyr::rename(total_ketones = total.ketones) %>%
  mutate(plate = floor(platepos / 100),
         omics_subset = omics_subset == "x",
         sex = factor(sub("(.)", "\\U\\1", sex, perl = TRUE),
                      levels = c("Female", "Male")),
         timepoint = factor(timepoint,
                            levels = c("SED", paste0(2 ^ (0:3), "W"))),
         insulin_iu = insulin * 0.023, # 1 pg/mL = 0.023 uIU/mL
         insulin_pm = insulin / 5.804) %>% # divide by molecular weight
  relocate(plate, .before = runseq) %>%
  relocate(bid, viallabel, .after = pid) %>%
  arrange(runseq)

# Alphabetically order analyte columns
ANALYTES <- ANALYTES %>%
  {cbind(.[, 1:9], .[, sort(colnames(.[, 10:ncol(.)]))])}

# Save
usethis::use_data(ANALYTES, internal = FALSE,
                  overwrite = TRUE, version = 3, compress = "bzip2")

