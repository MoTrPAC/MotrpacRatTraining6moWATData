library(dplyr)

## Data was downloaded from "Kinase_Substrate_Dataset.xlsx", as provided by PhosphoSitePlus: https://www.phosphosite.org/staticDownloads

PSP_KINASE_SUBSTRATE <-
  file.path("data-raw", "Kinase_Substrate_Dataset.xlsx") %>%
  readxl::read_xlsx() %>%
  filter(KIN_ORGANISM == "human",
         SUB_ORGANISM == KIN_ORGANISM,
         KIN_ACC_ID != SUB_ACC_ID) # remove instances of autophosphorylation

# Save
usethis::use_data(PSP_KINASE_SUBSTRATE, internal = FALSE,
                  overwrite = TRUE, version = 3, compress = "bzip2")

