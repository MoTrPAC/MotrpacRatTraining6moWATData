## Data was downloaded from "Kinase_Substrate_Dataset.xlsx" provided by PhosphoSitePlus and saved to an rda file.

# KS_data <- readxl::read_xlsx("./data/Kinase_Substrate_Dataset.xlsx")

## Reformat

# Filter to human kinases and substrates,
# remove instances of autophosphorylation
# KS_data <- KS_data %>%
#   filter(KIN_ORGANISM == "human",
#          KIN_ORGANISM == SUB_ORGANISM,
#          KIN_ACC_ID != SUB_ACC_ID)
#
# usethis::use_data(KS_data, internal = FALSE, overwrite = TRUE,
#                   version = 3, compress = "bzip2")

