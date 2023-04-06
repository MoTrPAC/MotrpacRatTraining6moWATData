library(dplyr)
library(tidyr)
library(readxl)

# Human MitoCarta3.0 database
MITOCARTA_HS <- file.path("data-raw", "Human.MitoCarta3.0.xls") %>%
  readxl::read_xls(sheet = "C MitoPathways") %>%
  .[, 2:4] %>%
  setNames(c("pathway", "hierarchy", "human_genes")) %>%
  filter(!is.na(pathway)) %>%
  mutate(human_genes = strsplit(human_genes, split = ", "))

# Save
usethis::use_data(MITOCARTA_HS, internal = FALSE,
                  overwrite = TRUE, version = 3, compress = "bzip2")
