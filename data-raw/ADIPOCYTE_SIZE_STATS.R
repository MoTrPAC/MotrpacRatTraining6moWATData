library(MotrpacRatTraining6moWATData)
library(MASS) # glm.nb
library(emmeans)
library(tidyverse)

# Bin adipocytes in 5 micron intervals by diameter
ADIPOCYTE_SIZE <- ADIPOCYTE_SIZE %>%
  mutate(diameter_bin = cut(diameter,
                            breaks = c(14.16, seq(20, 60, 5), 62.32),
                            include.lowest = TRUE, right = FALSE,
                            ordered_result = TRUE))

# Count adipocytes by bin and experimental group
count_summary <- ADIPOCYTE_SIZE %>%
  group_by(bid, sex, timepoint, diameter_bin) %>%
  summarise(binned_adipocytes = n()) %>%
  group_by(bid) %>%
  mutate(total_adipocytes = sum(binned_adipocytes),
         adipocyte_prop = binned_adipocytes / total_adipocytes) %>%
  ungroup()

# Negative Binomial GLM with log link and offset
nb_mod <- glm.nb(binned_adipocytes ~ sex*timepoint*diameter_bin +
                   offset(log(total_adipocytes)),
                 data = count_summary)
summary(nb_mod)

plot(nb_mod, which = 1)
plot(nb_mod, which = 2)
plot(nb_mod, which = 3)
plot(nb_mod, which = 4)

ADIPOCYTE_SIZE_STATS <- emmeans(nb_mod, specs = dunnett ~ timepoint,
                                by = c("diameter_bin", "sex")) %>%
  summary(which = 2, infer = TRUE, type = "response") %>%
  as.data.frame() %>%
  mutate(diameter_bin = ordered(diameter_bin))

# Save
usethis::use_data(ADIPOCYTE_SIZE_STATS, internal = FALSE,
                  overwrite = TRUE, version = 3, compress = "bzip2")
