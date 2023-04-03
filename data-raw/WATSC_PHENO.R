library(MotrpacRatTraining6moData)
library(MotrpacRatTraining6moWATData)
library(tidyverse)


x <- filter(PHENO, tissue == "WAT-SC")

## VO2max
vo2max <- x %>%
  dplyr::select(viallabel, starts_with("vo2")) %>%
  pivot_longer(cols = -viallabel,
               names_to = c(".value", "group"),
               names_pattern = "vo2\\.max\\.test\\.(.*)_(\\d{1})$") %>%
  mutate(group = ifelse(group == 1, "pre", "post"))

## NMR
nmr <- x %>%
  dplyr::select(viallabel, starts_with("nmr")) %>%
  pivot_longer(cols = -viallabel,
               names_to = c(".value", "group"),
               names_pattern = "nmr\\.testing\\.(.*)_(\\d{1})$") %>%
  mutate(group = ifelse(group == 1, "pre", "post"))

## Training info (treadmill speed, blood lactate, weight)
## Up to 40 days (5 days/wk for 8 weeks) per viallabel
training <- x %>%
  dplyr::select(viallabel, starts_with("training")) %>%
  pivot_longer(cols = -c(viallabel, training.days_visit),
               names_to = c("day", ".value"),
               names_pattern = "^training\\.day(\\d{1,2})_?(.*)$") %>%
  mutate(day = as.numeric(day)) %>%
  filter(!is.na(day))

## Familiarization
famil <- x %>%
  dplyr::select(viallabel, starts_with("familiarization")) %>%
  rename_all(~ sub("^familiarization\\.", "", .x)) %>%
  select(-c(siteid, compliant)) # one unique value

## Terminal weights
term <- x %>%
  dplyr::select(viallabel, starts_with("terminal")) %>%
  rename_all(~ sub("^terminal\\.weight\\.", "", .x))

## Calculated variables
calc_vars <- x %>%
  dplyr::select(viallabel, starts_with("calculated")) %>%
  rename_all(~ sub("^calculated\\.variables\\.", "", .x))

## Other useful info (sex, age, d_sacrifice, etc.)
other <- x %>%
  dplyr::select(viallabel, pid, bid,
                key.agegroup, key.d_arrive, key.d_sacrifice, timepoint = group,
                time_to_freeze, starts_with("tissue")) %>%
  mutate(timepoint = ifelse(timepoint == "control", "SED",
                            toupper(timepoint)),
         timepoint = factor(timepoint,
                            levels = c("SED", paste0(2^(0:3), "W")))) %>%
  rename_all(~ sub("^key\\.", "", .x))

## Registration
regis <- x %>%
  dplyr::select(viallabel, starts_with("registration")) %>%
  rename_all(~ sub("^registration\\.", "", .x)) %>%
  select(-c(siteid, comments)) %>%
  left_join(other, by = "viallabel") %>%
  select(-contains("d_arrive"))

# List of data.frames
WATSC_PHENO <- list(
  registration = regis,
  familiarization = famil,
  NMR = nmr,
  training = training,
  VO2max = vo2max,
  terminal = term,
  calculated = calc_vars
)

# Save
usethis::use_data(WATSC_PHENO, overwrite = TRUE,
                  compress = "bzip2", version = 3)
