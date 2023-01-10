library(MotrpacRatTraining6moWATData)
library(tidyverse)
library(emmeans)

theme_set(theme_minimal())

# -omics-only samples
PLASMA_ANALYTES <- filter(PLASMA_ANALYTES, omics_subset)

## Glycerol --------------------------------------------------------------------
ggplot(PLASMA_ANALYTES, aes(x = timepoint, y = log(glycerol))) +
  geom_point(na.rm = TRUE) +
  facet_grid(~ sex)
# There appears to be a sex*timepoint interaction effect

# WLS
gly.wt <- group_by(PLASMA_ANALYTES, sex, timepoint) %>%
  mutate(wt = 1 / var(log(glycerol), na.rm = T)) %>%
  pull(wt)

gly.m1 <- lm(log(glycerol) ~ sex*timepoint,
             weights = gly.wt,
             data = PLASMA_ANALYTES)
plot(gly.m1, which = 1)
plot(gly.m1, which = 2)
plot(gly.m1, which = 3)
plot(gly.m1, which = 4)
plot(gly.m1, which = 5)

summary(gly.m1)

gly.res <- emmeans(gly.m1, specs = dunnett ~ timepoint,
                   weights = "cells", by = "sex") %>%
  summary(which = 2, infer = TRUE) %>%
  as.data.frame() %>%
  mutate(analyte = "glycerol",
         formula = "log(glycerol) ~ sex*timepoint")

## Glucagon --------------------------------------------------------------------
ggplot(PLASMA_ANALYTES, aes(x = timepoint, y = log(glucagon))) +
  geom_point(na.rm = TRUE) +
  facet_grid(~ sex)
# Small outlier in F2W group

# WLS
gcg.wt <- group_by(PLASMA_ANALYTES, sex, timepoint) %>%
  mutate(wt = 1 / var(log(glucagon), na.rm = T)) %>%
  pull(wt)

gcg.m1 <- lm(log(glucagon) ~ sex*timepoint,
             weights = gcg.wt,
             data = PLASMA_ANALYTES)
plot(gcg.m1, which = 1) # obs. 6 may be outlier
plot(gcg.m1, which = 5)

PLASMA_ANALYTES[6, c("viallabel", "sex", "timepoint", "glucagon")]
#     viallabel    sex timepoint glucagon
# 6 90578013111 Female        2W      4.3
# This is the minimum glucagon value

# Remove obs. 6 and refit

gcg.wt2 <- filter(PLASMA_ANALYTES, glucagon != 4.3) %>%
  group_by(sex, timepoint) %>%
  mutate(wt = 1 / var(log(glucagon), na.rm = T)) %>%
  pull(wt)

gcg.m2 <- lm(log(glucagon) ~ sex*timepoint,
             weights = gcg.wt2,
             data = PLASMA_ANALYTES[-6, ])
plot(gcg.m2, which = 1)
plot(gcg.m2, which = 2)
# looks good

# Compare coefficients (% change)
(coef(summary(gcg.m1))[, 1] - coef(summary(gcg.m2))[, 1]) /
  coef(summary(gcg.m2))[, 1] * 100
# 2W coefficients change quite a bit. Report both results

# Hypothesis tests
gcg.res1 <- emmeans(gcg.m1, specs = dunnett ~ timepoint,
                    weights = "cells", by = "sex") %>%
  summary(which = 2, infer = TRUE) %>%
  as.data.frame() %>%
  mutate(analyte = "glucagon",
         formula = "log(glucagon) ~ sex*timepoint")

gcg.res2 <- emmeans(gcg.m2, specs = dunnett ~ timepoint,
                    weights = "cells", by = "sex") %>%
  summary(which = 2, infer = TRUE) %>%
  as.data.frame() %>%
  mutate(analyte = "glucagon",
         obs_removed = PLASMA_ANALYTES$viallabel[6],
         formula = "log(glucagon) ~ sex*timepoint")

## Glucose ---------------------------------------------------------------------
ggplot(PLASMA_ANALYTES, aes(x = timepoint, y = log(glucose))) +
  geom_point(na.rm = TRUE) +
  facet_grid(~ sex)

# WLS
glu.wt <- group_by(PLASMA_ANALYTES, sex, timepoint) %>%
  mutate(wt = 1 / var(log(glucose), na.rm = T)) %>%
  pull(wt)

glu.m1 <- lm(log(glucose) ~ sex*timepoint,
             weights = glu.wt,
             data = PLASMA_ANALYTES)
plot(glu.m1, which = 1)
plot(glu.m1, which = 2)
plot(glu.m1, which = 3)

# Hypothesis tests
glu.res <- emmeans(glu.m1, specs = dunnett ~ timepoint,
                   by = "sex") %>%
  summary(which = 2, infer = TRUE) %>%
  as.data.frame() %>%
  mutate(analyte = "glucose",
         formula = "log(glucose) ~ sex*timepoint")

## NEFA ------------------------------------------------------------------------
ggplot(PLASMA_ANALYTES, aes(x = timepoint, y = log(nefa))) +
  geom_point(na.rm = TRUE) +
  facet_grid(~ sex)

# WLS
nefa.wt <- group_by(PLASMA_ANALYTES, sex, timepoint) %>%
  mutate(wt = 1 / var(log(nefa), na.rm = T)) %>%
  pull(wt)

NEFA.m1 <- lm(log(nefa) ~ sex*timepoint,
              weights = nefa.wt,
              data = PLASMA_ANALYTES)
plot(NEFA.m1, which = 1)
plot(NEFA.m1, which = 2)
plot(NEFA.m1, which = 3)

# Hypothesis testing
NEFA.res <- emmeans(NEFA.m1, specs = dunnett ~ timepoint,
                    by = "sex") %>%
  summary(which = 2, infer = TRUE) %>%
  as.data.frame() %>%
  mutate(analyte = "NEFA",
         formula = "log(NEFA) ~ sex*timepoint")

## Leptin ----------------------------------------------------------------------
ggplot(PLASMA_ANALYTES, aes(x = timepoint, y = log(leptin))) +
  geom_point(na.rm = TRUE) +
  facet_grid(~ sex)

# WLS
lep.wt <- group_by(PLASMA_ANALYTES, sex, timepoint) %>%
  mutate(wt = 1 / var(log(leptin), na.rm = T)) %>%
  pull(wt)

lep.m1 <- lm(log(leptin) ~ sex*timepoint,
             weights = lep.wt,
             data = PLASMA_ANALYTES)
summary(lep.m1)

plot(lep.m1, which = 1)
plot(lep.m1, which = 2)
plot(lep.m1, which = 3)
plot(lep.m1, which = 4)
plot(lep.m1, which = 5)
plot(lep.m1, which = 6)
# looks fine

# Hypothesis tests
lep.res <- emmeans(lep.m1, specs = dunnett ~ timepoint,
                   by = "sex") %>%
  summary(which = 2, infer = TRUE) %>%
  as.data.frame() %>%
  mutate(analyte = "leptin",
         formula = "log(leptin) ~ sex*timepoint")

## HOMA-IR ---------------------------------------------------------------------
ggplot(PLASMA_ANALYTES, aes(x = timepoint, y = log(homa_ir))) +
  geom_point(na.rm = TRUE) +
  facet_grid(~ sex)
# Likely a sex effect, but no interaction

# WLS
homa.wt <- group_by(PLASMA_ANALYTES, sex, timepoint) %>%
  mutate(wt = 1 / var(log(homa_ir), na.rm = T)) %>%
  pull(wt)

homa.m1 <- lm(log(homa_ir) ~ sex,
              weights = homa.wt,
              data = PLASMA_ANALYTES)
plot(homa.m1, which = 1)
plot(homa.m1, which = 2)
plot(homa.m1, which = 3)
plot(homa.m1, which = 4)
plot(homa.m1, which = 5)
plot(homa.m1, which = 6)

summary(homa.m1)

homa.res <- emmeans(homa.m1, specs = pairwise ~ sex) %>%
  summary(which = 2, infer = TRUE) %>%
  as.data.frame() %>%
  mutate(analyte = "HOMA-IR",
         formula = "log(HOMA-IR) ~ sex")



## Insulin ---------------------------------------------------------------------
ggplot(PLASMA_ANALYTES, aes(x = timepoint, y = log(insulin))) +
  geom_point(
    na.rm = TRUE,
    position = ggbeeswarm::position_beeswarm(cex = 3)
  ) +
  facet_grid(~ sex)
# likely a sex effect, but no timepoint effect

# WLS
ins.wt <- group_by(PLASMA_ANALYTES, sex, timepoint) %>%
  mutate(wt = 1 / var(log(insulin), na.rm = T)) %>%
  pull(wt)

ins.m1 <- lm(log(insulin) ~ sex,
             weights = ins.wt,
             data = PLASMA_ANALYTES)
plot(ins.m1, which = 1)
plot(ins.m1, which = 2)
plot(ins.m1, which = 3)
plot(ins.m1, which = 4)
plot(ins.m1, which = 5)
plot(ins.m1, which = 6)

summary(ins.m1)

ins.res <- emmeans(ins.m1, specs = pairwise ~ sex) %>%
  summary(which = 2, infer = TRUE) %>%
  as.data.frame() %>%
  mutate(analyte = "insulin",
         formula = "log(insulin) ~ sex")


## Insulin/glucagon ratio ------------------------------------------------------
ggplot(PLASMA_ANALYTES, aes(x = timepoint, y = log(ins_gcg_ratio))) +
  geom_point(na.rm = TRUE) +
  facet_grid(~ sex)

# WLS
ig.wt <- PLASMA_ANALYTES %>%
  group_by(sex, timepoint) %>%
  mutate(wt = 1 / var(log(ins_gcg_ratio), na.rm = T)) %>%
  pull(wt)

ig.m0 <- lm(log(ins_gcg_ratio) ~ sex + timepoint,
            weights = ig.wt,
            data = PLASMA_ANALYTES)
ig.m1 <- update(ig.m0, formula = . ~ sex*timepoint)
anova(ig.m0, ig.m1)
AIC(ig.m0, ig.m1)
# include interaction

plot(ig.m1, which = 1)
plot(ig.m1, which = 2)
plot(ig.m1, which = 3)
plot(ig.m1, which = 4)
plot(ig.m1, which = 5)
plot(ig.m1, which = 6)
# obs. 6 looks like an outlier. This one was an outlier in glucagon.
# We will remove and report both results

ig.wt2 <- PLASMA_ANALYTES %>%
  slice(-6) %>% # remove outlier
  group_by(sex, timepoint) %>%
  mutate(wt = 1 / var(log(ins_gcg_ratio), na.rm = T)) %>%
  pull(wt)

ig.m2 <- lm(log(ins_gcg_ratio) ~ sex*timepoint,
            weights = ig.wt2,
            data = PLASMA_ANALYTES[-6, ])


ig.res1 <- emmeans(ig.m1, specs = dunnett ~ timepoint,
                   by = "sex") %>%
  summary(which = 2, infer = TRUE) %>%
  as.data.frame() %>%
  mutate(analyte = "insulin/glucagon",
         formula = "log(insulin/glucagon) ~ sex*timepoint")

ig.res2 <- emmeans(ig.m2, specs = dunnett ~ timepoint,
                   by = "sex") %>%
  summary(which = 2, infer = TRUE) %>%
  as.data.frame() %>%
  mutate(analyte = "insulin/glucagon",
         obs_removed = PLASMA_ANALYTES$viallabel[6],
         formula = "log(insulin/glucagon) ~ sex*timepoint")


## Combined results ------------------------------------------------------------
PLASMA_ANALYTE_STATS <- list(
  ins.res,  # insulin
  homa.res, # HOMA-IR
  gly.res,  # glycerol
  gcg.res1, # glucagon
  gcg.res2, # glucagon (1 outlier removed)
  glu.res,  # glucose
  lep.res,  # leptin
  NEFA.res, # NEFA
  ig.res1,  # insulin/glucagon
  ig.res2   # insulin/glucagon molar ratio (1 outlier removed)
) %>%
  data.table::rbindlist(fill = TRUE) %>%
  mutate(signif = cut(p.value, include.lowest = TRUE, right = FALSE,
                      breaks = c(0, 1e-3, 1e-2, 0.05, 1),
                      labels = c("***", "**", "*", ""))) %>%
  select(analyte, obs_removed, formula, sex, contrast, everything())

# Save
usethis::use_data(PLASMA_ANALYTE_STATS, internal = FALSE,
                  overwrite = TRUE, version = 3, compress = "bzip2")

