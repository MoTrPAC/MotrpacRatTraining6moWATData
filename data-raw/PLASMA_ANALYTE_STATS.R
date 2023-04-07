library(MotrpacRatTraining6moWATData)
library(tidyverse)
library(emmeans)
library(gamlss)

theme_set(theme_minimal())

# -omics-only samples
PLASMA_ANALYTES <- filter(PLASMA_ANALYTES, omics_subset)

## Glycerol --------------------------------------------------------------------
ggplot(PLASMA_ANALYTES, aes(x = timepoint, y = glycerol)) +
  geom_point(na.rm = TRUE) +
  facet_grid(~ sex)
# There appears to be a sex * timepoint interaction effect

gly.m1 <- gamlss(glycerol ~ sex * timepoint + 0,
                 family = GG(mu.link = "log"),
                 data = PLASMA_ANALYTES)
plot(gly.m1)

summary(gly.m1)

gly.res <- emmeans(gly.m1, specs = dunnett ~ timepoint,
                   by = "sex", infer = TRUE, type = "response") %>%
  summary(which = 2) %>%
  as.data.frame() %>%
  mutate(analyte = "glycerol",
         formula = "glycerol ~ sex * timepoint + 0",
         family = "GG()")

## Glucagon --------------------------------------------------------------------
ggplot(PLASMA_ANALYTES, aes(x = timepoint, y = glucagon)) +
  geom_point(na.rm = TRUE) +
  facet_grid(~ sex)
# Small outlier in F2W group

gcg.m1 <- gamlss(glucagon ~ sex * timepoint + 0,
                 family = GG(mu.link = "log"),
                 method = RS(100),
                 data = PLASMA_ANALYTES)
plot(gcg.m1)

gcg.res <- emmeans(gcg.m1, specs = dunnett ~ timepoint,
                   by = "sex", infer = TRUE, type = "response") %>%
  summary(which = 2) %>%
  as.data.frame() %>%
  mutate(analyte = "glucagon",
         formula = "glucagon ~ sex * timepoint + 0",
         family = "GG()")

## Glucose ---------------------------------------------------------------------
ggplot(PLASMA_ANALYTES, aes(x = timepoint, y = glucose)) +
  geom_point(na.rm = TRUE) +
  facet_grid(~ sex)

glu.m1 <- gamlss(glucose ~ sex * timepoint + 0,
                 family = GA(),
                 data = PLASMA_ANALYTES)
plot(glu.m1)

# Hypothesis tests
glu.res <- emmeans(glu.m1, specs = dunnett ~ timepoint,
                   by = "sex", infer = TRUE, type = "response") %>%
  summary(which = 2) %>%
  as.data.frame() %>%
  mutate(analyte = "glucose",
         formula = "glucose ~ sex * timepoint + 0",
         family = "GA()")


## NEFA ------------------------------------------------------------------------
ggplot(PLASMA_ANALYTES, aes(x = timepoint, y = nefa)) +
  geom_point(na.rm = TRUE) +
  facet_grid(~ sex)

NEFA.m1 <- gamlss(nefa ~ sex * timepoint + 0,
                  family = GA(),
                  data = PLASMA_ANALYTES)
plot(NEFA.m1)

# Hypothesis testing
NEFA.res <- emmeans(NEFA.m1, specs = dunnett ~ timepoint,
                    by = "sex", infer = TRUE, type = "response") %>%
  summary(which = 2) %>%
  as.data.frame() %>%
  mutate(analyte = "NEFA",
         formula = "nefa ~ sex * timepoint + 0",
         family = "GA()")

# Confidence intervals for group means
emmeans(NEFA.m1, specs = ~ timepoint | sex, type = "response") %>%
  summary(which = 1)


## Leptin ----------------------------------------------------------------------
ggplot(PLASMA_ANALYTES, aes(x = timepoint, y = leptin)) +
  geom_point(na.rm = TRUE) +
  facet_grid(~ sex)

lep.m1 <- gamlss(leptin ~ sex * timepoint + 0,
                 family = GG(),
                 method = RS(100),
                 data = PLASMA_ANALYTES)
plot(lep.m1)
summary(lep.m1)

lep.res <- emmeans(lep.m1, specs = dunnett ~ timepoint,
                   by = "sex", infer = TRUE, type = "response") %>%
  summary(which = 2) %>%
  as.data.frame() %>%
  mutate(analyte = "leptin",
         formula = "leptin ~ sex * timepoint + 0",
         family = "GG()")

## HOMA-IR ---------------------------------------------------------------------
ggplot(PLASMA_ANALYTES, aes(x = timepoint, y = homa_ir)) +
  geom_point(na.rm = TRUE) +
  facet_grid(~ sex)
# Likely a sex effect, but no timepoint effect

homa.m0 <- gamlss(homa_ir ~ sex + 0,
                  family = GG(),
                  data = PLASMA_ANALYTES)

homa.m1 <- update(homa.m0, formula = . ~ sex + timepoint + 0)
gamlss::LR.test(null = homa.m0, alternative = homa.m1) # not better

plot(homa.m0)
summary(homa.m0)

homa.res <- emmeans(homa.m0, specs = pairwise ~ sex,
                    infer = TRUE, type = "response") %>%
  summary(which = 2) %>%
  as.data.frame() %>%
  mutate(analyte = "HOMA-IR",
         formula = "homa_ir ~ sex + 0",
         family = "GG()")

## Insulin ---------------------------------------------------------------------
ggplot(PLASMA_ANALYTES, aes(x = timepoint, y = insulin)) +
  geom_point(
    na.rm = TRUE,
    position = ggbeeswarm::position_beeswarm(cex = 3)
  ) +
  facet_grid(~ sex)
# likely a sex effect, but no timepoint effect

# Likely a sex effect, but no timepoint effect
ins.m0 <- gamlss(insulin ~ sex + 0,
                 family = GG(),
                 data = PLASMA_ANALYTES)

ins.m1 <- update(ins.m0, formula = . ~ sex + timepoint + 0)
gamlss::LR.test(null = ins.m0, alternative = ins.m1) # not better

plot(ins.m0)
summary(ins.m0)

ins.res <- emmeans(ins.m0, specs = pairwise ~ sex,
                   infer = TRUE, type = "response") %>%
  summary(which = 2) %>%
  as.data.frame() %>%
  mutate(analyte = "insulin",
         formula = "insulin ~ sex + 0",
         family = "GG()")


## Insulin/glucagon ratio ------------------------------------------------------
ggplot(PLASMA_ANALYTES, aes(x = timepoint, y = ins_gcg_ratio)) +
  geom_point(na.rm = TRUE) +
  facet_grid(~ sex)
# sex effect, but likely no timepoint effect

ig.m0 <- gamlss(insulin / 5.804 ~ glucagon + 0,
                weights = 1 / glucagon,
                family = GA("log"),
                data = PLASMA_ANALYTES)
ig.m1 <- update(ig.m0, formula = . ~ sex * glucagon + 0)
gamlss::LR.test(null = ig.m0, alternative = ig.m1) # model with sex is better

ig.m2 <- update(ig.m0, formula = . ~ sex * glucagon * timepoint + 0)
gamlss::LR.test(null = ig.m1, alternative = ig.m2) # not better

plot(ig.m1)
summary(ig.m1)

# How to interpret results? 60.281 is the mean of glucagon (default)
# Use min and max of glucagon instead
ig.res <- emmeans(ig.m1, specs = pairwise ~ sex,
                  by = "glucagon", cov.reduce = range,
                  infer = TRUE, type = "response") %>%
  summary(which = 2) %>%
  as.data.frame() %>%
  mutate(analyte = "insulin/glucagon",
         formula = "insulin / 5.804 ~ sex * glucagon + 0",
         family = "GA()",
         weights = "1 / glucagon")

# Change in the range of glucagon by sex
emmip(ig.m1, sex ~ glucagon, cov.reduce = range)

# Confidence intervals for group means
emmeans(ig.m1, specs = pairwise ~ sex * glucagon,
        type = "response") %>%
  summary(which = 1)


## Combined results ------------------------------------------------------------
PLASMA_ANALYTE_STATS <- list(
  ins.res,  # insulin
  homa.res, # HOMA-IR
  gly.res,  # glycerol
  gcg.res,  # glucagon
  glu.res,  # glucose
  lep.res,  # leptin
  NEFA.res, # NEFA
  ig.res   # insulin/glucagon
) %>%
  map(function(x) {
    x %>%
      dplyr::rename_with(.fn = ~ sub("asymp\\.L", "lower\\.", .x)) %>%
      dplyr::rename_with(.fn = ~ sub("asymp\\.U", "upper\\.", .x))
  }) %>%
  data.table::rbindlist(fill = TRUE) %>%
  mutate(glucagon = round(as.numeric(as.character(glucagon)), 2),
         stat_type = ifelse(is.infinite(df), "z", "t"),
         stat = ifelse(stat_type == "z", z.ratio, t.ratio),
         signif = cut(p.value, include.lowest = TRUE, right = FALSE,
                      breaks = c(0, 1e-3, 1e-2, 0.05, 1),
                      labels = c("***", "**", "*", ""))) %>%
  dplyr::select(-c(z.ratio, t.ratio)) %>%
  dplyr::select(analyte, sex, everything())

# Save
usethis::use_data(PLASMA_ANALYTE_STATS, internal = FALSE,
                  overwrite = TRUE, version = 3, compress = "bzip2")

