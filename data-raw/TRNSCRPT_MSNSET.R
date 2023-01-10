library(MotrpacRatTraining6moData)
library(dplyr)
library(edgeR)
library(MSnbase)


# Phenodata
p_data <- PHENO %>%
  filter(tissue == "WAT-SC") %>%
  select(pid:viallabel, sex, timepoint = group) %>%
  mutate(sex = factor(stringr::str_to_title(sex),
                      levels = c("Female", "Male")),
         timepoint = ifelse(timepoint == "control", "SED",
                            toupper(timepoint)),
         timepoint = factor(timepoint,
                            levels = c("SED", paste0(2^(0:3), "W"))),
         exp_group = interaction(substr(sex, 1, 1), timepoint, sep = "_")) %>%
  arrange(sex, timepoint) %>%
  mutate(exp_group = factor(exp_group, levels = unique(exp_group))) %>%
  # Add columns used in differential analysis
  inner_join(select(TRNSCRPT_META, viallabel, rin = RIN,
                    pct_globin, pct_umi_dup, median_5_3_bias),
             by = "viallabel") %>%
  `rownames<-`(.[["viallabel"]])

# Normalized transcriptomics data
count_mat <- TRNSCRPT_WATSC_RAW_COUNTS %>%
  select(feature_ID, where(is.numeric)) %>%
  tibble::column_to_rownames("feature_ID") %>%
  as.matrix() %>%
  .[, rownames(p_data)]

# Feature data
f_data <- FEATURE_TO_GENE %>%
  filter(feature_ID %in% rownames(count_mat)) %>%
  select(feature_ID, gene_symbol, entrez_gene) %>%
  distinct() %>%
  # Some transcripts have more than one gene ID
  group_by(feature_ID) %>%
  # For each transcript, remove genes that start with
  # "LOC" or "NEWGENE" unless there are no other genes
  filter(!(grepl("^LOC|^NEWGENE", gene_symbol) &
             !all(grepl("^LOC|^NEWGENE", gene_symbol)))) %>%
  # If all genes start with "LOC" or "NEWGENE",
  # only keep those with Entrez IDs unless none of them have Entrez IDs
  filter(!(all(grepl("^LOC|^NEWGENE", gene_symbol)) &
             is.na(entrez_gene) & !all(is.na(entrez_gene)))) %>%
  # Collapse duplicates
  summarise(across(c(gene_symbol, entrez_gene),
                   ~ ifelse(all(is.na(.x)), NA_character_,
                            paste(.x, collapse = ";")))) %>%
  as.data.frame() %>%
  `rownames<-`(.[["feature_ID"]]) %>%
  .[rownames(count_mat), ] # reorder features


## Following 10.12688/f1000research.9005.3
# Filter lowly expressed genes
dge_raw <- DGEList(counts = count_mat,
                   samples = p_data,
                   genes = f_data,
                   group = p_data$exp_group,
                   remove.zeros = TRUE)
keep <- filterByExpr(dge_raw) # do not normalize first
dge_raw <- dge_raw[keep, , keep.lib.sizes = FALSE]
dim(dge_raw) # 16404  50

# Calculate normalization factors
dge_raw <- calcNormFactors(dge_raw, method = "TMM")

# How many transcripts have more than one gene? About 1.2%
table(grepl(";", dge_raw$genes$gene_symbol))
# FALSE  TRUE
# 16212   192

# Check for extreme outliers
CPM <- cpm(dge_raw, log = TRUE)

plotMDS(CPM, top = 1000,
        label = dge_raw$samples$bid, dim.plot = c(1, 3))
# 90423 and 90410 are extreme outliers. We will discard these.

# Remove 2 outlier samples and recalculate normalization factors
dge_raw <- dge_raw[, !colnames(dge_raw) %in% OUTLIERS$viallabel]

# Remove unnecessary columns; mean impute, center, and scale others
dge_raw[["samples"]] <- dge_raw[["samples"]] %>%
  select(-c(group, lib.size, norm.factors)) %>%
  # Use code from MotrpacRatTraining6mo::fix_covariates to
  # mean impute, center, and scale
  mutate(
    across(c(rin, pct_globin, pct_umi_dup, median_5_3_bias), function(cov) {
      cov[is.na(cov)] <- mean(cov)
      cov <- scale(cov)[, 1]
      return(cov)
    }))

# Create MSnset from DGEList
TRNSCRPT_MSNSET <- with(dge_raw,
                        MSnSet(exprs = counts, fData = genes, pData = samples))

# Save
usethis::use_data(TRNSCRPT_MSNSET, internal = FALSE, overwrite = TRUE,
                  version = 3, compress = "bzip2")

