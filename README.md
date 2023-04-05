
# MotrpacRatTraining6moWATData

<!-- badges: start -->

![R package
version](https://img.shields.io/github/r-package/v/PNNL-Comp-Mass-Spec/MotrpacRatTraining6moWATData?label=R%20package)
[![R-CMD-check](https://github.com/PNNL-Comp-Mass-Spec/MotrpacRatTraining6moWATData/actions/workflows/check-and-pkgdown.yaml/badge.svg)](https://github.com/PNNL-Comp-Mass-Spec/MotrpacRatTraining6moWATData/actions/workflows/check-and-pkgdown.yaml)
![Last
commit](https://img.shields.io/github/last-commit/PNNL-Comp-Mass-Spec/MotrpacRatTraining6moWATData/master)

<!-- badges: end -->

Contains data and code necessary to reproduce analyses that appear in
“Sexual dimorphism and the multi-omic response to exercise training in
rat subcutaneous white adipose tissue”
(<https://doi.org/10.1101/2023.02.03.527012>).

- Differential analysis of transcriptomics, proteomics,
  phosphoproteomics, and metabolomics subcutaneous white adipose tissue
  -omics datasets.
- Enrichment analyses (FGSEA or KSEA) for each -ome.
- Weighted Gene Co-expression Network Analysis for all -omes excluding
  proteomics.
- Over-representation analysis (ORA) performed on WGCNA
  modules/clusters.

The
[data-raw/](https://github.com/PNNL-Comp-Mass-Spec/MotrpacRatTraining6moWATData/tree/master/data-raw)
folder contains code to reproduce all analyses. Most R script names
match the objects they create, though some scripts create multiple
objects (such as WATSC_DA.R), and so are given more general names. Code
to reproduce figures will likely be added as vignettes at a later date.

## Installation

``` r
devtools::install_github("PNNL-Comp-Mass-Spec/MotrpacRatTraining6moWATData")
```

## MotrpacRatTraining6moWAT

The <a
href="https://pnnl-comp-mass-spec.github.io/MotrpacRatTraining6moWAT/"
target="_blank">MotrpacRatTraining6moWAT</a> R package contains helper
functions for analysis and visualization of the data contained in
MotrpacRatTraining6moWATData.

## Getting Help

For questions, bug reporting, and data requests for this package, please
<a
href="https://github.com/PNNL-Comp-Mass-Spec/MotrpacRatTraining6moWATData/issues"
target="_blank">submit a new issue</a> and include as many details as
possible.

If the concern is related to functions provided in the <a
href="https://github.com/PNNL-Comp-Mass-Spec/MotrpacRatTraining6moWAT"
target="_blank">MotrpacRatTraining6moWAT</a> package, please submit an
issue <a
href="https://github.com/PNNL-Comp-Mass-Spec/MotrpacRatTraining6moWAT/issues"
target="_blank">here</a> instead.

## Acknowledgements

MoTrPAC is supported by the National Institutes of Health (NIH) Common
Fund through cooperative agreements managed by the National Institute of
Diabetes and Digestive and Kidney Diseases (NIDDK), National Institute
of Arthritis and Musculoskeletal Diseases (NIAMS), and National
Institute on Aging (NIA).

Specifically, the MoTrPAC Study is supported by NIH grants U24OD026629
(Bioinformatics Center), U24DK112349, U24DK112342, U24DK112340,
U24DK112341, U24DK112326, U24DK112331, U24DK112348 (Chemical Analysis
Sites), U01AR071133, U01AR071130, U01AR071124, U01AR071128, U01AR071150,
U01AR071160, U01AR071158 (Clinical Centers), U24AR071113 (Consortium
Coordinating Center), U01AG055133, U01AG055137 and U01AG055135
(PASS/Animal Sites).

## Data Use Agreement

Recipients and their Agents agree that in publications using **any**
data from MoTrPAC public-use data sets they will acknowledge MoTrPAC as
the source of data, including the version number of the data sets used,
e.g.:

- Data used in the preparation of this article were obtained from the
  Molecular Transducers of Physical Activity Consortium (MoTrPAC)
  database, which is available for public access at
  [motrpac-data.org](motrpac-data.org). Specific datasets used are
  \[version numbers\].
- Data used in the preparation of this article were obtained from the
  Molecular Transducers of Physical Activity Consortium (MoTrPAC)
  MotrpacRatTraining6moWATData R package \[version number\].
