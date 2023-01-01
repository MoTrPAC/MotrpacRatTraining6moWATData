#' @title MSigDB Gene Ontology sets
#'
#' @description Gene sets from the Molecular Signatures Database. Used
#'   for Fast Gene Set Enrichment Analysis (FGSEA).
#'
#' @format A `data.frame` with 6552 rows (pathways) and 5 columns:
#' \describe{
#'   \item{gs_subcat}{character; the database to which the gene set belongs. One
#'   of "GO:BP", "GO:MF", or "GO:CC".}
#'   \item{gs_exact_source}{character; unique gene set identifier.}
#'   \item{gs_description}{character; description of each gene set.}
#'   \item{entrez_gene}{numeric list; Entrez gene identifiers in each set.}
#'   \item{set_size}{numeric; gene set size. Length of `entrez_gene`.}
#' }
#'
#' @details Gene Ontology (Biological Process, Cellular Component, Molecular
#'   Function) pathways from version 7.5.1 of the Molecular Signatures Database
#'   (MSigDB). Obtained via \code{\link[msigdbr]{msigdbr}} (`organism = "Rattus
#'   norvegicus"`) and reformatted. Results filtered to pathways with at least
#'   15 and no more than 300 genes. Filtering by size was done because smaller
#'   gene sets tend to be less reliable, while larger gene sets tend to be less
#'   interpretable.
#'
#' @references Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M.,
#'   Mesirov, J. P., & Tamayo, P. (2015). The Molecular Signatures Database
#'   (MSigDB) hallmark gene set collection. *Cell systems, 1*(6), 417–425.
#'   \url{https://doi.org/10.1016/j.cels.2015.12.004}
#'
#'   Dolgalev, I. msigdbr: MSigDB Gene Sets for Multiple Organisms in a Tidy
#'   Data Format. R package version 7.5.1,
#'   \url{https:://igordot.github.io/msigdbr}
#'
#'   Ashburner, M., Ball, C. A., Blake, J. A., Botstein, D., Butler, H., Cherry,
#'   J. M., Davis, A. P., Dolinski, K., Dwight, S. S., Eppig, J. T., Harris, M.
#'   A., Hill, D. P., Issel-Tarver, L., Kasarskis, A., Lewis, S., Matese, J. C.,
#'   Richardson, J. E., Ringwald, M., Rubin, G. M., & Sherlock, G. (2000). Gene
#'   ontology: tool for the unification of biology. The Gene Ontology
#'   Consortium. *Nature genetics, 25*(1), 25–29.
#'   \url{https://doi.org/10.1038/75556}
#'
#'   Gene Ontology Consortium (2021). The Gene Ontology resource: enriching a
#'   GOld mine. *Nucleic acids research, 49*(D1), D325–D334.
#'   \url{https://doi.org/10.1093/nar/gkaa1113}
#'
#' @keywords datasets
#'
#' @md
"MSIGDB_PATHWAYS"


#' @title PhosphoSitePlus kinase-substrate relationships
#'
#' @description Kinase-substrate relationship data used to infer kinase activity
#'   from changes in phosphorylation sites.
#'
#' @details The "Kinase_Substrate_Dataset.xlsx" file, downloaded from
#'   v6.6.0.4 of PhosphoSitePlus on 2022-06-05
#'   (\url{https://www.phosphosite.org/staticDownloads}), was filtered to human
#'   kinases and substrates. Additionally, instances of auto-phosphorylation
#'   were removed.
#'
#' @references Hornbeck, P. V., Zhang, B., Murray, B., Kornhauser, J. M.,
#'   Latham, V., & Skrzypek, E. (2015). PhosphoSitePlus, 2014: mutations, PTMs
#'   and recalibrations. *Nucleic acids research, 43*(Database issue),
#'   D512--D520. \url{https://doi.org/10.1093/nar/gku1267}
#'
#' @keywords datasets
#'
#' @md
"PSP_KINASE_SUBSTRATE"


#' @title Mitochondrial DNA qPCR
#'
#' @description qPCR data.
#'
#' @format A `data.frame` with 30 rows and 8 columns:
#' \describe{
#'   \item{bid}{numeric; }
#'   \item{sex}{factor; the sex of the rat with 2 levels "Female" and "Male".}
#'   \item{timepoint}{factor; exercise training group. Either "SED"
#'   (sedentary) or the number of weeks of training ("1W", "2W", "4W", "8W").}
#'   \item{exp_group}{factor; unique combination of `sex` (first letter) and
#'   `timepoint`.}
#'   \item{mean_delta_CT}{numeric; mean of the duplicate \eqn{\Delta C_T} values
#'   from each sample.}
#'   \item{SE_delta_CT}{numeric; standard error of `mean_delta_CT`.}
#'   \item{delta_delta_CT}{numeric; `mean_delta_CT` values centered on the mean
#'   of the "F_SED" `exp_group`.}
#'   \item{relative_expr}{numeric; \eqn{2^{-\Delta \Delta C_T}} values.}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' head(MITO_DNA)
#'
#' @md
"MITO_DNA"


## MSnSets ---------------------------------------------------------------------

#' @title Proteomics MSnSet
#'
#' @description An \code{\link[MSnbase]{MSnSet-class}} object containing
#'   proteomics data from rat WAT provided by MoTrPAC.
#'
#' @format An \code{\link[MSnbase]{MSnSet-class}} object with 9964 features and
#'   60 samples.
#'
#'   ## `assayData`
#'
#'   A `matrix` with 9964 rows and 60 columns.
#'
#'   ## `featureData`
#'
#'   A `data.frame` with 9964 rows and 3 variables containing additional feature
#'   information:
#'   \describe{
#'     \item{feature_ID}{character; Reference Sequence (RefSeq) protein
#'     identifier.}
#'     \item{gene_symbol}{character list; the gene symbol associated with each
#'     protein.}
#'     \item{entrez_gene}{numeric list; the Entrez gene identifier associated
#'     with each protein.}
#'   }
#'
#'   ## `phenoData`
#'
#'   A `data.frame` with 60 rows and 7 variables containing sample information:
#'   \describe{
#'     \item{pid}{numeric; }
#'     \item{bid}{numeric; }
#'     \item{labelid}{numeric; }
#'     \item{viallabel}{character; }
#'     \item{sex}{factor; the sex of the rat with 2 levels "Female" and "Male".}
#'     \item{timepoint}{factor; exercise training group. Either "SED"
#'     (sedentary) or the number of weeks of training ("1W", "2W", "4W", "8W").}
#'     \item{exp_group}{factor; unique combination of `sex` (first letter) and
#'     `timepoint`.}
#'   }
#'
#' @examples
#' PROT_MSNSET # summary
#' library(MSnbase)
#' exprs(PROT_MSNSET)[1:6, 1:5] # assayData
#' head(fData(PROT_MSNSET))     # featureData
#' head(pData(PROT_MSNSET))     # phenoData
#'
#' @keywords datasets
#'
#' @md
"PROT_MSNSET"


#' @title Phosphoproteomics MSnSet
#'
#' @description An \code{\link[MSnbase]{MSnSet-class}} object containing
#'   phosphoproteomics data from rat WAT provided by MoTrPAC.
#'
#' @format An \code{\link[MSnbase]{MSnSet-class}} object with 9964 features and
#'   60 samples.
#'
#'   ## `assayData`
#'
#'   A `matrix` with 30304 rows and 60 columns.
#'
#'   ## `featureData`
#'
#'   A `data.frame` with 30304 rows and 7 variables containing additional
#'   feature information:
#'   \describe{
#'     \item{feature_ID}{character; }
#'     \item{gene_symbol}{character list; the gene symbol associated with each
#'     protein.}
#'     \item{entrez_gene}{numeric list; the Entrez gene identifier associated
#'     with each protein.}
#'     \item{site}{character; }
#'     \item{human_feature_ID}{character; }
#'     \item{human_uniprot}{character; }
#'     \item{human_site}{character; }
#'   }
#'
#'   ## `phenoData`
#'
#'   A `data.frame` with 60 rows and 7 variables containing sample information:
#'   \describe{
#'     \item{pid}{numeric; }
#'     \item{bid}{numeric; }
#'     \item{labelid}{numeric; }
#'     \item{viallabel}{character; }
#'     \item{sex}{factor; the sex of the rat with 2 levels "Female" and "Male".}
#'     \item{timepoint}{factor; exercise training group. Either "SED"
#'     (sedentary) or the number of weeks of training ("1W", "2W", "4W", "8W").}
#'     \item{exp_group}{factor; unique combination of `sex` (first letter) and
#'     `timepoint`.}
#'   }
#'
#' @examples
#' PHOSPHO_MSNSET # summary
#' library(MSnbase)
#' exprs(PHOSPHO_MSNSET)[1:5, 1:5] # assayData
#' head(fData(PHOSPHO_MSNSET))     # featureData
#' head(pData(PHOSPHO_MSNSET))     # phenoData
#'
#' @keywords datasets
#'
#' @md
"PHOSPHO_MSNSET"


#' @title Transcriptomics MSnSet
#'
#' @description An \code{\link[MSnbase]{MSnSet-class}} object containing
#'   transcriptomics data from rat WAT provided by MoTrPAC.
#'
#' @keywords datasets
#'
#' @format An \code{\link[MSnbase]{MSnSet-class}} object with 16443 features and
#'   48 samples.
#'
#'   ## `assayData`
#'
#'   A `matrix` of counts with 16443 rows and 48 columns.
#'
#'   ## `featureData`
#'
#'   A `data.frame` with 30304 rows and 7 variables containing additional
#'   feature information:
#'   \describe{
#'     \item{feature_ID}{character; }
#'     \item{gene_symbol}{character list; the gene symbol associated with each
#'     protein.}
#'     \item{entrez_gene}{numeric list; the Entrez gene identifier associated
#'     with each protein.}
#'     \item{site}{character; }
#'     \item{human_feature_ID}{character; }
#'     \item{human_uniprot}{character; }
#'     \item{human_site}{character; }
#'   }
#'
#'   Prior to collapsing the `entrez_id` and `gene_symbol` columns into lists,
#'   any rows where `gene_symbol` began with "LOC" were removed unless that
#'   would remove entire features. This reduces the impact of the one-to-many
#'   feature-to-gene mapping on FGSEA results.
#'
#'   ## `phenoData`
#'
#'   A `data.frame` with 48 rows and 9 variables containing sample information:
#'   \describe{
#'     \item{viallabel}{character; unique sample identifier.}
#'     \item{rin}{numeric; RNA integrity number (RIN).}
#'     \item{pct_globin}{numeric; percent of reads mapping to globin.}
#'     \item{pct_umi_dup}{numeric; percent of PCR duplicates as quantified with
#'     Unique Molecular Identifiers (UMIs).}
#'     \item{median_5_3_bias}{numeric; median 5'-3' bias.}
#'     \item{bid}{character; Unique 5 digit numeric identifier of all samples
#'     collected for an acute test/sample collection
#'     period. All samples collected during that period will have the same BID.
#'     Same as first 5 digits of `viallabel`.}
#'     \item{sex}{factor; the sex of the rat with 2 levels "Female" and "Male".}
#'     \item{timepoint}{factor; exercise training group. Either "SED"
#'     (sedentary) or the number of weeks of training ("1W", "2W", "4W", "8W").}
#'     \item{exp_group}{factor; unique combination of `sex` (first letter) and
#'     `timepoint`.}
#'   }
#'
#' @details Two outlier samples ("90423", "90410") were removed (48 remain). Then,
#'  `pct_globin`, `rin`, `pct_umi_dup`, and `median_5_3_bias` were mean-imputed,
#'  centered, and scaled.
#'
#' @examples
#' TRNSCRPT_MSNSET # summary
#' library(MSnbase)
#' exprs(TRNSCRPT_MSNSET)[1:5, 1:5] # assayData
#' head(fData(TRNSCRPT_MSNSET))     # featureData
#' head(pData(TRNSCRPT_MSNSET))     # phenoData
#'
#' @keywords datasets
#'
#' @md
"TRNSCRPT_MSNSET"


#' @title Metabolomics MSnSet
#'
#' @description An MSnSet object containing the subcutaneous white adipose
#'   tissue (scWAT) metabolomics data provided by the MoTrPAC Bioinformatics
#'   Center.
#'
#' @format An \code{\link[MSnbase]{MSnSet-class}} object with 1063 non-redundant features and 50 samples.
#'
#' @keywords datasets
"METAB_MSNSET"


## Differential analysis results -----------------------------------------------

#' @title Proteomics differential analysis
#'
#' @description Differential analysis results of 9964 proteins measured in rat
#'   scWAT.
#'
#' @format A named list of `data.frames`:
#'
#'     "trained_vs_SED": all trained timepoints compared to their sex-matched
#'     sedentary control group. Example: F_1W - F_SED. Total of 8 contrasts and
#'     9964*8 = 79,712 rows.
#'
#'     "MvF": sedentary males compared to sedentary females. Only 1 contrast,
#'     9964 rows.
#'
#'     "DD": comparisons of the male and female training responses. Example:
#'     (M_8W - M_SED) - (F_8W - F_SED). Total of 4 contrasts and 9964*4 = 39,856
#'     rows.
#'
#'     Each data.frame has 11 columns, most of which are from
#'     \code{\link[limma]{topTable}}:
#'
#'     \describe{
#'       \item{feature_ID}{character; Reference Sequence (RefSeq) protein
#'       identifier.}
#'       \item{gene_symbol}{character; gene symbol.}
#'       \item{entrez_gene}{numeric; Entrez gene identifier.}
#'       \item{logFC}{numeric; difference in the weighted mean \eqn{log_2} relative
#'       abundances of the groups specified in the "contrast" column.}
#'       \item{AveExpr}{numeric; mean \eqn{log_2} relative abundance of each
#'       protein.}
#'       \item{t}{numeric; LIMMA moderated t-statistic.}
#'       \item{P.Value}{numeric; p-value.}
#'       \item{adj.P.Val}{numeric; Benjamini-Hochberg adjusted p-value. P-values
#'       are adjusted across all contrasts.}
#'       \item{feature}{character; same as "feature_ID".}
#'       \item{contrast}{factor; the contrast(s) being tested.}
#'       \item{se}{numeric; standard error (SE) of the "logFC" column.}
#'     }
#'
#' @keywords datasets
"PROT_DA"


#' @title Phosphoproteomics differential analysis
#'
#' @description Differential analysis results of 30,304 phosphorylation sites
#'   measured in rat scWAT.
#'
#' @format A named list of `data.frames`:
#'
#'     "trained_vs_SED": all trained timepoints compared to their sex-matched
#'     sedentary control group. Example: F_1W - F_SED. Total of 8 contrasts and
#'     30,304*8 = 242,432 rows.
#'
#'     "MvF": sedentary males compared to sedentary females. Only 1 contrast,
#'     30,304 rows.
#'
#'     "DD": comparisons of the male and female training responses. Example:
#'     (M_8W - M_SED) - (F_8W - F_SED). Total of 4 contrasts and 30,304*4 =
#'     121,216 rows.
#'
#'     Each data.frame has 15 columns, most of which are from
#'     \code{\link[limma]{topTable}}:
#'
#'     \describe{
#'       \item{feature_ID}{character; combination of RefSeq protein identifier
#'       and phosphorylation sites.}
#'       \item{gene_symbol}{character; gene symbol.}
#'       \item{entrez_gene}{numeric; Entrez gene identifier.}
#'       \item{site}{character; phosphorylation site(s) extracted from the
#'       "feature_ID" column and reformatted slightly. Multiple instances of
#'       phosphorylation on the same protein are separated by semicolons (;).}
#'       \item{human_feature_ID}{character; combination of human UniProt protein
#'       identifier and phosphorylation sites. Similar to the "feature_ID"
#'       column.}
#'       \item{human_uniprot}{character; human UniProt protein identifiers.}
#'       \item{human_site}{character; equivalent human phosphorylation site.}
#'       \item{logFC}{numeric; difference in the weighted mean \eqn{log_2} relative
#'       abundances of the groups specified in the "contrast" column.}
#'       \item{AveExpr}{numeric; mean \eqn{log_2} relative abundance of each
#'       protein.}
#'       \item{t}{numeric; LIMMA moderated t-statistic.}
#'       \item{P.Value}{numeric; p-value.}
#'       \item{adj.P.Val}{numeric; Benjamini-Hochberg adjusted p-value. P-values
#'       are adjusted across all contrasts.}
#'       \item{feature}{character; same as "feature_ID".}
#'       \item{contrast}{factor; the contrast(s) being tested.}
#'       \item{se}{numeric; standard error (SE) of the "logFC" column.}
#'     }
#'
#' @keywords datasets
"PHOSPHO_DA"


#' @title Transcriptomics differential analysis
#'
#' @description Differential analysis results of 16,547 transcripts measured in
#'   rat scWAT.
#'
#' @format A named list of `data.frames`:
#'
#'     "trained_vs_SED": all trained timepoints compared to their sex-matched
#'     sedentary control group. Example: F_1W - F_SED. Total of 8 contrasts and
#'     16,547*8 = 132,376 rows.
#'
#'     "MvF": sedentary males compared to sedentary females. Only 1 contrast,
#'     16,547 rows.
#'
#'     "DD": comparisons of the male and female training responses. Example:
#'     (M_8W - M_SED) - (F_8W - F_SED). Total of 4 contrasts and 16,547*4 =
#'     66,188 rows.
#'
#'     Each data.frame has 11 columns, most of which are from
#'     \code{\link[limma]{topTable}}:
#'
#'     \describe{
#'       \item{feature_ID}{character; transcript gene identifier.}
#'       \item{gene_symbol}{character; gene symbol. For 195 of the transcripts,
#'       there is a one-to-many mapping between transcripts and genes. Multiple
#'       genes are separated by semicolons (;).}
#'       \item{entrez_gene}{numeric; Entrez gene identifier. For 180 of the
#'       transcripts, there is a one-to-many mapping between transcripts and
#'       genes. Multiple genes are separated by semicolons (;).}
#'       \item{logFC}{numeric; difference in the weighted mean \eqn{log_2} relative
#'       abundances of the groups specified in the "contrast" column.}
#'       \item{AveExpr}{numeric; mean \eqn{log_2} relative abundance of each
#'       protein.}
#'       \item{t}{numeric; LIMMA moderated t-statistic.}
#'       \item{P.Value}{numeric; p-value.}
#'       \item{adj.P.Val}{numeric; Benjamini-Hochberg adjusted p-value. P-values
#'       are adjusted across all contrasts.}
#'       \item{feature}{character; same as "feature_ID".}
#'       \item{contrast}{factor; the contrast(s) being tested.}
#'       \item{se}{numeric; standard error (SE) of the "logFC" column.}
#'     }
#'
#' @keywords datasets
"TRNSCRPT_DA"


## WGCNA results -----------------------------------------------------------

#' @title Proteomics WGCNA
#'
#' @keywords datasets
"PROT_WGCNA"


#' @title Transcriptomics WGCNA
#'
#' @keywords datasets
"TRNSCRPT_WGCNA"


#' @title Metabolomics WGCNA
#'
#' @keywords datasets
"METAB_WGCNA"


## FGSEA results -------------------------------------------------------------

#' @title Proteomics Gene Ontology FGSEA
#'
#' @keywords datasets
"PROT_FGSEA"

#' @title Proteomics MitoCarta FGSEA
#'
#' @keywords datasets
"PROT_MITOCARTA_FGSEA"


#' @title Transcriptomics Gene Ontology FGSEA
#'
#' @keywords datasets
"TRNSCRPT_FGSEA"


#' @title Metabolomics RefMet chemical subclass FGSEA
#'
#' @keywords datasets
"METAB_FGSEA"


#' @title Phosphoproteomics KSEA
#'
#' @keywords datasets
"PHOSPHO_KSEA"


#' @title Metabolomics differential analysis
#'
#' @keywords datasets
"METAB_DA"


#' @title Metabolomics triacylglyceride concentration
#'
#' @keywords datasets
"METAB_TG_CONC"


