## MSnSet Objects ------------------------------------------------------------

#' @title Proteomics MSnSet
#'
#' @description An object of class
#'   \code{\link[MSnbase:MSnSet-class]{MSnbase::MSnSet}} containing
#'   median-MAD-normalized log\eqn{_2}-transformed protein expression data for
#'   differential analysis and WGCNA.
#'
#' @format Object of class \code{\link[MSnbase:MSnSet-class]{MSnbase::MSnSet}}
#'   with 9964 features and 60 samples.
#'
#'   \describe{
#'     \item{\bold{assayData}}{
#'       Object of class \code{matrix} with 9964 rows and 60 columns. See
#'       \link[MotrpacRatTraining6moData]{PROT_WATSC_NORM_DATA} for more
#'       details.
#'     }
#'     \item{\bold{featureData}}{
#'       Object of class \code{data.frame} with 9964 rows and 3 variables
#'       containing additional feature information:
#'
#'       \describe{
#'          \item{\code{feature_ID}}{character; Reference Sequence (RefSeq)
#'          protein identifier.}
#'          \item{\code{gene_symbol}}{character; official gene symbol.}
#'          \item{\code{entrez_gene}}{integer; Entrez gene identifier.}
#'       }
#'     }
#'     \item{\bold{phenoData}}{
#'       Object of class \code{data.frame} with 60 rows and 7 variables
#'       containing additional sample information:
#'
#'       \describe{
#'         \item{\code{pid}}{integer; randomly generated 8-digit identifier used
#'         in linkage to phenotypic data. All samples from the same animal have
#'         the same PID.}
#'         \item{\code{bid}}{integer; unique 5 digit identifier of all
#'         samples collected for an acute test/sample collection period. All
#'         samples collected during that period will have the same BID.}
#'         \item{\code{labelid}}{integer; unique 11 digit specimen label
#'         identifier, originating at the collection site, that provides
#'         a link to specimen processing and is used for shipments to the
#'         biorepository. Same as \code{viallabel} only in instances
#'         where aliquots are not further processed at the biorepository.}
#'         \item{\code{viallabel}}{character; unique 11 digit sample vial
#'         identifier. Starts with the \code{bid}.}
#'         \item{\code{sex}}{factor; the sex of the rat with levels "Female" and
#'         "Male".}
#'         \item{\code{timepoint}}{factor; exercise training group. Either "SED"
#'         (sedentary) or the number of weeks of training ("1W", "2W", "4W",
#'         "8W").}
#'         \item{\code{exp_group}}{factor; experimental group. Unique
#'         combination of \code{sex} (first letter) and \code{timepoint}.}
#'       }
#'     }
#'   }
#'
#' @details \link[MotrpacRatTraining6moData]{FEATURE_TO_GENE} for details of
#'   feature-to-gene mapping.
#'
#' @seealso \code{\link[MSnbase]{MSnSet-class}},
#'   \link[MotrpacRatTraining6moData]{PROT_WATSC_NORM_DATA},
#'   \link[MotrpacRatTraining6moData]{FEATURE_TO_GENE}
#'
#' @examples
#' library(MSnbase)
#' PROT_MSNSET                   # display summary
#' exprs(PROT_MSNSET)[1:10, 1:5] # assayData (first 10 rows and 5 columns)
#' head(fData(PROT_MSNSET))      # featureData
#' head(pData(PROT_MSNSET))      # phenoData
#'
#' @keywords datasets
"PROT_MSNSET"


#' @title Phosphoproteomics MSnSet
#'
#' @description An object of class
#'   \code{\link[MSnbase:MSnSet-class]{MSnbase::MSnSet}} containing
#'   median-MAD-normalized log\eqn{_2}-transformed protein phosphorylation data
#'   for differential analysis.
#'
#' @format Object of class \code{\link[MSnbase:MSnSet-class]{MSnbase::MSnSet}}
#'   with 30304 features and 60 samples.
#'
#'   \describe{
#'     \item{\bold{assayData}}{
#'       Object of class \code{matrix} with 30304 rows and 60 columns. See
#'       \link[MotrpacRatTraining6moData]{PHOSPHO_WATSC_NORM_DATA} for
#'       more details.
#'     }
#'     \item{\bold{featureData}}{
#'       Object of class \code{data.frame} with 30304 rows and 7 variables
#'       containing additional feature information:
#'
#'       \describe{
#'          \item{\code{feature_ID}}{character; Reference Sequence (RefSeq)
#'          protein identifier followed by an underscore and position(s) of
#'          phosphorylation.}
#'          \item{\code{gene_symbol}}{character; official gene symbol.}
#'          \item{\code{entrez_gene}}{integer; Entrez gene identifier.}
#'          \item{\code{site}}{character; phosphorylation site(s) extracted from
#'          \code{feature_ID} column and stripped of trailing lowercase amino
#'          acid letters (s, t, y). Multiple sites are separated by semicolons.}
#'          \item{\code{human_feature_ID}}{character; human UniProt accession
#'          followed by an underscore and positions of phosphorylation.}
#'          \item{\code{human_uniprot}}{character; human UniProt accession
#'          extracted from \code{human_feature_ID} column.}
#'          \item{\code{human_site}}{character; human phosphorylation site(s)
#'          extracted from \code{human_feature_ID} column and stripped of
#'          trailing lowercase amino acid letters (s, t, y). Multiple sites are
#'          separated by semicolons.}
#'       }
#'     }
#'     \item{\bold{phenoData}}{
#'       Object of class \code{data.frame} with 60 rows and 7 variables
#'       containing additional sample information:
#'
#'       \describe{
#'         \item{\code{pid}}{integer; randomly generated 8-digit identifier used
#'         in linkage to phenotypic data. All samples from the same animal have
#'         the same PID.}
#'         \item{\code{bid}}{integer; unique 5 digit identifier of all
#'         samples collected for an acute test/sample collection period. All
#'         samples collected during that period will have the same BID.}
#'         \item{\code{labelid}}{integer; unique 11 digit specimen label
#'         identifier, originating at the collection site, that provides
#'         a link to specimen processing and is used for shipments to the
#'         biorepository. Same as \code{viallabel} only in instances
#'         where aliquots are not further processed at the biorepository.}
#'         \item{\code{viallabel}}{character; unique 11 digit sample vial
#'         identifier. Starts with the \code{bid}.}
#'         \item{\code{sex}}{factor; the sex of the rat with levels "Female" and
#'         "Male".}
#'         \item{\code{timepoint}}{factor; exercise training group. Either "SED"
#'         (sedentary) or the number of weeks of training ("1W", "2W", "4W",
#'         "8W").}
#'         \item{\code{exp_group}}{factor; experimental group. Unique
#'         combination of \code{sex} (first letter) and \code{timepoint}.}
#'       }
#'     }
#'   }
#'
#' @details \link[MotrpacRatTraining6moData]{RAT_TO_HUMAN_PHOSPHO} for details
#'   of mapping between rat and human phosphorylation sites.
#'
#'   \link[MotrpacRatTraining6moData]{FEATURE_TO_GENE} for details of
#'   feature-to-gene mapping.
#'
#' @seealso \code{\link[MSnbase]{MSnSet-class}},
#'   \link[MotrpacRatTraining6moData]{RAT_TO_HUMAN_PHOSPHO},
#'   \link[MotrpacRatTraining6moData]{FEATURE_TO_GENE}
#'
#' @examples
#' library(MSnbase)
#' PHOSPHO_MSNSET                   # summary
#' exprs(PHOSPHO_MSNSET)[1:10, 1:5] # assayData (first 10 rows and 5 columns)
#' head(fData(PHOSPHO_MSNSET))      # featureData
#' head(pData(PHOSPHO_MSNSET))      # phenoData
#'
#' @keywords datasets
"PHOSPHO_MSNSET"


#' @title Transcriptomics MSnSet
#'
#' @description An object of class
#'   \code{\link[MSnbase:MSnSet-class]{MSnbase::MSnSet}} containing
#'   transcriptomics data for differential analysis and (once transformed to
#'   log\eqn{_2} counts per million reads) WGCNA.
#'
#' @format Object of class \code{\link[MSnbase:MSnSet-class]{MSnbase::MSnSet}}
#'   with 16404 features and 48 samples.
#'
#'   \describe{
#'     \item{\bold{assayData}}{
#'       Object of class \code{matrix} with 16404 rows and 48 columns containing
#'       RSEM gene counts.
#'     }
#'     \item{\bold{featureData}}{
#'       Object of class \code{data.frame} with 16404 rows and 3 variables
#'       containing additional feature information:
#'
#'       \describe{
#'          \item{\code{feature_ID}}{character; Reference Sequence (RefSeq)
#'          protein identifier.}
#'          \item{\code{gene_symbol}}{character; official gene symbol.}
#'          \item{\code{entrez_gene}}{integer; Entrez gene identifier.}
#'       }
#'     }
#'     \item{\bold{phenoData}}{
#'       Object of class \code{data.frame} with 48 rows and 13 variables
#'       containing additional sample information:
#'
#'       \describe{
#'         \item{\code{pid}}{integer; randomly generated 8-digit identifier used
#'         in linkage to phenotypic data. All samples from the same animal have
#'         the same PID.}
#'         \item{\code{bid}}{integer; unique 5 digit identifier of all samples
#'         collected for an acute test/sample collection period. All samples
#'         collected during that period will have the same BID.}
#'         \item{\code{labelid}}{integer; unique 11 digit specimen label
#'         identifier, originating at the collection site, that provides a link
#'         to specimen processing and is used for shipments to the
#'         biorepository. Same as \code{viallabel} only in instances where
#'         aliquots are not further processed at the biorepository.}
#'         \item{\code{viallabel}}{character; unique 11 digit sample vial
#'         identifier. Starts with the \code{bid}.}
#'         \item{\code{sex}}{factor; the sex of the rat with levels "Female" and
#'         "Male".}
#'         \item{\code{timepoint}}{factor; exercise training group. Either "SED"
#'         (sedentary) or the number of weeks of training ("1W", "2W", "4W",
#'         "8W").}
#'         \item{\code{exp_group}}{factor; experimental group. Unique
#'         combination of \code{sex} (first letter) and \code{timepoint}.}
#'         \item{\code{rin}}{numeric; RNA integrity number (RIN).}
#'         \item{\code{pct_globin}}{numeric; percent of reads mapping to
#'         globin.}
#'         \item{\code{pct_umi_dup}}{numeric; percent of PCR duplicates as
#'         quantified with Unique Molecular Identifiers (UMIs).}
#'         \item{\code{median_5_3_bias}}{numeric; median 5'-3' bias.}
#'       }
#'     }
#'   }
#'
#' @details The original 32883 transcripts were filtered with
#'   \code{\link[edgeR]{filterByExpr}} using \code{group =
#'   TRNSCRPT_MSNSET[["exp_group"]]}. Then, library sizes were TMM normalized
#'   with \code{\link[edgeR]{calcNormFactors}}. Two outlier samples (90423017005
#'   and 90410017005) were removed before variables \code{pct_globin},
#'   \code{rin}, \code{pct_umi_dup}, and \code{median_5_3_bias} were
#'   mean-imputed, centered, and scaled.
#'
#'   \link[MotrpacRatTraining6moData]{OUTLIERS} for sample outlier information.
#'
#'   \link[MotrpacRatTraining6moData]{FEATURE_TO_GENE} for details of
#'   feature-to-gene mapping.
#'
#' @seealso \code{\link[MSnbase]{MSnSet-class}},
#'   \link[MotrpacRatTraining6moData]{TRNSCRPT_WATSC_RAW_COUNTS},
#'   \link[MotrpacRatTraining6moData]{FEATURE_TO_GENE}, \link[edgeR]{DGEList},
#'   \link[edgeR]{filterByExpr}, \link[edgeR]{calcNormFactors}
#'
#' @examples
#' library(MSnbase)
#'
#' TRNSCRPT_MSNSET                   # summary
#'
#' exprs(TRNSCRPT_MSNSET)[1:10, 1:5] # assayData (first 10 rows and 5 columns)
#'
#' head(fData(TRNSCRPT_MSNSET))      # featureData
#'
#' head(pData(TRNSCRPT_MSNSET))      # phenoData
#'
#' @keywords datasets
"TRNSCRPT_MSNSET"


#' @title Metabolomics MSnSet
#'
#' @description An object of class
#'   \code{\link[MSnbase:MSnSet-class]{MSnbase::MSnSet}} containing
#'   log\eqn{_2}-transformed metabolite data for differential analysis and
#'   WGCNA.
#'
#' @format Object of class \code{\link[MSnbase:MSnSet-class]{MSnbase::MSnSet}}
#'   with 1063 features and 50 samples.
#'
#'   \describe{
#'     \item{\bold{assayData}}{
#'       Object of class \code{matrix} with 1063 rows and 50 columns. See
#'       \link[MotrpacRatTraining6moData]{METAB_NORM_DATA_NESTED} for
#'       more details.
#'     }
#'     \item{\bold{featureData}}{
#'       Object of class \code{data.frame} with 1063 rows and 13 variables
#'       containing additional feature information:
#'
#'       \describe{
#'          \item{\code{feature_ID}}{character; metabolite RefMet identifier.}
#'          \item{\code{dataset}}{character; metabolomics platform in which the
#'          feature was measured.}
#'          \item{\code{name_in_figures}}{character; alternative feature
#'          identifier primarily used in visualizations.}
#'          \item{\code{refmet_super_class}}{character; RefMet chemical super
#'          class.}
#'          \item{\code{refmet_main_class}}{character; RefMet chemical main
#'          class.}
#'          \item{\code{refmet_sub_class}}{character; RefMet chemical sub
#'          class.}
#'          \item{\code{lipid_class}}{character; (lipids only) the class of the
#'          lipid. There are 43 unique classes.}
#'          \item{\code{chain_length}}{integer; (lipids only) total number of
#'          carbons that comprise the fatty acid chains.}
#'          \item{\code{double_bond}}{integer; total number of C=C double bonds
#'          in the fatty acid chains.}
#'          \item{\code{rt}}{numeric; retention time. Used to assign
#'          \code{lipid_class}.}
#'          \item{\code{mz}}{numeric; mass-to-charge ratio (m/z).}
#'          \item{\code{neutral_mass}}{numeric; neutral mass (g/mol).}
#'          \item{\code{formula}}{character; chemical formula.}
#'       }
#'     }
#'     \item{\bold{phenoData}}{
#'       Object of class \code{data.frame} with 50 rows and 6 variables
#'       containing additional sample information:
#'
#'       \describe{
#'         \item{\code{pid}}{integer; randomly generated 8-digit identifier used
#'         in linkage to phenotypic data. All samples from the same animal have
#'         the same PID.}
#'         \item{\code{bid}}{integer; unique 5 digit identifier of all
#'         samples collected for an acute test/sample collection period. All
#'         samples collected during that period will have the same BID.}
#'         \item{\code{labelid}}{integer; unique 11 digit specimen label
#'         identifier, originating at the collection site, that provides
#'         a link to specimen processing and is used for shipments to the
#'         biorepository. Same as \code{viallabel} only in instances
#'         where aliquots are not further processed at the biorepository.}
#'         \item{\code{sex}}{factor; the sex of the rat with levels "Female" and
#'         "Male".}
#'         \item{\code{timepoint}}{factor; exercise training group. Either "SED"
#'         (sedentary) or the number of weeks of training ("1W", "2W", "4W",
#'         "8W").}
#'         \item{\code{exp_group}}{factor; experimental group. Unique
#'         combination of \code{sex} (first letter) and \code{timepoint}.}
#'       }
#'     }
#'   }
#'
#' @seealso \link[MSnbase]{MSnSet-class},
#'   \link[MotrpacRatTraining6moData]{METAB_FEATURE_ID_MAP},
#'   \link[MotrpacRatTraining6moData]{METAB_NORM_DATA_NESTED},
#'   \href{https://www.metabolomicsworkbench.org/databases/refmet/index.php}{Metabolomics
#'   Workbench RefMet Database}
#'
#' @references Fahy, E., & Subramaniam, S. (2020). RefMet: A reference
#'   nomenclature for metabolomics. \emph{Nature Methods, 17}(12), 1173--1174.
#'   \url{https://doi.org/10.1038/s41592-020-01009-y}

#'
#' @examples
#' library(MSnbase)
#' METAB_MSNSET                   # summary
#' exprs(METAB_MSNSET)[1:10, 1:5] # assayData (first 10 rows and 5 columns)
#' head(fData(METAB_MSNSET))      # featureData
#' head(pData(METAB_MSNSET))      # phenoData
#'
#' @keywords datasets
"METAB_MSNSET"


## Differential Analysis -------------------------------------------------------

#' @title Differential analysis results
#'
#' @description Differential analysis of proteins, transcripts, phosphorylation
#'   sites, and metabolites measured in rat scWAT.
#'
#' @usage
#' PROT_DA     # proteomics
#' PHOSPHO_DA  # phosphoproteomics
#' TRNSCRPT_DA # transcriptomics
#' METAB_DA    # metabolomics
#'
#' @format A named list of 3 \code{data.frame} objects:
#'
#'   \describe{
#'     \item{"trained_vs_SED"}{All trained timepoints compared to their
#'     sex-matched sedentary control group. Example: F_1W - F_SED. Total of 8
#'     contrasts.}
#'     \item{"MvF_SED"}{Sedentary males compared to sedentary females. Only 1
#'     contrast.}
#'     \item{"MvF_exercise_response"}{Comparisons of the male and female
#'     training responses. Example: (M_8W - M_SED) - (F_8W - F_SED). Total of 4
#'     contrasts.}
#'   }
#'
#'   In addition to all columns present in the \code{MSnbase::fData} tables,
#'   each \code{data.frame} contains the following 8 variables:
#'
#'   \describe{
#'     \item{\code{logFC}}{numeric; difference in the weighted mean log\eqn{_2}
#'     values of the groups specified in the \code{contrast} column.}
#'     \item{\code{AveExpr}}{numeric; mean log\eqn{_2} relative abundance of
#'     each feature.}
#'     \item{\code{t}}{numeric; \acronym{LIMMA} moderated t-statistic.}
#'     \item{\code{P.Value}}{numeric; p-value.}
#'     \item{\code{adj.P.Val}}{numeric; Benjamini-Hochberg adjusted p-value.
#'     P-values are adjusted across all contrasts.}
#'     \item{\code{feature}}{character; uniquely defines each feature that was
#'     tested. Same as \code{feature_ID}.}
#'     \item{\code{contrast}}{factor; the contrast(s) being tested.}
#'     \item{\code{se}}{numeric; standard error (\acronym{SE}) of the
#'     \code{logFC} column.}
#'   }
#'
#' @details Differential analysis was performed with
#'   \code{\link[MotrpacRatTraining6moWAT:limma_full]{MotrpacRatTraining6moWAT::limma_full}}
#'   on the \code{MSnSet} objects.
#'
#' @seealso \code{\link[limma]{topTable}},
#'   \code{\link[MotrpacRatTraining6moWAT]{limma_full}},
#'   \code{\link{PROT_MSNSET}}, \code{\link{PHOSPHO_MSNSET}},
#'   \code{\link{TRNSCRPT_MSNSET}}, \code{\link{METAB_MSNSET}}
#'
#' @examples
#' ## Number of differential features (-ome-wide FDR < 0.05)
#' # Convenience function
#' f1 <- function(x) {
#'   purrr::map(x, with, table(contrast, adj.P.Val < 0.05))
#' }
#'
#' f1(PROT_DA)     # Proteins
#' f1(PHOSPHO_DA)  # Phosphosites
#' f1(TRNSCRPT_DA) # Transcripts
#' f1(METAB_DA)    # Metabolites
#'
#' @keywords datasets
#'
#' @name WATSC_DA
"PROT_DA"

#' @rdname WATSC_DA
#' @format NULL
#' @usage NULL
"PHOSPHO_DA"

#' @rdname WATSC_DA
#' @format NULL
#' @usage NULL
"TRNSCRPT_DA"

#' @rdname WATSC_DA
#' @format NULL
#' @usage NULL
"METAB_DA"


# Enrichment Analysis ----------------------------------------------------------

#' @title MSigDB Gene Ontology sets
#'
#' @description Gene sets from the Molecular Signatures Database (MSigDB). Used
#'   for Fast Gene Set Enrichment Analysis (FGSEA).
#'
#' @format A \code{data.frame} with 5106 rows and 5 columns:
#' \describe{
#'   \item{\code{gs_subcat}}{character; the database to which the gene set
#'   belongs. One of "GO:BP", "GO:MF", or "GO:CC".}
#'   \item{\code{gs_exact_source}}{character; unique gene set identifier.}
#'   \item{\code{gs_description}}{character; gene set description.}
#'   \item{\code{entrez_gene}}{numeric list; Entrez gene identifiers in each
#'   set.} \item{\code{set_size}}{numeric; gene set size. Lengths of
#'   \code{entrez_gene}.}
#' }
#'
#' @details Gene Ontology (Biological Process, Cellular Component, Molecular
#'   Function) pathways from version 7.5.1 of the Molecular Signatures Database
#'   (MSigDB). Obtained via \link[msigdbr]{msigdbr} (\code{organism =
#'   "Rattus norvegicus"}) and reformatted. Prefiltered to pathways with at
#'   least 15 and no more than 300 genes. Filtering by size was done because
#'   smaller gene sets tend to be less reliable, while larger gene sets tend to
#'   be less interpretable.
#'
#' @references Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M.,
#'   Mesirov, J. P., & Tamayo, P. (2015). The Molecular Signatures Database
#'   (MSigDB) hallmark gene set collection. \emph{Cell systems, 1}(6), 417--425.
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
#'   Consortium. \emph{Nature genetics, 25}(1), 25--29.
#'   \url{https://doi.org/10.1038/75556}
#'
#'   Gene Ontology Consortium (2021). The Gene Ontology resource: enriching a
#'   GOld mine. \emph{Nucleic acids research, 49}(D1), D325--D334.
#'   \url{https://doi.org/10.1093/nar/gkaa1113}
#'
#' @keywords datasets
"MSIGDB_PATHWAYS"


#' @title PhosphoSitePlus kinase-substrate relationships
#'
#' @description Kinase-substrate relationship data used to infer kinase activity
#'   from changes in phosphorylation sites.
#'
#' @format An object of class \code{data.frame}.
#'
#' @details The "Kinase_Substrate_Dataset.xlsx" file, downloaded from v6.6.0.4
#'   of PhosphoSitePlus on 2022-06-05
#'   (\url{https://www.phosphosite.org/staticDownloads}), was filtered to human
#'   kinases and substrates. Additionally, instances of auto-phosphorylation
#'   were removed.
#'
#' @references Hornbeck, P. V., Zhang, B., Murray, B., Kornhauser, J. M.,
#'   Latham, V., & Skrzypek, E. (2015). PhosphoSitePlus, 2014: mutations, PTMs
#'   and recalibrations. \emph{Nucleic acids research, 43}(Database issue),
#'   D512--D520. \url{https://doi.org/10.1093/nar/gku1267}
#'
#' @keywords datasets
"PSP_KINASE_SUBSTRATE"


#' @title Gene Ontology FGSEA results
#'
#' @usage
#' TRNSCRPT_FGSEA # transcriptomics
#' PROT_FGSEA     # proteomics
#'
#' @keywords datasets
#'
#' @name GO_FGSEA
"TRNSCRPT_FGSEA"

#' @rdname GO_FGSEA
#' @format NULL
#' @usage NULL
"PROT_FGSEA"


#' @title Proteomics MitoCarta FGSEA
#'
#' @keywords datasets
"PROT_MITOCARTA_FGSEA"


#' @title Phosphoproteomics KSEA
#'
#' @keywords datasets
"PHOSPHO_KSEA"


#' @title Metabolomics RefMet chemical subclass FGSEA
#'
#' @keywords datasets
"METAB_FGSEA"


## WGCNA -----------------------------------------------------------------------

#' @title WGCNA results
#'
#' @usage
#' PROT_WGCNA     # proteomics
#' TRNSCRPT_WGCNA # transcriptomics
#' METAB_WGCNA    # metabolomics
#'
#' @details See \code{\link[MotrpacRatTraining6moWAT]{run_WGCNA}} for a
#'   description of the output.
#'
#' @keywords datasets
#'
#' @name WATSC_WGCNA
"PROT_WGCNA"

#' @rdname WATSC_WGCNA
#' @format NULL
#' @usage NULL
"TRNSCRPT_WGCNA"

#' @rdname WATSC_WGCNA
#' @format NULL
#' @usage NULL
"METAB_WGCNA"


## Miscellaneous ---------------------------------------------------------------

#' @title Metabolomics triacylglyceride concentration
#'
#' @keywords datasets
"METAB_TG_CONC"


#' @title Mitochondrial DNA qPCR
#'
#' @description qPCR data.
#'
#' @format \code{data.frame} with 30 rows and 8 columns:
#'
#' \describe{
#'   \item{bid}{numeric; }
#'   \item{sex}{factor; the sex of the rat with 2 levels "Female" and "Male".}
#'   \item{timepoint}{factor; exercise training group. Either "SED"
#'   (sedentary) or the number of weeks of training ("1W", "2W", "4W", "8W").}
#'   \item{exp_group}{factor; unique combination of \code{sex} (first letter)
#'   and \code{timepoint}.}
#'   \item{mean_delta_CT}{numeric; mean of the duplicate \eqn{\Delta C_T} values
#'   from each sample.}
#'   \item{SE_delta_CT}{numeric; standard error of \code{mean_delta_CT}.}
#'   \item{delta_delta_CT}{numeric; \code{mean_delta_CT} values centered on the
#'   mean of the "F_SED" \code{exp_group}.}
#'   \item{relative_expr}{numeric; \eqn{2^{-\Delta \Delta C_T}} values.}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' head(MITO_DNA, 10)
"MITO_DNA"


#' @title Plasma clinical analytes
#'
#' @description A set of 9 common clinical analytes measured in plasma.
#'
#' @format An object of class \code{data.frame} with 146 rows and 21 columns.
#'
#' \describe{
#'   \item{\code{plate}}{integer; 96-well plate identifier (1--4).}
#'   \item{\code{runseq}}{integer; run order (1--146).}
#'   \item{\code{platepos}}{integer; position in each 96-well plate.}
#'   \item{\code{pid}}{integer; randomly generated 8-digit identifier used in
#'   linkage to phenotypic data. All samples from the same animal have the same
#'   PID.}
#'   \item{\code{bid}}{integer; unique 5 digit identifier of all samples
#'   collected for an acute test/sample collection period. All samples collected
#'   during that period will have the same BID.}
#'   \item{\code{labelid}}{integer; unique 11 digit specimen label identifier,
#'   originating at the collection site, that provides a link to specimen
#'   processing and is used for shipments to the biorepository. Same as
#'   \code{viallabel} only in instances where aliquots are not further processed
#'   at the biorepository.}
#'   \item{\code{omics_subset}}{logical; whether the sample was selected for
#'   -omics analysis.}
#'   \item{\code{sex}}{factor; the sex of the rat with levels "Female" and
#'   "Male".}
#'   \item{\code{timepoint}}{factor; exercise training group. Either "SED"
#'   (sedentary) or the number of weeks of training ("1W", "2W", "4W", "8W").}
#'   \item{\code{total_ketones}}{numeric; total ketone bodies (\eqn{\mu}mol/L).}
#'   \item{\code{nefa}}{numeric; non-esterified fatty acids (NEFA) (mmol/L).}
#'   \item{\code{lactate}}{numeric; lactate (mmol/L).}
#'   \item{\code{glucose}}{numeric; glucose (mg/dL).}
#'   \item{\code{glycerol}}{numeric; glycerol (mg/dL).}
#'   \item{\code{glucagon}}{numeric; glucagon (pmol).}
#'   \item{\code{insulin}}{numeric; insulin (pg/mL).}
#'   \item{\code{insulin_iu}}{numeric; insulin (\eqn{\mu}IU/mL) using the
#'   conversion 1 pg/mL = 0.023 \eqn{\mu}IU/mL provided in the package insert.}
#'   \item{\code{leptin}}{numeric; leptin (pg/mL).}
#'   \item{\code{corticosterone}}{numeric; corticosterone (ng/mL).}
#'   \item{\code{ins_gcg_ratio}}{numeric; molar insulin/glucagon ratio.
#'   Calculated as (insulin (pg/mL) / 5.804 (Da)) / glucagon (pmol).}
#'   \item{\code{homa_ir}}{numeric; Homeostatic Model Assessment for Insulin
#'   Resistance (HOMA-IR). Insulin (\eqn{\mu}IU/mL) x glucose (mg/dL) / 405.}
#' }
#'
#' @details
#' \bold{Fujifilm Wako (Osaka, Japan):}
#' \itemize{
#'   \item Total ketones: \emph{415-73301, 411-73401}
#'   \item NEFA: \emph{995-34791, 999-34691, 993-35191, 991-34891}
#' }
#'
#' \bold{Beckman (Brea, CA):}
#' \itemize{
#'   \item Lactate: \emph{A95550}
#'   \item Glucose: \emph{B24985}
#'   \item Glycerol: \emph{445850}
#' }
#'
#' \bold{Meso Scale Discovery (Rockville, MD):}
#' \itemize{
#'   \item Glucagon: \emph{K1535YK}
#'   \item Insulin, Leptin: \emph{K15158C (multiplex assay)}
#' }
#'
#' \bold{Alpco (Salem, NH):}
#' \itemize{
#'   \item Corticosterone: \emph{55-CORMS-E01}
#' }
#'
#' @keywords datasets
"PLASMA_ANALYTES"


#' @title Plasma analyte statistical analyses
#'
#' @description Statistical analyses of most clinical analytes described in \code{\link[MotrpacRatTraining6moWATData]{PLASMA_ANALYTES}}
#'
#' @usage
#' PLASMA_ANALYTE_STATS
#'
#' @keywords datasets
"PLASMA_ANALYTE_STATS"


#' @title Measures of adipocyte size
#'
#' @description Measures of adipocyte size.
#'
#' @usage
#' ADIPOCYTE_SIZE
#'
#' @format A \code{data.frame} with 55825 rows and 8 variables.
#'
#' @keywords datasets
"ADIPOCYTE_SIZE"


#' @title Analysis of adipocyte size differences
#'
#' @usage
#' ADIPOCYTE_SIZE_STATS
#'
#' @format A \code{data.frame} with 40 rows and 11 variables.
#'
#' @keywords datasets
"ADIPOCYTE_SIZE_STATS"


#' @title Over-representation analysis of WGCNA modules
#'
#' @description Results of performing Gene Ontology (GO) over-representation
#'   analysis (Fisher Exact / Hypergeometric test) on the modules produced by
#'   weighted gene co-expression network analysis (WGCNA).
#'
#' @keywords datasets
#'
#' @usage
#' PROT_MODULE_ORA     # proteomics
#' TRNSCRPT_MODULE_ORA # transcriptomics
#' METAB_MODULE_ORA    # metabolomics
#'
#' @name module_ORA
"PROT_MODULE_ORA"

#' @rdname module_ORA
#' @format NULL
#' @usage NULL
"TRNSCRPT_MODULE_ORA"

#' @rdname module_ORA
#' @format NULL
#' @usage NULL
"METAB_MODULE_ORA"
