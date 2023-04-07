## ExpressionSet Objects =======================================================

#' @title Proteomics ExpressionSet
#'
#' @description An object of class
#'   \code{\link[Biobase:ExpressionSet-class]{Biobase::ExpressionSet}}
#'   containing median-MAD-normalized log\eqn{_2}-transformed protein expression
#'   data for differential analysis and WGCNA.
#'
#' @format Object of class
#'   \code{\link[Biobase:ExpressionSet-class]{Biobase::ExpressionSet}} with 9964
#'   features and 60 samples.
#'
#' @details \link[MotrpacRatTraining6moData]{FEATURE_TO_GENE} for details of
#'   feature-to-gene mapping.
#'
#' @seealso \code{\link[Biobase]{ExpressionSet-class}},
#'   \link[MotrpacRatTraining6moData]{PROT_WATSC_NORM_DATA},
#'   \link[MotrpacRatTraining6moData]{FEATURE_TO_GENE}
#'
#' @examples
#' library(Biobase)
#'
#' PROT_EXP                   # display summary
#' exprs(PROT_EXP)[1:10, 1:5] # assayData (first 10 rows and 5 columns)
#' head(fData(PROT_EXP))      # featureData
#' head(pData(PROT_EXP))      # phenoData
#'
#' @keywords datasets
"PROT_EXP"


#' @title Phosphoproteomics ExpressionSet
#'
#' @description An object of class
#'   \code{\link[Biobase:ExpressionSet-class]{Biobase::ExpressionSet}}
#'   containing median-MAD-normalized log\eqn{_2}-transformed protein
#'   phosphorylation data.
#'
#' @format Object of class
#'   \code{\link[Biobase:ExpressionSet-class]{Biobase::ExpressionSet}} with
#'   30304 features and 60 samples.
#'
#' @details \link[MotrpacRatTraining6moData]{RAT_TO_HUMAN_PHOSPHO} for details
#'   of mapping between rat and human phosphorylation sites.
#'
#'   \link[MotrpacRatTraining6moData]{FEATURE_TO_GENE} for details of
#'   feature-to-gene mapping.
#'
#' @seealso \code{\link[Biobase]{ExpressionSet-class}},
#'   \link[MotrpacRatTraining6moData]{RAT_TO_HUMAN_PHOSPHO},
#'   \link[MotrpacRatTraining6moData]{FEATURE_TO_GENE}
#'
#' @examples
#' library(Biobase)
#'
#' PHOSPHO_EXP                   # summary
#' exprs(PHOSPHO_EXP)[1:10, 1:5] # assayData (first 10 rows and 5 columns)
#' head(fData(PHOSPHO_EXP))      # featureData
#' head(pData(PHOSPHO_EXP))      # phenoData
#'
#' @keywords datasets
"PHOSPHO_EXP"


#' @title Transcriptomics ExpressionSet
#'
#' @description An object of class
#'   \code{\link[Biobase:ExpressionSet-class]{Biobase::ExpressionSet}}
#'   containing transcriptomics data for differential analysis and (once
#'   transformed to log\eqn{_2} counts per million reads) WGCNA.
#'
#' @format Object of class
#'   \code{\link[Biobase:ExpressionSet-class]{Biobase::ExpressionSet}} with
#'   16404 features and 48 samples.
#'
#' @details The original 32883 transcripts were filtered with
#'   \code{\link[edgeR]{filterByExpr}} using \code{group =
#'   TRNSCRPT_EXP[["exp_group"]]}. Then, library sizes were TMM normalized
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
#' @seealso \code{\link[Biobase]{ExpressionSet-class}},
#'   \link[MotrpacRatTraining6moData]{TRNSCRPT_WATSC_RAW_COUNTS},
#'   \link[MotrpacRatTraining6moData]{FEATURE_TO_GENE}, \link[edgeR]{DGEList},
#'   \link[edgeR]{filterByExpr}, \link[edgeR]{calcNormFactors}
#'
#' @examples
#' library(Biobase)
#'
#' TRNSCRPT_EXP                   # summary
#' exprs(TRNSCRPT_EXP)[1:10, 1:5] # assayData (first 10 rows and 5 columns)
#' head(fData(TRNSCRPT_EXP))      # featureData
#' head(pData(TRNSCRPT_EXP))      # phenoData
#'
#' @keywords datasets
"TRNSCRPT_EXP"


#' @title Metabolomics ExpressionSet
#'
#' @description An object of class
#'   \code{\link[Biobase:ExpressionSet-class]{Biobase::ExpressionSet}}
#'   containing log\eqn{_2}-transformed metabolite data for differential
#'   analysis and WGCNA.
#'
#' @format Object of class
#'   \code{\link[Biobase:ExpressionSet-class]{Biobase::ExpressionSet}} with 1063
#'   features and 50 samples.
#'
#' @seealso \code{\link[Biobase]{ExpressionSet-class}},
#'   \link[MotrpacRatTraining6moData]{METAB_FEATURE_ID_MAP},
#'   \link[MotrpacRatTraining6moData]{METAB_NORM_DATA_NESTED},
#'   \href{https://www.metabolomicsworkbench.org/databases/refmet/index.php}{Metabolomics
#'   Workbench RefMet Database}
#'
#' @references Fahy, E., & Subramaniam, S. (2020). RefMet: A reference
#'   nomenclature for metabolomics. \emph{Nature Methods, 17}(12), 1173–1174.
#'   \url{https://doi.org/10.1038/s41592-020-01009-y}
#'
#' @examples
#' library(Biobase)
#'
#' METAB_EXP                   # summary
#' exprs(METAB_EXP)[1:10, 1:5] # assayData (first 10 rows and 5 columns)
#' head(fData(METAB_EXP))      # featureData
#' head(pData(METAB_EXP))      # phenoData
#'
#' @keywords datasets
"METAB_EXP"


## Differential Analysis =======================================================

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
#'   In addition to all columns present in the \code{Biobase::fData} tables,
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
#'   on the \code{ExpressionSet} objects.
#'
#' @seealso \code{\link[limma]{topTable}},
#'   \code{\link[MotrpacRatTraining6moWAT]{limma_full}},
#'   \code{\link{PROT_EXP}}, \code{\link{PHOSPHO_EXP}},
#'   \code{\link{TRNSCRPT_EXP}}, \code{\link{METAB_EXP}}
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


# Enrichment Analysis ==========================================================

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
#' @references Ashburner, M., Ball, C. A., Blake, J. A., Botstein, D., Butler,
#'   H., Cherry, J. M., Davis, A. P., Dolinski, K., Dwight, S. S., Eppig, J. T.,
#'   Harris, M. A., Hill, D. P., Issel-Tarver, L., Kasarskis, A., Lewis, S.,
#'   Matese, J. C., Richardson, J. E., Ringwald, M., Rubin, G. M., & Sherlock,
#'   G. (2000). Gene ontology: tool for the unification of biology. The Gene
#'   Ontology Consortium. \emph{Nature genetics, 25}(1), 25–29.
#'   \url{https://doi.org/10.1038/75556}
#'
#'   Dolgalev, I. msigdbr: MSigDB Gene Sets for Multiple Organisms in a Tidy
#'   Data Format. R package version 7.5.1,
#'   \url{https:://igordot.github.io/msigdbr}
#'
#'   Gene Ontology Consortium (2021). The Gene Ontology resource: enriching a
#'   GOld mine. \emph{Nucleic acids research, 49}(D1), D325–D334.
#'   \url{https://doi.org/10.1093/nar/gkaa1113}
#'
#'   Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M., Mesirov, J. P.,
#'   & Tamayo, P. (2015). The Molecular Signatures Database (MSigDB) hallmark
#'   gene set collection. \emph{Cell systems, 1}(6), 417–425.
#'   \url{https://doi.org/10.1016/j.cels.2015.12.004}
#'
#' @keywords datasets
"MSIGDB_PATHWAYS"


#' @title PhosphoSitePlus kinase-substrate relationships
#'
#' @description Kinase-substrate relationship data used to infer kinase activity
#'   from changes in phosphorylation sites.
#'
#' @format An object of class \code{data.frame} with 12684 rows and 16 columns:
#'
#' \describe{
#'   \item{GENE}{character; gene symbol.}
#'   \item{KINASE}{character; kinase.}
#'   \item{KIN_ACC_ID}{character; kinase UniProt accession.}
#'   \item{KIN_ORGANISM}{character; kinase organism. Only "human".}
#'   \item{SUBSTRATE}{character; substrate.}
#'   \item{SUB_GENE_ID}{integer; substrate Entrez gene ID.}
#'   \item{SUB_ACC_ID}{character; substrate UniProt accession.}
#'   \item{SUB_GENE}{character; substrate gene.}
#'   \item{SUB_ORGANISM}{character; substrate organism. Only "human".}
#'   \item{SUB_MOD_RSD}{character; modification site (S, T, or Y followed by
#'   amino acid position.)}
#'   \item{SITE_GRP_ID}{character; unique identifier for a modification site and
#'   its homologous sites in all proteoforms and species.}
#'   \item{SITE_+/-7_AA}{character; site flanking sequence. The amino acid for
#'   the site itself is lowercase. If the site is near the start or end of the
#'   protein sequence, the flanking sequence will be padded with underscores.}
#'   \item{DOMAIN}{character; pfam domain for the modification site.}
#'   \item{IN_VIVO_RXN}{character; in vivo kinase assay/reaction.}
#'   \item{IN_VITRO_RXN}{character; in vitro kinase assay/reaction.}
#'   \item{CST_CAT#}{character; unknown.}
#' }
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
#'   D512–D520. \url{https://doi.org/10.1093/nar/gku1267}
#'
#' @keywords datasets
"PSP_KINASE_SUBSTRATE"


#' @title Human MitoCarta3.0 Pathways
#'
#' @description Pathways from the human MitoCarta3.0 database.
#'
#' @usage MITOCARTA_HS
#'
#' @format A \code{data.frame} with 149 rows (pathways) and 3 variables:
#'
#' \describe{
#'   \item{pathway}{character; pathway description.}
#'   \item{hierarchy}{character; pathway hierarchy.}
#'   \item{human_genes}{list; character list of human genes in each pathway.}
#' }
#'
#' @source The "HumanMitoCarta3.0.xls" file was downloaded from
#'   \url{https://www.broadinstitute.org/files/shared/metabolism/mitocarta/human.mitocarta3.0.html}.
#'   Columns 2–4 were extracted from the "C MitoPathways" sheet and renamed.
#'
#' @seealso
#' \url{https://personal.broadinstitute.org/scalvo/MitoCarta3.0/human.mitopathways3.0.html}
#'
#' @references Rath, S., Sharma, R., Gupta, R., Ast, T., Chan, C., Durham, T.
#'   J., Goodman, R. P., Grabarek, Z., Haas, M. E., Hung, W. H. W., Joshi, P.
#'   R., Jourdain, A. A., Kim, S. H., Kotrys, A. V., Lam, S. S., McCoy, J. G.,
#'   Meisel, J. D., Miranda, M., Panda, A., Patgiri, A., … Mootha, V. K. (2021).
#'   MitoCarta3.0: an updated mitochondrial proteome now with sub-organelle
#'   localization and pathway annotations. \emph{Nucleic acids research,
#'   49}(D1), D1541–D1547. \url{https://doi.org/10.1093/nar/gkaa1011}
#'
#' @keywords datasets
#'
"MITOCARTA_HS"


#' @title Gene Ontology FGSEA results
#'
#' @description Fast Gene Set Enrichment Analysis (FGSEA) of transcriptomics and
#'   proteomics differential analysis results using Gene Ontology terms from the
#'   Molecular Signatures Database (MSigDB v7.5.1).
#'
#' @usage
#' TRNSCRPT_FGSEA # transcriptomics
#' PROT_FGSEA     # proteomics
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
#'   Each \code{data.frame} contains the following 12 variables:
#'
#'   \describe{
#'     \item{pathway}{character; Gene Ontology pathway identifier.}
#'     \item{gs_subcat}{character; GO:BF (Biological Processes), GO:MF
#'     (Molecular Functions), or GO:CC (Cellular Components). See
#'     \url{http://geneontology.org/docs/ontology-documentation/}}
#'     \item{gs_description}{character; Gene Ontology set description. See
#'     \code{\link[MotrpacRatTraining6moWAT]{update_GO_names}} for details.}
#'     \item{pval}{numeric; enrichment p-value.}
#'     \item{padj}{numeric; BH-adjusted p-value. P-values are adjusted across
#'     all contrasts within each ontology (GO:BP, GO:CC, GO:MF).}
#'     \item{log2err}{numeric; the expected error for the standard deviation of
#'     the p-value logarithm.}
#'     \item{ES}{numeric; gene set enrichment score.}
#'     \item{NES}{numeric; normalized enrichment score. Accounts for differences
#'     in gene set size. Calculated as the quotient of the ES to the absolute
#'     mean of the permutation enrichment scores that have the same sign as the
#'     ES. See \url{https://doi.org/10.1101/060012} for details. Values in [-1,
#'     1] are not interesting.}
#'     \item{size}{integer; number of genes in each pathway after it has been
#'     filtered to those that were present in the differential analysis
#'     results.}
#'     \item{contrast}{character; the comparison being made.}
#'     \item{leadingEdge}{list; the Entrez gene IDs that drive the enrichment.
#'     If the set has a negative ES, the genes are arranged from most to least
#'     negative; otherwise, they are arranged from most to least positive. See
#'     \url{http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_a_Leading}.}
#'     \item{leadingEdge_genes}{list; the \code{leadingEdge} after mapping
#'     Entrez genes to more human-readable gene symbols. Any genes that could
#'     not be mapped were discarded, so \code{lengths(leadingEdge_genes)} may
#'     not match \code{lengths(leadingEdge)}.}
#'   }
#'
#' @seealso \code{\link[fgsea]{fgseaMultilevel}},
#'   \code{\link[MotrpacRatTraining6moWAT]{fgsea2}},
#'   \code{\link[MotrpacRatTraining6moWAT]{msigdbr2}},
#'   \code{\link[MotrpacRatTraining6moWAT]{enrichmat}}, \code{\link{PROT_DA}},
#'   \code{\link{TRNSCRPT_DA}},
#'   \code{\link{MSIGDB_PATHWAYS}}
#'
#' @references Ashburner, M., Ball, C. A., Blake, J. A., Botstein, D., Butler, H., Cherry,
#'   J. M., Davis, A. P., Dolinski, K., Dwight, S. S., Eppig, J. T., Harris, M.
#'   A., Hill, D. P., Issel-Tarver, L., Kasarskis, A., Lewis, S., Matese, J. C.,
#'   Richardson, J. E., Ringwald, M., Rubin, G. M., & Sherlock, G. (2000). Gene
#'   ontology: tool for the unification of biology. The Gene Ontology
#'   Consortium. \emph{Nature genetics, 25}(1), 25–29.
#'   \url{https://doi.org/10.1038/75556}
#'
#'   Dolgalev, I. msigdbr: MSigDB Gene Sets for Multiple Organisms in a Tidy
#'   Data Format. R package version 7.5.1,
#'   \url{https:://igordot.github.io/msigdbr}
#'
#'   Gene Ontology Consortium (2021). The Gene Ontology resource: enriching a
#'   GOld mine. \emph{Nucleic acids research, 49}(D1), D325–D334.
#'   \url{https://doi.org/10.1093/nar/gkaa1113}
#'
#'   Korotkevich, G., Sukhov, V., Budin, N., Shpak, B., Artyomov, M. N., &
#'   Sergushichev, A. (2021). Fast gene set enrichment analysis. \emph{BioRxiv}.
#'   \url{https://doi.org/10.1101/060012}
#'
#'   Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M.,
#'   Mesirov, J. P., & Tamayo, P. (2015). The Molecular Signatures Database
#'   (MSigDB) hallmark gene set collection. \emph{Cell systems, 1}(6), 417–425.
#'   \url{https://doi.org/10.1016/j.cels.2015.12.004}
#'
#'   Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L.,
#'   Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E.
#'   S., & Mesirov, J. P. (2005). Gene set enrichment analysis: A
#'   knowledge-based approach for interpreting genome-wide expression profiles.
#'   \emph{Proceedings of the National Academy of Sciences of the United States
#'   of America, 102}(43), 15545–15550.
#'   \url{https://doi.org/10.1073/pnas.0506580102}
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
#' @description FGSEA of the proteomics differential analysis using terms from
#'   the MitoCarta3.0 database.
#'
#' @usage PROT_MITOCARTA_FGSEA
#'
#' @format A named list of 3 \code{data.frame} objects. See
#'   \code{\link{GO_FGSEA}} for details.
#'
#' @seealso \code{\link{MITOCARTA_HS}}
#'
#' @references Korotkevich, G., Sukhov, V., Budin, N., Shpak, B., Artyomov, M.
#'   N., & Sergushichev, A. (2021). Fast gene set enrichment analysis.
#'   \emph{BioRxiv}. \url{https://doi.org/10.1101/060012}
#'
#'   Rath, S., Sharma, R., Gupta, R., Ast, T., Chan, C., Durham, T. J., Goodman,
#'   R. P., Grabarek, Z., Haas, M. E., Hung, W. H. W., Joshi, P. R., Jourdain,
#'   A. A., Kim, S. H., Kotrys, A. V., Lam, S. S., McCoy, J. G., Meisel, J. D.,
#'   Miranda, M., Panda, A., … Mootha, V. K. (2020). MitoCarta3.0: An updated
#'   mitochondrial proteome now with sub-organelle localization and pathway
#'   annotations. \emph{Nucleic Acids Research, 49}(D1), D1541–D1547.
#'   \url{https://doi.org/10.1093/nar/gkaa1011}
#'
#'   Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L.,
#'   Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E.
#'   S., & Mesirov, J. P. (2005). Gene set enrichment analysis: A
#'   knowledge-based approach for interpreting genome-wide expression profiles.
#'   \emph{Proceedings of the National Academy of Sciences of the United States
#'   of America, 102}(43), 15545–15550.
#'   \url{https://doi.org/10.1073/pnas.0506580102}
#'
#' @keywords datasets
"PROT_MITOCARTA_FGSEA"


#' @title Phosphoproteomics KSEA
#'
#' @description Kinase—Substrate Enrichment Analysis (KSEA) of the
#'   phosphoproteomics differential analysis results using kinase sets from
#'   PhosphoSitePlus.
#'
#' @usage PHOSPHO_KSEA
#'
#' @format A named list of 3 \code{data.frame} objects. Similar to what is in
#'   \code{\link{GO_FGSEA}}.
#'
#' @details See \code{\link{PSP_KINASE_SUBSTRATE}} for details on the kinase
#'   sets.
#'
#'   Rat phosphorylation sites in \code{\link{PHOSPHO_DA}} were mapped to their
#'   human orthologs and multiply-phosphorylated sites were split and aggregated
#'   to the single-site level by taking the arithmetic mean of the
#'   -log\eqn{_{10}}(p-value) signed by the log\eqn{_2} fold-change. FGSEA was
#'   then performed with \code{\link[MotrpacRatTraining6moWAT]{fgsea2}} by
#'   grouping substrates according to their kinases and treating these kinase
#'   groups as gene sets. Only those kinases with at least 3 substrates were
#'   tested for enrichment.
#'
#' @references Hornbeck, P. V., Zhang, B., Murray, B., Kornhauser, J. M.,
#'   Latham, V., & Skrzypek, E. (2015). PhosphoSitePlus, 2014: mutations, PTMs
#'   and recalibrations. \emph{Nucleic acids research, 43}(Database issue),
#'   D512–D520. \url{https://doi.org/10.1093/nar/gku1267}
#'
#' @seealso \code{\link{PHOSPHO_DA}}, \code{\link{PSP_KINASE_SUBSTRATE}},
#'   \code{\link[MotrpacRatTraining6moWAT]{fgsea2}}
#'
#' @keywords datasets
"PHOSPHO_KSEA"


#' @title Metabolomics RefMet chemical subclass FGSEA
#'
#' @description FGSEA of the metabolomics differential analysis results using
#'   the Metabolomics Workbench RefMet chemical subclasses.
#'
#' @usage METAB_FGSEA
#'
#' @format A named list of 3 \code{data.frame} objects. Similar to what is in
#'   \code{\link{GO_FGSEA}}.
#'
#' @details Metabolites were grouped according to their RefMet subclass and
#'   FGSEA was performed the same as for gene sets. Only those subclasses with
#'   at least 10 metabolites, after filtering to what was detected in WAT, were
#'   tested for enrichment.
#'
#' @seealso
#'   \href{https://www.metabolomicsworkbench.org/databases/refmet/index.php}{Metabolomics
#'   Workbench RefMet Database}
#'
#' @references Fahy, E., & Subramaniam, S. (2020). RefMet: A reference
#'   nomenclature for metabolomics. \emph{Nature Methods, 17}(12), 1173–1174.
#'   \url{https://doi.org/10.1038/s41592-020-01009-y}
#'
#' @keywords datasets
"METAB_FGSEA"


## WGCNA =======================================================================

#' @title WGCNA results
#'
#' @description Weighted Gene Co-expression Network Analysis (WGCNA) applied to
#'   subcutaneous white adipose tissue proteomics, transcriptomics, and
#'   metabolomics datasets to derive clusters of similarly-behaving features for
#'   further examination.
#'
#' @usage
#' PROT_WGCNA     # proteomics
#' TRNSCRPT_WGCNA # transcriptomics
#' METAB_WGCNA    # metabolomics
#'
#' @details See \code{\link[MotrpacRatTraining6moWAT]{run_WGCNA}} for a
#'   description of the output.
#'
#' @seealso \code{\link[MotrpacRatTraining6moWAT]{run_WGCNA}}
#'
#' @references Horvath, S. (2011). \emph{Weighted Network Analysis}. Springer New
#'   York. \url{https://doi.org/10.1007/978-1-4419-8819-5}
#'
#'   Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted
#'   correlation network analysis. \emph{BMC Bioinformatics, 9}(1), 1–13.
#'   \url{https://doi.org/10.1186/1471-2105-9-559}
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


#' @title Over-representation analysis of WGCNA modules
#'
#' @description Results of performing Gene Ontology over-representation analysis
#'   (Fisher Exact / Hypergeometric test) on each of the modules produced by
#'   Weighted Gene Co-expression Network Analysis (WGCNA).
#'
#' @usage
#' TRNSCRPT_MODULE_ORA # transcriptomics
#' PROT_MODULE_ORA     # proteomics
#' METAB_MODULE_ORA    # metabolomics
#'
#' @details The same gene sets and RefMet subclasses that were used to generate
#'   the FGSEA results were also used here.
#'
#' @seealso \code{\link[MotrpacRatTraining6moWAT]{fora2}},
#'   \code{\link{WATSC_WGCNA}}, \code{\link{MSIGDB_PATHWAYS}}
#'
#' @references Ashburner, M., Ball, C. A., Blake, J. A., Botstein, D., Butler,
#'   H., Cherry, J. M., Davis, A. P., Dolinski, K., Dwight, S. S., Eppig, J. T.,
#'   Harris, M. A., Hill, D. P., Issel-Tarver, L., Kasarskis, A., Lewis, S.,
#'   Matese, J. C., Richardson, J. E., Ringwald, M., Rubin, G. M., & Sherlock,
#'   G. (2000). Gene ontology: tool for the unification of biology. The Gene
#'   Ontology Consortium. \emph{Nature genetics, 25}(1), 25–29.
#'   \url{https://doi.org/10.1038/75556}
#'
#'   Dolgalev, I. msigdbr: MSigDB Gene Sets for Multiple Organisms in a Tidy
#'   Data Format. R package version 7.5.1,
#'   \url{https:://igordot.github.io/msigdbr}
#'
#'   Gene Ontology Consortium (2021). The Gene Ontology resource: enriching a
#'   GOld mine. \emph{Nucleic acids research, 49}(D1), D325–D334.
#'   \url{https://doi.org/10.1093/nar/gkaa1113}
#'
#'   Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M., Mesirov, J. P.,
#'   & Tamayo, P. (2015). The Molecular Signatures Database (MSigDB) hallmark
#'   gene set collection. \emph{Cell systems, 1}(6), 417–425.
#'   \url{https://doi.org/10.1016/j.cels.2015.12.004}
#'
#' @keywords datasets
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


## Miscellaneous ===============================================================

#' @title scWAT Phenotype Data
#'
#' @description Phenotype data for scWAT samples from the MoTrPAC endurance
#'   exercise training study in 6-month-old-rats: animal registration and
#'   familiarization data, NMR body composition, VO\eqn{_2}max testing, training
#'   dates and duration, terminal muscle weights, and additional calculated
#'   variables.
#'
#' @usage WATSC_PHENO
#'
#' @format A list of 7 \code{data.frame} objects, each containing various
#'   phenotype measures in a longer table format.
#'
#' @details See \code{\link[MotrpacRatTraining6moData]{PHENO}} for descriptions
#'   of each of the columns. Note that some columns may have been modified
#'   slightly in order to convert the data into a longer format.
#'
#' @examples
#' # NMR data
#' head(WATSC_PHENO[["NMR"]])
#'
#' @keywords datasets
"WATSC_PHENO"


#' @title Metabolomics triglyceride concentration
#'
#' @description Standardized median concentrations of triglyceride (normalized
#'   to internal standard) detected in subcutaneous white adipose tissue.
#'
#' @usage METAB_TG_CONC
#'
#' @format An object of class \code{ExpressionSet} with 190 features
#'   (triglyceride) and 10 samples (each combination of sex and training
#'   timepoint).
#'
#' @keywords datasets
"METAB_TG_CONC"


#' @title Mitochondrial DNA qPCR
#'
#' @description Mitochondrial DNA quantitative polymerase chain reaction (qPCR)
#'   data.
#'
#' @format \code{data.frame} with 30 rows and 8 columns:
#'
#' \describe{
#'   \item{bid}{integer; unique 5 digit identifier of all samples collected for
#'   an acute test/sample collection period. All samples collected during that
#'   period will have the same BID.}
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
#' @details Quantification of mitochondrial DNA (mtDNA) was performed and
#'   described by Amar, \emph{et al.}
#'   (\url{https://doi.org/10.1101/2023.01.13.523698}). Briefly, real-time
#'   quantitative PCR was performed in duplicate for each of the scWAT samples
#'   selected for -omics analysis. The \eqn{2^{-\Delta \Delta C_T}} method
#'   (\url{https://doi.org/10.1006/meth.2001.1262}) was then applied to estimate
#'   the relative expression of the mitochondrial D-loop. Since both target
#'   (D-loop) and internal control (\eqn{\beta}-actin) were amplified in the
#'   same well, \eqn{\Delta C_T} was calculated as the mean of
#'   (\eqn{C_{T,\beta-loop} - C_{T,\beta-actin}}) for each sample. Then,
#'   \eqn{\Delta \Delta C_T} values were obtained by subtracting each
#'   \eqn{\Delta C_T} value by the mean \eqn{\Delta C_T} of the sedentary female
#'   group (the calibrator).
#'
#' @references Amar, D., Gay, N. R., Jimenez-Morales, D., Beltran, P. M. J.,
#'   Ramaker, M. E., Raja, A. N., Zhao, B., Sun, Y., Marwaha, S., Gaul, D.,
#'   Hershman, S. G., Xia, A., Lanza, I., Fernandez, F. M., Montgomery, S. B.,
#'   Hevener, A. L., Ashley, E. A., Walsh, M. J., Sparks, L. M., … The MoTrPAC
#'   Study Group (2023). The mitochondrial multi-omic response to exercise
#'   training across tissues. \emph{BioRxiv}.
#'   \url{https://doi.org/10.1101/2023.01.13.523698}
#'
#'   Livak, K. J., & Schmittgen, T. D. (2001). Analysis of Relative Gene
#'   Expression Data Using Real-Time Quantitative PCR and the 2−ΔΔCT Method.
#'   \emph{Methods, 25}(4), 402–408.
#'   \url{https://doi.org/10.1006/meth.2001.1262}
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
#'   \item{\code{plate}}{integer; 96-well plate identifier (1–4).}
#'   \item{\code{runseq}}{integer; run order (1–146).}
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
#' @md
#'
#' @keywords datasets
"PLASMA_ANALYTES"


#' @title Statistical analyses of clinical analytes in plasma
#'
#' @description Statistical analyses of most clinical analytes described in
#'   \code{\link{PLASMA_ANALYTES}}.
#'
#' @usage PLASMA_ANALYTE_STATS
#'
#' @keywords datasets
"PLASMA_ANALYTE_STATS"


#' @title Measures of adipocyte size
#'
#' @description Measures of adipocyte size: diameter, area, volume.
#'
#' @usage
#' ADIPOCYTE_SIZE
#'
#' @format A \code{data.frame} with 55825 rows and 7 variables:
#'
#' \describe{
#'   \item{bid}{integer; unique 5 digit identifier of all samples collected for
#'   an acute test/sample collection period. All samples collected during that
#'   period will have the same BID.}
#'   \item{sex}{factor; the sex of the rat with 2 levels "Female" and "Male".}
#'   \item{timepoint}{factor; exercise training group. Either "SED"
#'   (sedentary) or the number of weeks of training ("1W", "2W", "4W", "8W").}
#'   \item{diameter}{numeric; adipocyte diameter.}
#'   \item{area}{numeric; adipocyte area.}
#'   \item{volume}{volume; adipocyte volume.}
#'   \item{count}{integer; total number of adipocytes across all images for a
#'   particular BID.}
#' }
#'
#' @details Adipocyte area was calculated using CellProfiler. Diameter and
#'   volume were derived from area, under the assumption that each adipocyte was
#'   a perfect sphere.
#'
#' @keywords datasets
"ADIPOCYTE_SIZE"


#' @title Analysis of differences in adipocyte diameter
#'
#' @description Comparisons of adipocyte diameter between timepoints using a
#'   negative binomial regression model.
#'
#' @usage ADIPOCYTE_SIZE_STATS
#'
#' @format A \code{data.frame} with 40 rows and 11 variables:
#'
#' \describe{
#'   \item{contrast}{factor; the comparison being made ("4W / SED" or
#'   "8W / SED").}
#'   \item{diameter_bin}{ordered factor; adipocyte diameter binned in
#'   approximately 5 micron intervals.}
#'   \item{sex}{factor; the sex of the rat with 2 levels "Female" and "Male".}
#'   \item{ratio}{numeric; ratio of the trained adipocyte rate to the sedentary
#'   adipocyte rate. See details.}
#'   \item{SE}{numeric; standard error for the ratio.}
#'   \item{df}{integer; degrees of freedom for the test.}
#'   \item{asymp.LCL}{numeric; lower bounds of the 95{\%} confidence interval.}
#'   \item{asymp.UCL}{numeric; upper bounds of the 95{\%} confidence interval.}
#'   \item{null}{integer; value for the null hypothesis: the \code{ratio}
#'   specified by \code{contrast} is 1.}
#'   \item{z.ratio}{numeric; Z-statistic.}
#'   \item{p.value}{numeric; Dunnett p-value.}
#' }
#'
#' @details Adipocytes are binned in approximately 5 micron intervals according
#'   to their diameter (\code{diameter_bin}). The number of adipocytes in each
#'   bin per BID is tallied, as well as the total number of adipocytes per BID.
#'   The log of the total adipocytes per BID is included as an offset term in a
#'   negative binomial regression model that contains all main effects and
#'   interactions between sex, timepoint, and \code{diameter_bin} as predictors
#'   of the number of binned adipocytes. The offset term allows for modeling of
#'   the number of binned adipocytes per total number of adipocytes by BID
#'   (modeled as a rate).
#'
#'   Within each combination of sex and diameter bin, the Dunnett procedure is
#'   applied to test for differences between the trained timepoints (4W, 8W) and
#'   the sedentary group (SED).
#'
#' @seealso \code{\link{ADIPOCYTE_SIZE}}, \code{\link[MASS]{glm.nb}},
#'   \code{\link[emmeans]{emmeans}}
#'
#' @keywords datasets
"ADIPOCYTE_SIZE_STATS"

