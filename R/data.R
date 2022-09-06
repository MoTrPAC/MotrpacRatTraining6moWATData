#' @title Pathways from MSigDB v7.5.1
#'
#' @description Gene sets / pathways from the Molecular Signatures Database
#'   (MSigDB).
#'
#' @format A `data.frame` with 6552 rows (pathways) and 5 variables:
#' \describe{
#'   \item{gs_exact_source}{character; unique pathway identifier.}
#'   \item{gs_subcat}{character; the database to which the pathway belongs. One
#'   of "GO:BP", "GO:MF", "GO:CC", "CP:REACTOME", or "CP:KEGG".}
#'   \item{gs_description}{character; description of each pathway.}
#'   \item{entrez_gene}{character list; Entrez gene IDs for each pathway.}
#'   \item{num_genes}{numeric; pathway size.}
#' }
#'
#' @details Reactome, Gene Ontology (Biological Process, Cellular Component,
#'   Molecular Function), and KEGG pathways from version 7.5.1 of the Molecular
#'   Signatures Database (MSigDB). Obtained via \code{\link[msigdbr]{msigdbr}}
#'   and reformatted. Results filtered to pathways with at least 15 and no more
#'   than 500 genes. Filtering by size was done because smaller gene sets tend
#'   to be less reliable, while larger gene sets tend to be less interpretable.
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
#'   Croft, D., O'Kelly, G., Wu, G., Haw, R., Gillespie, M., Matthews, L.,
#'   Caudy, M., Garapati, P., Gopinath, G., Jassal, B., Jupe, S., Kalatskaya,
#'   I., Mahajan, S., May, B., Ndegwa, N., Schmidt, E., Shamovsky, V., Yung, C.,
#'   Birney, E., Hermjakob, H., … Stein, L. (2011). Reactome: a database of
#'   reactions, pathways and biological processes. **Nucleic acids research,
#'   39**(Database issue), D691–D697. \url{https://doi.org/10.1093/nar/gkq1018}
#'
#'   Kanehisa, M., & Goto, S. (2000). KEGG: kyoto encyclopedia of genes and
#'   genomes. *Nucleic acids research, 28*(1), 27–30.
#'   \url{https://doi.org/10.1093/nar/28.1.27}
#'
#' @usage data("msigdb_pathways")
#'
#' @keywords datasets
#'
#' @md
"msigdb_pathways"



#' @title PhosphoSitePlus Kinase-Substrate Relationship Data
#'
#' @description Kinase-substrate relationship data used to infer kinase activity
#'   from changes in phosphorylation sites.
#'
#' @usage data("KS_data")
#'
#' @details The "Kinase_Substrate_Dataset.xlsx" file, downloaded from
#'   PhosphoSitePlus on 2022-06-05
#'   (\url{https://www.phosphosite.org/staticDownloads}), was filtered to human
#'   kinases and substrates. Additionally, instances of auto-phosphorylation
#'   were removed.
#'
#' @references Hornbeck, P. V., Zhang, B., Murray, B., Kornhauser, J. M.,
#'   Latham, V., & Skrzypek, E. (2015). PhosphoSitePlus, 2014: mutations, PTMs
#'   and recalibrations. *Nucleic acids research, 43*(Database issue),
#'   D512--D520. \url{https://doi.org/10.1093/nar/gku1267}
#'
#' @md
#'
"KS_data"




#' @title scWAT Proteomics MSnSet
#'
#' @description An MSnSet object containing the subcutaneous white adipose
#'   tissue (scWAT) proteomics data provided by the MoTrPAC Bioinformatics
#'   Center.
#'
#' @format An \code{\link[MSnbase]{MSnSet-class}} object.
#'
#' @usage data("m_prot")
#'
#' @keywords datasets
"m_prot"




#' @title scWAT Phospho-proteomics MSnSet
#'
#' @description An MSnSet object containing the subcutaneous white adipose
#'   tissue (scWAT) phospho-proteomics data provided by the MoTrPAC Bioinformatics
#'   Center.
#'
#' @format An \code{\link[MSnbase]{MSnSet-class}} object.
#'
#' @usage data("m_phospho")
#'
#' @keywords datasets
"m_phospho"




#' @title scWAT Transcriptomics MSnSet
#'
#' @description An MSnSet object containing the subcutaneous white adipose
#'   tissue (scWAT) transcriptomics data provided by the MoTrPAC Bioinformatics
#'   Center.
#'
#' @usage data("m_trnscrpt")
#'
#' @format An \code{\link[MSnbase]{MSnSet-class}} object with 16443 features and
#'   48 samples.
#'
#'   # exprs
#'
#'   `exprs(m_trnscrpt)` is a count matrix that has been filtered to remove 321
#'   low-count features with \code{\link[edgeR]{filterByExpr}}.
#'
#'   # fData
#'
#'   `fData(m_trnscrpt)` is a \code{\link[base]{data.frame}} with 16443 rows and
#'   3 variables:
#'   \describe{
#'     \item{entrez_id}{character list; the Entrez gene(s) associated with each
#'     feature (Ensembl gene).}
#'     \item{gene_symbol}{character list; the gene symbol(s) associated with
#'     each feature.}
#'     \item{n}{integer; number of genes annotated to each feature.}
#'   }
#'
#'   Prior to collapsing the `entrez_id` and `gene_symbol` columns into lists,
#'   any rows where `gene_symbol` began with "LOC" were removed unless that
#'   would remove entire features. This reduces the impact of the one-to-many
#'   feature-to-gene mapping on FGSEA results.
#'
#'   # pData
#'
#'   `pData(m_trnscrpt)` is a \code{\link[base]{data.frame}} with 48 rows
#'   (samples) and 9 variables:
#'   \describe{
#'     \item{viallabel}{character; unique sample identifier.}
#'     \item{pct_globin}{numeric; percent of reads mapping to globin.}
#'     \item{rin}{numeric; RNA integrity number (RIN).}
#'     \item{pct_umi_dup}{numeric; percent of PCR duplicates as quantified with Unique Molecular Identifiers (UMIs).}
#'     \item{median_5_3_bias}{numeric; median 5'-3' bias.}
#'     \item{sex}{character; "M" for males, "F" for females.}
#'     \item{timepoint}{character; "SED" for sedentary control animals. Otherwise, the weeks of exercise training that the animals underwent ("1W", "2W", "4W", "8W").}
#'     \item{exp_group}{character; the experimental group (each combination of `sex` and `timepoint`).}
#'     \item{bid}{character; Unique 5 digit numeric identifier of all samples collected for an acute test/sample collection period. All samples collected during that period will have the same BID. Same as first 5 digits of `viallabel`.}
#'   }
#'
#' @details Two outlier samples ("90423", "90410") were removed (48 remain). Then,
#'  `pct_globin`, `rin`, `pct_umi_dup`, and `median_5_3_bias` were mean-imputed,
#'  centered, and scaled.
"m_trnscrpt"




#' @title scWAT Metabolomics MSnSet
#'
#' @description An MSnSet object containing the subcutaneous white adipose
#'   tissue (scWAT) metabolomics data provided by the MoTrPAC Bioinformatics
#'   Center.
#'
#' @format An \code{\link[MSnbase]{MSnSet-class}} object with 1063 non-redundant features and 50 samples.
#'
#' @usage data("m_metab")
#'
#' @keywords datasets
"m_metab"


