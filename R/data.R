#' @title Pathways from MSigDB v7.5.1
#'
#' @description Reactome, Gene Ontology (Biological Process, Cellular Component,
#'   Molecular Function), and KEGG pathways from version 7.5.1 of the Molecular
#'   Signatures Database (MSigDB). Obtained via \code{\link[msigdbr]{msigdbr}}
#'   and reformatted. Results filtered to pathways with at least 10 genes.
#'
#' @format A `data.frame` with 8652 rows (pathways) and 4 variables: \describe{
#'   \item{gs_exact_source}{character; unique pathway identifier.}
#'   \item{gs_subcat}{character; the database to which the pathway belongs. One
#'   of "GO:BP", "GO:MF", "GO:CC", "CP:REACTOME", or "CP:KEGG".}
#'   \item{gs_description}{character; description of each pathway.}
#'   \item{entrez_gene}{character list; Entrez gene IDs for each pathway.} }
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
#' @usage data(msigdb_pathways)
#'
#' @author Tyler Sagendorf
#'
#' @keywords datasets
#'
#' @md
"msigdb_pathways"


#' @title scWAT Proteomics MSnSet
#'
#' @description An MSnSet object containing the subcutaneous white adipose tissue (scWAT) proteomics data provided by the MoTrPAC Bioinformatics Center.
#'
#' @format An \code{\link[MSnbase]{MSnSet-class}} object.
#'
#' @usage data(m_prot)
#'
#' @author Tyler Sagendorf
#'
#' @keywords datasets
"m_prot"

