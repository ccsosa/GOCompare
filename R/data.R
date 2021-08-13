#' @title A thaliana functional enrichment analysis results for
#'  "AID","DCE","RCD","SPS" cancer-hallmarks
#' @name A_thaliana_compress
#' @docType data
#' @description This dataset is a subset of the original dataset obtained
#'  for Clavijo-Buriticá (In preparation)
#' @format A data frame with 2000 rows and 6 variables:
#' \describe{
#'   \item{Enrichment_FDR}{Numeric: False discovery rate values for the GO term}
#'   \item{Genes_in_list}{numeric: Number of genes in the list of genes for a given GO term}
#'   \item{Total_genes}{numeric: Number of genes in the genome of a species for a given GO term}
#'   \item{Functional_Category}{character: GO term name or GO term id}
#'   \item{Genes}{character: Genes found fot a given GO term}
#'   \item{feature}{character: A column representing the belonging of a group of comparison}
#' }
#' @references
#' Clavijo-Buriticá, Sosa, C.C., Mosquera, A.J. Álvarez, A., Medina, J. Quimbaya,
#' M.A. A systematic comparison of the molecular machinery associated with Cancer-Hallmarks
#' between plants and humans reveals Arabidopsis thaliana as a useful model to understand
#' specific carcinogenic events (to be submitted, Target journal: Plos Biology)
#' @source \url{https://data.mendeley.com/datasets/myyy2wxd59/1}
"A_thaliana_compress"

#' @title H. sapiens functional enrichment analysis results for
#'  "AID","DCE","RCD","SPS" cancer-hallmarks
#' @name H_sapiens_compress
#' @docType data
#' @description This dataset is a subset of the original dataset obtained
#' for Clavijo-Buriticá (In preparation)
#' @format A data frame with 2000 rows and 6 variables:
#' \describe{
#'   \item{Enrichment_FDR}{Numeric: False discovery rate values for the GO term}
#'   \item{Genes_in_list}{numeric: Number of genes in the list of genes for a given GO term}
#'   \item{Total_genes}{numeric: Number of genes in the genome of a species for a given GO term}
#'   \item{Functional_Category}{character: GO term name or GO term id}
#'   \item{Genes}{character: Genes found fot a given GO term}
#'   \item{feature}{character: A column representing the belonging of a group of comparison}
#' }
#' @references
#' Clavijo-Buriticá, Sosa, C.C., Mosquera, A.J. Álvarez, A., Medina, J. Quimbaya,
#' M.A. A systematic comparison of the molecular machinery associated with Cancer-Hallmarks
#' between plants and humans reveals Arabidopsis thaliana as a useful model to understand
#' specific carcinogenic events (to be submitted, Target journal: Plos Biology)
#' @source \url{https://data.mendeley.com/datasets/myyy2wxd59/1}
"H_sapiens_compress"

#' @title A thaliana functional enrichment analysis  of 2224 ortholog genes related
#'  to cancer-hallmarks
#' @name A_thaliana
#' @docType data
#' @description This dataset is the original dataset obtained for Clavijo-Buriticá (In preparation)
#' @format A data frame with 4063 rows and 6 variables:
#' \describe{
#'   \item{Enrichment_FDR}{Numeric: False discovery rate values for the GO term}
#'   \item{Genes_in_list}{numeric: Number of genes in the list of genes for a given GO term}
#'   \item{Total_genes}{numeric: Number of genes in the genome of a species for a given GO term}
#'   \item{Functional_Category}{character: GO term name or GO term id}
#'   \item{Genes}{character: Genes found fot a given GO term}
#'   \item{feature}{character: A column representing the belonging of a group of comparison}
#' }
#' @references
#' Clavijo-Buriticá, Sosa, C.C., Mosquera, A.J. Álvarez, A., Medina, J. Quimbaya,
#' M.A. A systematic comparison of the molecular machinery associated with Cancer-Hallmarks
#' between plants and humans reveals Arabidopsis thaliana as a useful model to understand
#' specific carcinogenic events (to be submitted, Target journal: Plos Biology)
#' @source \url{https://data.mendeley.com/datasets/myyy2wxd59/1}
"A_thaliana"

#' @title H. sapiens functional enrichment analysis of 5494 genes related to cancer-hallmarks
#' @name H_sapiens
#' @docType data
#' @description This dataset is a subset of the original dataset obtained for Clavijo-Buriticá (In preparation)
#' @format A data frame with 2000 rows and 6 variables:
#' \describe{
#'   \item{Enrichment_FDR}{Numeric: False discovery rate values for the GO term}
#'   \item{Genes_in_list}{numeric: Number of genes in the list of genes for a given GO term}
#'   \item{Total_genes}{numeric: Number of genes in the genome of a species for a given GO term}
#'   \item{Functional_Category}{character: GO term name or GO term id}
#'   \item{Genes}{character: Genes found fot a given GO term}
#'   \item{feature}{character: A column representing the belonging of a group of comparison}
#' }
#' @references
#' Clavijo-Buriticá, Sosa, C.C., Mosquera, A.J. Álvarez, A., Medina, J. Quimbaya,
#' M.A. A systematic comparison of the molecular machinery associated with Cancer-Hallmarks
#' between plants and humans reveals Arabidopsis thaliana as a useful model to understand
#' specific carcinogenic events (to be submitted, Target journal: Plos Biology)
#' @source \url{https://data.mendeley.com/datasets/myyy2wxd59/1}
"H_sapiens"

#' @title Functional enrichment analysis comparison between H. sapiens and A. thaliana
#' for "AID","DCE","RCD","SPS" cancer-hallmarks
#' @name comparison_ex_compress
#' @docType data
#' @description This dataset is the results of running the compareGOspecies species and
#'  it is composed of four slots:
#' \describe{
#'   \item{graphics}{PCoA graphics}
#'   \item{distance}{numeric: Jaccard distance matrix }
#'   \item{shared_GO_list}{data.frame with shared GO terms between species}
#'   \item{unique_GO_list}{data.frame with unique GO terms and their belonging two each species}
#' }
#' @references
#' Clavijo-Buriticá, Sosa, C.C., Mosquera, A.J. Álvarez, A., Medina, J. Quimbaya,
#' M.A. A systematic comparison of the molecular machinery associated with Cancer-Hallmarks
#' between plants and humans reveals Arabidopsis thaliana as a useful model to understand
#' specific carcinogenic events (to be submitted, Target journal: Plos Biology)
#' @source \url{https://data.mendeley.com/datasets/myyy2wxd59/1}
"comparison_ex_compress"
