#' @title Comprehensive visual and graph comparison between two species and a series of categories
#' @name compareGOspecies
#' @description compareGOspecies provides a simple workflow to compare  results of functional enrichment analysis
#'  for two species and a column which represent features to be compared.
#'  This function will extract the unique GO terms for two matrices matrices provided and will generate a presence-
#'  absence matrix where rows will represent a combination of categories and species(e.g H.sapies AID) and columns will be
#'  the GO terms analyzed.Further, this code will calculate Jaccard distances and will provide as outputs a principal
#'  coordinates analysis (PCoA), the jaccard distance matrix, the list of shared GO terms between species and finally,
#'  a list of the unique GO terms and the belongig to the species.
#' @param df1 A data frame with the results of a functional enrichment analysis for the species 1 with an extra column "feature" with
#'  the features to be compared
#' @param df2 A data frame with the results of a functional enrichment analysis for the species 2 with an extra column "feature" with
#'  the features to be compared
#' @param GOterm_field This is a string with the column name of the GO terms (e.g; "Functional.Category")
#' @param species1 This is a string with the species name for the species 1 (e.g; "H. sapiens")
#' @param species2 This is a string with the species name for the species 2 (e.g; "A. thaliana")
#' @param GOterm_field This is a string with the column name of the GO terms (e.g; "Functional.Category")
#' @return This function will return a list with four slots: graphic, distance shared_GO_list,and unique_GO_list
#' @examples
#'
#' #Loading example datasets
#' data(H_sapiens_compress)
#' data(A_thaliana_compress)
#' #Defining the column with the GO terms to be compared
#' GOterm_field <- "Functional.Category"
#' #Defining the species names
#' species1 <- "H. sapiens"
#' species2 <- "A. thaliana"
#' #Running function
#' x <- compareGOspecies(H_sapiens_compress,A_thaliana_compress,GOterm_field,species1,species2)
#' #Displaying PCoA results
#' x$graphic
#' # Checking shared GO terms between species
#' print(tapply(x$shared_GO_list$feature,x$shared_GO_list$feature,length))
#' # Checking unique GO terms for each species
#' print(tapply(x$unique_GO_list$species,x$unique_GO_list$species,length))
#' @export
#' @importFrom vegan vegdist
#' @importFrom ape pcoa
#' @importFrom ggplot2 ggplot geom_polygon geom_point
#' @importFrom ggrepel geom_text_repel


dir2 <-
  "D:/PROGRAMAS/Dropbox/shared/CANCER_HALLMARKS/GOCompare/GOCompare/data"
outdir <- "D:"
#load(paste0(dir2,"/","A_thaliana.rda"))
load(paste0(dir2, "/", "H_sapiens.rda"))


df <- H_sapiens
GOterm_field <- "Functional_Category"


graphGOspecies <- function(df, GOterm_field, option = 2, saveGraph = FALSE,outdir = NULL) {
  if (is.null(option)) {
    stop("Please use a valid option")
  }

  if (isTRUE(saveGraph) & is.null(outdir)) {
    stop("Please add a valid pathway to save your graph")
  }

  features_list <- unique(df[, "feature"])
  GO_list <- unique(df[, GOterm_field])


  if (isTRUE(option == 1)) {
    message("Using features as edge")
    pb <-
      utils::txtProgressBar(min = 0,
                            max = length(features_list),
                            style = 3)

    option1 <- lapply(1:length(features_list), function(i) {
      utils::setTxtProgressBar(pb, i)


      x <- df[, GOterm_field][which(df$feature == features_list[[i]])]
      if (length(x) > 1) {
        x <- data.frame(t(combn(x, 2)), GO = features_list[[i]])
        colnames(x) <- c("SOURCE", "TARGET", "GO")
      } else {
        message(paste("No records for", features_list[[i]]))
      }

      return(x)
    })

    close(pb)
    option1 <- do.call(rbind, option1)
    return(option1)
  } else {
    ##option 2

    message("Using GO terms as edges")
    pb <-
      utils::txtProgressBar(min = 0,
                            max = length(GO_list),
                            style = 3)

    option2 <- lapply(1:length(GO_list), function(i) {
      utils::setTxtProgressBar(pb, i)


      x <- df$feature[which(df[, GOterm_field] == GO_list[[i]])]
      if (length(x) > 1) {
        x <- data.frame(t(combn(x, 2)), GO = GO_list[[i]])
        colnames(x) <- c("SOURCE", "TARGET", "FEATURE")
}
      return(x)
    })

    close(pb)
    option2 <- do.call(rbind, option2)
    return(option2)
  }

  if (isTRUE(saveGraph)) {
    if (isTRUE(option == 1)) {
      x1 <- graph_from_data_frame(option1, directed = FALSE)
      write.graph(x1,
                  file = paste0(outdir, "/", "option1.graphml"),
                  format = "graphml")
    } else {
      x2 <- graph_from_data_frame(option2, directed = FALSE)
      write.graph(x2,
                  file = paste0(outdir, "/", "option2.graphml"),
                  format = "graphml")
    }
  }

}



#x <- graphGOspecies(df, GOterm_field, option = 2, saveGraph=FALSE,outdir = NULL)
