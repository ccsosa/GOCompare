#' @title network representation in format graphml for the results of functional enrichment analysis for one species
#' @name graphGOspecies
#' @description graphGOspecies is a function to create undirected graphs using two options:
#' 1.) Nodes are GO terms such as biological processes and the edges are features.
#' 2.) Nodes are features  and the edges are O terms such as biological processes.
#' @param df A data frame with the results of a functional enrichment analysis for a species with an extra column "feature" with
#'  the features to be compared
#' @param GOterm_field This is a string with the column name of the GO terms (e.g; "Functional.Category")
#' @param option  (values: 1 or 2). This option allows create either a graph where nodes are GO terms and edges are features or
#' a graph where nodes are features and edges are GO terms (default value=2).
#' @param saveGraph logical, if \code{TRUE} the function will allow save the graph in graphml format
#' @param outdir This parameter will allow save the graph file in a folder described here (e.g: "D:").This parameter only
#' works when saveGraph=TRUE
#' @return This function will return a table representing an edge list
#' @examples
#'
#' #Loading example datasets
#' data(H_sapiens)
#' #Defining the column with the GO terms to be compared
#' GOterm_field <- "Functional_Category"
#' #Running function
#' x <- graphGOspecies(df=H_sapiens,
#'                      GOterm_field=GOterm_field,
#'                      option = 2,
#'                      saveGraph=FALSE,
#'                      outdir = NULL)
#' #Displaying results
#' head(x)
#' @importFrom utils combn setTxtProgressBar txtProgressBar
#' @importFrom igraph graph_from_data_frame write.graph
#' @export

graphGOspecies <- function(df, GOterm_field, option = 2, saveGraph = FALSE,outdir = NULL) {
  if (is.null(option)) {
    stop("Please use a valid option")
  }

  if (isTRUE(saveGraph) & is.null(outdir)) {
    stop("Please add a valid path to save your graph")
  }

  features_list <- unique(df[, "feature"])
  GO_list <- unique(df[, GOterm_field])


  if (isTRUE(option == 1)) {
    message("Using features as edges")
    pb <-
      utils::txtProgressBar(min = 0,
                            max = length(features_list),
                            style = 3)

    option1 <- lapply(1:length(features_list), function(i) {
      utils::setTxtProgressBar(pb, i)


      x <- df[, GOterm_field][which(df$feature == features_list[[i]])]
      if (length(x) > 1) {
        x <- data.frame(t(utils::combn(x, 2)), GO = features_list[[i]])
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
        x <- data.frame(t(utils::combn(x, 2)), GO = GO_list[[i]])
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
      x1 <- igraph::graph_from_data_frame(option1, directed = FALSE)
      igraph::write.graph(x1,
                  file = paste0(outdir, "/", "option1.graphml"),
                  format = "graphml")
    } else {
      x2 <- igraph::graph_from_data_frame(option2, directed = FALSE)
      igraph::write.graph(x2,
                  file = paste0(outdir, "/", "option2.graphml"),
                  format = "graphml")
    }
  }

}
