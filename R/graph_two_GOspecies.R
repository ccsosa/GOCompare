#' @title Comprehensive visual and graph comparison between two species and a series of categories
#' @name graph_two_GOspecies
#' @description graph_two_GOspecies is a function to create undirected graphs to compare GO terms between two species using two options:
#' 1.) Nodes are GO terms such as biological processes and the edges are features and species are edges attributes
#' 2.) Nodes are GO terms such as biological processes and species status (e.g. A.thaliana, H. sapiens or shared) are edgess
#' @param x is a list of running the comparegOspecies species.
#' @param GOterm_field This is a string with the column name of the GO terms (e.g; "Functional_Category")
#' @param species1 This is a string with the species name for the species 1 (e.g; "H. sapiens")
#' @param species2 This is a string with the species name for the species 2 (e.g; "A. thaliana")
#' @param option  (values: 1 or 2). This option allows create either a graph where nodes are GO terms and edges are features and species are edges arributes or
#'  a graph where nodes are GO terms and edges are species belonging  (default value=2).
#' @param saveGraph logical, if \code{TRUE} the function will allow save the graph in graphml format
#' @param outdir This parameter will allow save the graph file in a folder described here (e.g: "D:").This parameter only
#'  works when saveGraph=TRUE
#' @examples
#'
#' GOterm_field <- "Functional_Category"
#' data(comparison_example)
#' #Defining the species names
#' species1 <- "H. sapiens"
#' species2 <- "A. thaliana"
#' x_graph <- graph_two_GOspecies(x=comparison_example,
#'           species1=species1,
#'           species2=species2,
#'           GOterm_field=GOterm_field,
#'           saveGraph = FALSE,
#'           option=1,
#'           outdir = NULL)
#' head(x_graph)
#' @return This function will return a table representing an edge list
#' @importFrom utils combn setTxtProgressBar txtProgressBar
#' @importFrom igraph graph_from_data_frame write.graph
#' @export

graph_two_GOspecies <-
  function(x,
           species1,
           species2,
           GOterm_field,
           saveGraph = FALSE,
           option = 2,
           outdir = NULL) {
    join_db <- rbind(x$shared_GO_list, x$unique_GO_list)
    GO_list <- unique(join_db$GO)
    unique_sp <- unique(join_db$species)

    if (is.null(option)) {
      stop("Please use a valid option")
    }

    if (isTRUE(saveGraph) & is.null(outdir)) {
      stop("Please add a valid pathway to save your graph")
    }


    if (option == 1) {
      message("Using GO terms and species as edges")
      pb <- utils::txtProgressBar(min = 0,
                                  max = length(GO_list),
                                  style = 3)



      graph_db1 <- lapply(seq_len(length(GO_list)), function(i) {
        utils::setTxtProgressBar(pb, i)


        x_shared <-
          join_db[which(join_db$GO == GO_list[[i]] &
                          join_db$species == "Shared"),]

        if (nrow(x_shared) == 1) {
          x_shared <-
            data.frame(
              SOURCE = x_shared$feature,
              TARGET = x_shared$feature,
              GO = GO_list[[i]],
              species = "Shared"
            )
        } else if (nrow(x_shared) > 1) {
          x_shared <-
            data.frame(t(utils::combn(x_shared$feature, 2)), GO = GO_list[[i]], species =
                         "Shared")
        } else {
          x_shared <-  data.frame(
            SOURCE = NA,
            TARGET = NA,
            FEATURE = NA,
            SP = NA
          )
        }
        colnames(x_shared) <- c("SOURCE", "TARGET", "FEATURE", "SP")


        x_noshared <-
          join_db[which(join_db$GO == GO_list[[i]] &
                          join_db$species != "Shared"),]

        x_noshared_1 <-
          x_noshared[which(x_noshared$species == species1),]
        x_noshared_2 <-
          x_noshared[which(x_noshared$species == species2),]


        if (nrow(x_noshared_1) == 1) {
          x_noshared_1 <-
            data.frame(
              SOURCE = x_noshared_1$feature,
              TARGET = x_noshared_1$feature,
              GO = GO_list[[i]],
              species = species1
            )
        } else if (nrow(x_noshared_1) > 1) {
          x_noshared_1 <-
            data.frame(t(utils::combn(x_noshared_1$feature, 2)), GO = GO_list[[i]], species =
                         species1)
        } else {
          x_noshared_1 <-  data.frame(
            SOURCE = NA,
            TARGET = NA,
            FEATURE = NA,
            SP = NA
          )

        }
        colnames(x_noshared_1) <-
          c("SOURCE", "TARGET", "FEATURE", "SP")

        if (nrow(x_noshared_2) == 1) {
          x_noshared_2 <-
            data.frame(
              SOURCE = x_noshared_2$feature,
              TARGET = x_noshared_2$feature,
              GO = GO_list[[i]],
              species = species2
            )
          colnames(x_noshared_2) <-
            c("SOURCE", "TARGET", "FEATURE", "SP")
        } else if (nrow(x_noshared_2) > 1) {
          x_noshared_2 <-
            data.frame(t(utils::combn(x_noshared_2$feature, 2)), GO = GO_list[[i]], species = species2)
        } else {
          x_noshared_2 <-  data.frame(
            SOURCE = NA,
            TARGET = NA,
            FEATURE = NA,
            SP = NA
          )
        }
        colnames(x_noshared_2) <-
          c("SOURCE", "TARGET", "FEATURE", "SP")


        #x_shared$i <- i ###
        #x_noshared_1$i <- i ###
        #x_noshared_2$i <- i ###



        x_final <- rbind(x_shared, x_noshared_1, x_noshared_2)
        x_final <- x_final[which(!is.na(x_final$SP)),]

        return(x_final)
      })

      close(pb)
      graph_db1 <- as.data.frame(do.call(rbind, graph_db1))
      graph_db1 <-
        graph_db1[which(graph_db1$SOURCE != "character(0)"),]
      graph_db1$SOURCE <- as.character(graph_db1$SOURCE)
      graph_db1$TARGET <- as.character(graph_db1$TARGET)
      return(graph_db1)
    } else {
      message("Using GO terms as nodes and species as edges")
      pb <-
        utils::txtProgressBar(min = 0,
                              max = length(unique_sp),
                              style = 3)
      graph_db2 <- lapply(seq_len(length(unique_sp)), function(i) {
        utils::setTxtProgressBar(pb, i)
        x <-
          unique(join_db$GO[which(join_db$species == unique_sp[[i]])])
        x <- data.frame(t(combn(x, 2)), species = unique_sp[[i]])
        colnames(x) <- c("SOURCE", "TARGET", "SP")

        return(x)
      })

      close(pb)

      graph_db2 <- do.call(rbind, graph_db2)
      return(graph_db2)
    }
    ############################################################################
    ###########################################################################
    if (isTRUE(saveGraph)) {
      if (isTRUE(option == 1)) {
        x1 <- igraph::graph_from_data_frame(graph_db1, directed = FALSE)
        igraph::write.graph(
          x1,
          file = paste0(outdir, "/", "comparison_option1.graphml"),
          format = "graphml"
        )
      } else {
        x2 <- igraph::graph_from_data_frame(graph_db2, directed = FALSE)
        igraph::write.graph(
          x2,
          file = paste0(outdir, "/", "comparison_option2.graphml"),
          format = "graphml"
        )
      }
    }
  }
