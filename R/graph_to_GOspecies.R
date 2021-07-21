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


GOterm_field <- "Functional_Category"
#Defining the species names
species1 <- "H. sapiens"
species2 <- "A. thaliana"



dir2 <-
  "D:/PROGRAMAS/Dropbox/shared/CANCER_HALLMARKS/GOCompare/GOCompare/data"
outdir <- "D:"
load(paste0(dir2, "/", "A_thaliana.rda"))
load(paste0(dir2, "/", "H_sapiens.rda"))

df1 <- H_sapiens
df2 <- A_thaliana


#Loading example datasets
data(H_sapiens)
data(A_thaliana)
#Defining the column with the GO terms to be compared
GOterm_field <- "Functional_Category"
#Defining the species names
species1 <- "H. sapiens"
species2 <- "A. thaliana"
#Running function


x <-
  compareGOspecies(H_sapiens,
                   A_thaliana,
                   GOterm_field,
                   species1,
                   species2)
x$shared_GO_list$species <- "Shared"
outdir <- "D:"

x <- graphGOspecies(x, GOterm_field, option==2,saveGraph = FALSE, outdir = NA)

graphGOspecies <- function(x, GOterm_field, saveGraph = FALSE, outdir = NA) {
   join_db <- rbind(x$shared_GO_list, x$unique_GO_list)
   GO_list <- unique(join_db$GO)
   unique_sp <- unique(join_db$species)

   if (is.null(option)) {
     stop("Please use a valid option")
   }

   if (isTRUE(saveGraph) & is.null(outdir)) {
     stop("Please add a valid pathway to save your graph")
   }


   if(option==1){
    message("Using GO terms and species as edges")
    pb <-
      utils::txtProgressBar(min = 0,
                            max = length(GO_list),
                            style = 3)



    graph_db1 <- lapply(1:length(GO_list), function(i) {
      utils::setTxtProgressBar(pb, i)


      x_shared <-
        join_db[which(join_db$GO == GO_list[[i]] &
                        join_db$species == "Shared"), ]

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
          data.frame(t(combn(x_shared$feature, 2)), GO = GO_list[[i]], species =
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
                        join_db$species != "Shared"), ]

      x_noshared_1 <- x_noshared[which(x_noshared$species == species1), ]
      x_noshared_2 <- x_noshared[which(x_noshared$species == species2), ]


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
          data.frame(t(combn(x_noshared_1$feature, 2)), GO = GO_list[[i]], species =
                       species1)
      } else {
        x_noshared_1 <-  data.frame(
          SOURCE = NA,
          TARGET = NA,
          FEATURE = NA,
          SP = NA
        )

      }
      colnames(x_noshared_1) <- c("SOURCE", "TARGET", "FEATURE", "SP")

      if (nrow(x_noshared_2) == 1) {
        x_noshared_2 <-
          data.frame(
            SOURCE = x_noshared_2$feature,
            TARGET = x_noshared_2$feature,
            GO = GO_list[[i]],
            species = species2
          )
        colnames(x_noshared_2) <- c("SOURCE", "TARGET", "FEATURE", "SP")
      } else if (nrow(x_noshared_2) > 1) {
        x_noshared_2 <-
          data.frame(t(combn(x_noshared_2$feature, 2)), GO = GO_list[[i]], species = species2)
      } else {
        x_noshared_2 <-  data.frame(
          SOURCE = NA,
          TARGET = NA,
          FEATURE = NA,
          SP = NA
        )
      }
      colnames(x_noshared_2) <- c("SOURCE", "TARGET", "FEATURE", "SP")


      #x_shared$i <- i ###
      #x_noshared_1$i <- i ###
      #x_noshared_2$i <- i ###



      x_final <- rbind(x_shared, x_noshared_1, x_noshared_2)
      x_final <- x_final[which(!is.na(x_final$SP)), ]

      return(x_final)
    })

    close(pb)
    graph_db1 <- as.data.frame(do.call(rbind, graph_db1))
    graph_db1 <- graph_db[which(graph_db1$SOURCE != "character(0)"), ]
    graph_db1$SOURCE <- as.character(graph_db1$SOURCE)
    graph_db1$TARGET <- as.character(graph_db1$TARGET)
    return(graph_db1)
   } else {

    message("Using GO terms as nodes and species as edges")
    pb <-
      utils::txtProgressBar(min = 0,
                            max = length(unique_sp),
                            style = 3)
    graph_db2 <- lapply(1:length(unique_sp), function(i) {
      utils::setTxtProgressBar(pb, i)
      x <- unique(join_db$GO[which(join_db$species==unique_sp[[i]])])
      x <- data.frame(t(combn(x, 2)), species = unique_sp[[i]])
      colnames(x) <- c("SOURCE", "TARGET", "SP")

      return(x)
    })

    close(pb)

    graph_db2 <- do.call(rbind,graph_db2)
    return(graph_db2)
   }
    ############################################################################
    ###########################################################################
    if (isTRUE(saveGraph)) {

      if (isTRUE(option == 1)) {
        x1 <- graph_from_data_frame(graph_db1, directed = FALSE)
        write.graph(x1,
                    file = paste0(outdir, "/", "comparison_option1.graphml"),
                    format = "graphml")
      } else {
        x2 <- graph_from_data_frame(graph_db2, directed = FALSE)
        write.graph(x2,
                    file = paste0(outdir, "/", "comparison_option2.graphml"),
                    format = "graphml")
      }

      # x1 <- igraph::graph_from_data_frame(graph_db, directed = FALSE)
      # write.table(graph_db,paste0(outdir,"/","comparison_graph.csv"),quote = FALSE,na = "",row.names = FALSE,sep = ",")
      # igraph::write.graph(
      #  x1,
      #  file = paste0(outdir, "/", "comparison_graph.graphml"),
      #  format = "graphml"
      # )
      }

  }
