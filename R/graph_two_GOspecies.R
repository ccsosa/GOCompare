#' @title Undirected network representation  for the results of functional
#'  enrichment analysis to compare two species and a series of categories
#' @name graph_two_GOspecies
#' @description graph_two_GOspecies is a function to create undirected graphs
#'  to compare GO terms between two species using two options:
#'  1.) Nodes are GO terms such as biological processes and the edges represent
#'  features for a species since the method creates a graph per species
#'  as well as shared GO terms between them.
#'  Edge weights are calculated as the intersection where cat(U) n cat(V)
#'  represents categories where the GO terms U and V are. nBP is the total number
#'  of biological processes represented by the GO terms (1). Node weights are
#'   calculated as the sum of all w(e) where the node is a participant (2)
#'    in each species and a shared GO terms(k) graphs.
#'
#'  (Please be patient, it requires a long time to finish).
#'  \deqn{w(e) = \frac{|cat(U) n cat(V)|}{|nBP|}}{%
#'  w(e) = |cat(U) n cat(V)| / |nBP| (1)}
#'
#'  \deqn{K_w(U) =  \sum(\sum(w(U,V)k=,1,k))}{%
#'  K_w(U) =  sum(sum(w(U,V),k=,1,k)) (2)}
#'
#'
#'  2.) Nodes are features and edges are GO terms available in the set of graphs (k) which consist of each species graphs and a shared GO terms graph (k).
#'  Two edges weights are calculated. First, edges weights are calculated as number of BP in the feature in comparison with the total number of GO terms available (3).
#'  Second, a shared weight is calculated for interactions shared between two species. Finally, node weights are calculated as the sum of all w(e) where the node is a participant (2) in each  species and a shared GO terms(k) graphs
#'
#'  \deqn{w(e) = \frac{|BP(U) n BP(V)|}{|nBP|}}{%
#'  w(e) = |BP(U) n BP(V)| / |nBP| (3)}
#'
#'  \deqn{K_w(U) =  \sum(\sum(w(U,V)k=,1,k))}{%
#'  K_w(U) =  sum(sum(w(U,V),k=,1,k)) (4)}
#'
#'
#' @param x is a list obtained as output of the comparegOspecies function
#' @param GOterm_field This is a string with the column name of the GO terms (e.g; "Functional_Category")
#' @param species1 This is a string with the species name for species 1 (e.g; "H. sapiens")
#' @param species2 This is a string with the species name for species 2 (e.g; "A. thaliana")
#' @param option  (values: "Categories or "GO"). This option allows create either a graph
#'  where nodes are GO terms and edges are features and GO as well as species belonging are edges attributes or
#'  a graph where nodes are GO terms and edges are species belonging  (default value="Categories")
#' @param numCores numeric, Number of cores to use for the process (default value numCores=2). For the example below, only one core will be used
#' @param saveGraph logical, if \code{TRUE} the function will allow save the graph in graphml format
#' @param outdir This parameter will allow save the graph file in a folder described here (e.g: "D:").This parameter only
#'  works when saveGraph=TRUE
#' @examples
#'
#' GOterm_field <- "Functional_Category"
#' data(comparison_ex_compress_CH)
#' #Defining the species names
#' species1 <- "H. sapiens"
#' species2 <- "A. thaliana"
#' x_graph <- graph_two_GOspecies(x=comparison_ex_compress_CH,
#'           species1=species1,
#'           species2=species2,
#'           GOterm_field=GOterm_field,
#'           numCores=1,
#'           saveGraph = FALSE,
#'           option= "Categories",
#'           outdir = NULL)
#'
#' @return This function will return a list with two slots: edges and nodes. Edges represent an edge list and their weights and
#'  nodes which represent the nodes and their respective weights (weights, shared)
#' @importFrom utils combn setTxtProgressBar txtProgressBar
#' @importFrom parallel makeCluster parLapplyLB stopCluster detectCores
#' @importFrom igraph graph_from_data_frame write.graph
#'
#' @export

graph_two_GOspecies <-
  function(x,
           species1,
           species2,
           GOterm_field,
           saveGraph = FALSE,
           option = "Categories",
           numCores=2,
           outdir = NULL) {

    x_det <- NULL
    x_det <- parallel::detectCores()

    join_db <- rbind(x$shared_GO_list, x$unique_GO_list)
    GO_list <- unique(join_db$GO)
    CAT_list <- unique(join_db$feature)
    unique_sp <- unique(join_db$species)

    if (is.null(option) | !option %in% c("Categories","GO")) {
      stop("Please use a valid option")
    }

    if (isTRUE(saveGraph) & is.null(outdir)) {
      stop("Please add a valid pathway to save your graph")
    }

    if (isFALSE(saveGraph) & !is.null(outdir)) {
      stop("Please select saveGraph=TRUE to save your graph")
    }

    if (numCores > x_det) {
      stop("Number of cores exceed the maximum allowed by the machine, use a coherent number of cores such as four")
    }

    message("Obtaining features and GO terms")


    if (option == "Categories") {

      message("Using Categories as nodes and species as edges")

      cl <- parallel::makeCluster(numCores)
      parallel::clusterExport(cl, varlist=c("GO_list","join_db","species1","species2"),envir=environment())




      graph_db1 <- parallel::parLapplyLB(cl,
                                        X = seq_len(length(GO_list)),
                                        fun = function (i){
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
                                          x_final <- rbind(x_shared, x_noshared_1, x_noshared_2)
                                          x_final <- x_final[which(!is.na(x_final$SP)),]
                                        })

      parallel::stopCluster(cl)


      graph_db1 <- as.data.frame(do.call(rbind, graph_db1))
      graph_db1 <-
      graph_db1[which(graph_db1$SOURCE != "character(0)"),]
      graph_db1 <- graph_db1[which(graph_db1$FEATURE %in% GO_list),]

      graph_db1$SOURCE <- as.character(graph_db1$SOURCE)
      graph_db1$TARGET <- as.character(graph_db1$TARGET)
      ##########################
      opt <- unique(graph_db1[,c("SOURCE","TARGET")])
      rm(cl)


      message("Extracting edge weights")

      cl <- parallel::makeCluster(numCores)
      parallel::clusterExport(cl, varlist=c("GO_list","join_db","species1","species2","graph_db1","opt"),envir=environment())

      res <- parallel::parLapplyLB(cl,
                                  X = seq_len(nrow(opt)),
                                  fun = function (i){


                                    x <- graph_db1[,c(3,4)][graph_db1$SOURCE %in% opt$SOURCE[[i]] &
                                                              graph_db1$TARGET %in% opt$TARGET[[i]],]
                                    x_feat_count <- length(unique(x$FEATURE))
                                    x_feat <- paste(x$FEATURE,
                                                    collapse = ";")
                                    x_df <- data.frame(SOURCE = opt$SOURCE[[i]],
                                                       TARGET = opt$TARGET[[i]],
                                                       FEATURES_N = x_feat_count,
                                                       WEIGHT = round(x_feat_count/length(GO_list),3),
                                                       FEATURES = x_feat,
                                                       SP1=length(x$SP[which(x$SP==species1)]),
                                                       SP2=length(x$SP[which(x$SP==species2)]),
                                                       SHARED=length(x$SP[which(x$SP=="Shared")]),
                                                       SHARED_WEIGHT=round((length(x$SP[which(x$SP=="Shared")])/x_feat_count),3))

                                  })
      parallel::stopCluster(cl)
      rm(cl)
      res <- do.call(rbind,res)
      colnames(res)[c(3,5)] <- c("GO_N","GO")


      #Extracting node weights
      message("Extracting node weights")

      cl <- parallel::makeCluster(numCores)
      parallel::clusterExport(cl, varlist=c("CAT_list","res"),envir=environment())
      x_att <- parallel::parLapplyLB(cl,
                                   X = seq_len(length(CAT_list)),
                                   fun = function (i){
                                     x <-res[res$SOURCE %in%  trimws(CAT_list[[i]]) | res$TARGET %in%  trimws(CAT_list[[i]]),]

                                     x <- data.frame(GO = trimws(CAT_list[[i]]),
                                                     GO_WEIGHT = sum(x$WEIGHT),
                                                     SHARED_WEIGHT=sum(x$SHARED_WEIGHT)

                                     )

                                   })
      parallel::stopCluster(cl)

      x_att <- do.call(rbind,x_att)

      res <- list(nodes=x_att,edges=res)


    } else if(option=="GO"){
      message("Using GO terms as nodes and species as edges")

      cl <- parallel::makeCluster(numCores)
      parallel::clusterExport(cl, varlist=c("unique_sp","join_db","species1","species2"),envir=environment())

      graph_db2 <- parallel::parLapplyLB(cl,
                                        X = seq_len(length(unique_sp)),
                                        fun = function (i){
                                          x <- join_db[which(join_db$species == unique_sp[[i]]),]
                                          x_GO <- unique(x$GO)
                                          x_feat <- unique(x$feature)

                                          x_i <- lapply(seq_len(length(x_feat)),function(j){
                                            message(j)
                                            if(length(x$GO[which(x$feature==x_feat[[j]])])>1){
                                            x <- data.frame(t(combn(x$GO[which(x$feature==x_feat[[j]])], 2)),
                                                            feat =x_feat[[j]],
                                                            species = unique_sp[[i]])
                                            } else {
                                              x <- data.frame(t(combn(rep(x$GO[which(x$feature==x_feat[[j]])],2), 2)),
                                                              feat =x_feat[[j]],
                                                              species = unique_sp[[i]])
                                            }
                                            return(x)})
                                          x_i <- do.call(rbind,x_i)
                                          colnames(x_i) <- c("SOURCE", "TARGET","FEATURE", "SP")
                                          x <-  aggregate(FEATURE ~ paste(SOURCE," @ ",TARGET)  , data = x_i, length)
                                          x_names <- strsplit(x[,1],"@",fixed = FALSE)
                                          x$SOURCE <- trimws(lapply(x_names, `[[`, 1))
                                          x$TARGET <- trimws(lapply(x_names, `[[`, 2))
                                          x <- x[,c("SOURCE","TARGET","FEATURE")]
                                          x$SP <- unique_sp[[i]]
                                          x$WEIGHT <-x$FEATURE/length(x_GO)
                                          x <- x
                                        })
      parallel::stopCluster(cl)

      graph_db2 <- graph_db2[!sapply(graph_db2,is.null)]
      graph_db2 <- do.call(rbind, graph_db2)
      graph_db2 <- graph_db2[graph_db2$SP %in%unique_sp,]

      rm(cl)

      #Extracting node weights
      message("Extracting node weights")

      cl <- parallel::makeCluster(numCores)
      parallel::clusterExport(cl, varlist=c("GO_list","graph_db2"),envir=environment())
      x_att <- parallel::parLapplyLB(cl,
                                   X = seq_len(length(GO_list)),
                                   fun = function (i){
                                     x <-graph_db2[graph_db2$SOURCE %in%  trimws(GO_list[[i]]) | graph_db2$TARGET %in%  trimws(GO_list[[i]]),]

                                     x <- data.frame(GO = trimws(GO_list[[i]]),
                                                     GO_WEIGHT = sum(x$WEIGHT)

                                     )

                                   })
      parallel::stopCluster(cl)

      x_att <- do.call(rbind,x_att)

      #Saving in a list object
      message("Saving results in a list object")
      res <- list(nodes=x_att,edges=graph_db2)



    }

    if (isTRUE(saveGraph)) {
      if (isTRUE(option == "Categories")) {
        x1 <- igraph::graph_from_data_frame(res$edges, directed = FALSE,vertices = res$nodes)
        message(paste("saving ",paste0(outdir, "/", "CAT_TWO_SP.graphml")))

        igraph::write.graph(
          x1,
          file = paste0(outdir, "/", "CAT_TWO_SP.graphml"),
          format = "graphml"
        )
      } else if(isTRUE(option =="GO")){


        x2 <- igraph::graph_from_data_frame(res$edges, directed = FALSE,vertices = res$nodes)
        message(paste("saving ",paste0(outdir, "/", "GO_TWO_SP.graphml")))
        igraph::write.graph(
          x2,
          file = paste0(outdir, "/", "GO_TWO_SP.graphml"),
          format = "graphml"
        )
      }
    }

    return(res)
  }

