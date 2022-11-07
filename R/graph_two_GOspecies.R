#' @title Undirected network representation  for the results of functional
#'  enrichment analysis to compare two species and a series of categories
#' @name graph_two_GOspecies
#' @description graph_two_GOspecies is a function to create undirected graphs
#'  \loadmathjax
#'
#'  The graph_two_GOspecies is an analog of the graphGOspecies function, and it has the same options (" Categories " and " GO ").
#'  Nevertheless, the way in which the edge and node weights are calculated is slightly different. Since two species are compared,
#'  three possible graphs are available \({G}_1,\, {G}_2\), and \({G}_3 \). \({G}_1\), and \({G}_2 \)represent each of the species analyzed
#'  and \({G}_3\) is a subgraph of \({G}_1,\, {G}_2\), which contains the GO terms or Categories co-ocurring between both species.
#'
#'  Categories option:
#'  (Weight): The nodes (V) represent groups of gene lists (categories), and the edges (E) represent GO terms co-occurring between pairs of categories
#'  and the weight of the nodes provides a measure of how a GO term is conserved between two species and a series of categories but it is biased
#'  to categories.
#'
#'  \mjsdeqn{\widehat{K}_w(u)=\sum_{v \epsilon V_1}^{}w(u,v) + \sum_{v \epsilon V_2}^{}w(u,v)} (5)
#'
#' (shared weight): The nodes (V) represent groups of gene lists (categories), and the edges (E) represent GO terms co-occurring between pairs of categories that are only
#'  shared between species. This node weight  \({K}_s \) is computed from a shared weight of edges \({s}\), where \({N}1\) and \({N}2\) are the set of GO terms associated
#'   with the edge \(e = (u,v) \) for species 1 and 2, respectively. Therefore the node shared weight \({K}_s(u)\)  is the sum of \({s}\).
#'
#'
#'
#'  \mjsdeqn{s(e) = \frac{\mid {N1} \ n \ {N2} \mid}{\mid {N1} \bigcup {N2} \mid}} (6)
#'
#'  \mjsdeqn{{K}_s(u)=\sum_{v \epsilon (V_1  \bigcup V_2) }^{}{s(u,v)}} (7)
#'
#'
#'  (combined weight): This node weight  \({K}_c(u)\) is a combination of the weight and the shared weight. The idea of this combined weight is
#'   to find categories with more frequent GO terms co-ocurring in order to observe functional similarities between two species with a balance of GO terms co-occurring among
#'    gene lists (categories) and the two species. This node weight varies from -1 (categories with GO terms found only in one species and few categories) to 1
#'    (categories with GO terms shared widely between species and among other categories). the combined node weight \({K}_c\) is defined as the sum of the min-max normalized weights
#'     \mjseqn{\widehat{K}_w} and \({K}_s\) minus 1.
#'
#'  \mjsdeqn{minmax(y)=\frac{y-min(y)}{max(y)-min(y)}} (8)
#'  \mjsdeqn{{K}_c(u)= minmax(\widehat{K}_w(u)) + minmax({K}_s(u)) - 1 } (9)
#'
#'
#'  GO option:
#'  Given there are three possible graphs are available \({G}_1,\, {G}_2\), and \({G}_3 \). \({G}_1\), and \({G}_2 \) represent each of the species analyzed
#'  and \({G}_3\) is a subgraph of \({G}_1,\, {G}_2\), which contains the GO terms or Categories co-ocurring between both species. For this case, Nodes are GO terms and edges are
#'  categories where a GO terms is co-ocurring. This weight is similar to the GO weight calculated for graphGOspecies function. it is calculated as the equation 5.
#'
#'  \mjsdeqn{\widehat{K}_w(u)=\sum_{v \epsilon V_1}^{}w(u,v) + \sum_{v \epsilon V_2}^{}w(u,v)} (5)
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
#' @param filename The name of the graph filename to be saved in the outdir detailed by the user.This parameter only
#' works when saveGraph=TRUE
#'
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
#'           outdir = NULL,
#'           filename= NULL)
#'
#' @return This function will return a list with two slots: edges and nodes.
#'  (Categories):
#'   Edges list columns:
#'   \tabular{rr}{
#'   Column \tab Description \cr
#'   SOURCE and TARGET \tab The source and target categories (Nodes in the edge) \cr
#'   GO_N \tab The number of GO terms between the categories \cr
#'   WEIGHT \tab Edge weight \cr
#'   GO \tab GO terms available for both nodes \cr
#'   SP1 \tab Number of GO terms for the species 1 \cr
#'   SP2 \tab Number of GO terms for the species 2 \cr
#'   SHARED \tab Number of GO terms shared or co-ocurring between the categories \cr
#'   SHARED_WEIGHT \tab Shared weight for the edge \cr
#'   }
#'
#'   Node list columns:
#'   \tabular{rr}{
#'   Column \tab Description \cr
#'   CAT\tab Category name \cr
#'   CAT_WEIGHT  \tab Node weight \cr
#'   SHARED_WEIGHT  \tab Shared weight for the node \cr
#'   COMBINED_WEIGHT \tab Combined weight for the node \cr
#'   }
#'
#'  (GO):
#'
#'   Edges list columns:
#'   \tabular{rr}{
#'   Column \tab Description \cr
#'   SOURCE and TARGET \tab The source and target GO terms (Nodes in the edge) \cr
#'   FEATURE \tab The number of Categories where both GO Terms were found \cr
#'   SP \tab Species where the GO terms was found (Species 1, Species 2 or Shared) \cr
#'   WEIGHT \tab Edge weight \cr
#'   }
#'   Node list columns:
#'
#'   \tabular{rr}{
#'   Column \tab Description \cr
#'   GO \tab GO term node name \cr
#'   GO_WEIGHT \tab Node weight \cr
#'   }
#'
#' @importFrom utils combn setTxtProgressBar txtProgressBar
#' @importFrom parallel makeCluster parLapplyLB stopCluster detectCores
#' @importFrom igraph graph_from_data_frame write.graph
#' @import mathjaxr
#' @export

graph_two_GOspecies <-
  function(x,
           species1,
           species2,
           GOterm_field,
           saveGraph = FALSE,
           option = "Categories",
           numCores=2,
           outdir = NULL,
           filename = NULL) {

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


    if (isTRUE(saveGraph) & !is.null(outdir) & is.null(filename)) {
      message(paste("Your filename will be saved as GO_TWO_SP.graphml or CAT_TWO_SP.graphml
            dependent of the option selected in the path:",outdir))
      message("please provide a filename if you want to avoid this behavior")
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


      # COMBINED WEIGHT

      GO_WEIGHT_N <- (x_att$GO_WEIGHT-min(x_att$GO_WEIGHT,na.rm = T))/
          (max(x_att$GO_WEIGHT,na.rm = T)-min(x_att$GO_WEIGHT,na.rm = T))

      SHARED_WEIGHT_N <-  (x_att$SHARED_WEIGHT-min(x_att$SHARED_WEIGHT,na.rm = T))/
           (max(x_att$SHARED_WEIGHT,na.rm =T)-min(x_att$SHARED_WEIGHT,na.rm = T))


      x_att$COMBINED_WEIGHT <- (GO_WEIGHT_N + SHARED_WEIGHT_N) - 1


        #((x_att$GO_WEIGHT/max(x_att$GO_WEIGHT,na.rm = T))+
        #(x_att$SHARED_WEIGHT/max(x_att$SHARED_WEIGHT,na.rm = T)))-1
      colnames(x_att) <- c("CAT","CAT_WEIGHT","SHARED_WEIGHT","COMBINED_WEIGHT")

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
        if(is.null(filename)){

          message(paste("saving ",paste0(outdir, "/", "CAT_TWO_SP.graphml")))
          x1 <- igraph::graph_from_data_frame(res$edges, directed = FALSE,vertices =  res$nodes)
          igraph::write.graph(x1,
                              file = paste0(outdir, "/", "CAT_TWO_SP.graphml"),
                              format = "graphml")
        } else {

          message(paste("saving ",paste0(outdir, "/",filename,".graphml")))
          x1 <- igraph::graph_from_data_frame(res$edges, directed = FALSE,vertices =  res$nodes)
          igraph::write.graph(x1,
                              file = paste0(outdir, "/",filename,".graphml"),
                              format = "graphml")

        }

      } else if(isTRUE(option =="GO")){
        if(is.null(filename)){

          message(paste("saving ",paste0(outdir, "/", "GO_TWO_SP.graphml")))
          x2 <- igraph::graph_from_data_frame(res$edges, directed = FALSE,vertices =  res$nodes)
          igraph::write.graph(x2,
                              file = paste0(outdir, "/", "GO_TWO_SP.graphml"),
                              format = "graphml")
        } else {

          message(paste("saving ",paste0(outdir, "/",filename,".graphml")))
          x2 <- igraph::graph_from_data_frame(res$edges, directed = FALSE,vertices =  res$nodes)
          igraph::write.graph(x2,
                              file = paste0(outdir, "/",filename,".graphml"),
                              format = "graphml")

        }
      }
    } else {
      message("Saving as edges and nodes lists")
    }

    return(res)
  }

