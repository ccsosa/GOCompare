#' @title Undirected network representation  for the results of functional
#'  enrichment analysis for one species
#' @name graphGOspecies
#' @description graphGOspecies is a function to create undirected graphs using two options:
#'  \loadmathjax
#'
#'  Categories option:
#'
#'  The nodes \mjseqn{(V)} represent groups of gene lists (categories), and the edges \mjseqn{(E)} represent GO terms co-occurring between pairs of categories. More specifically,
#'  Two categories: \mjseqn{u,v \epsilon V } are connected by an edge \mjseqn{e=(u,v)}.the edge weights \mjseqn{w(e)} are defined as the ratio of the number of GO terms co-occurring
#'  between two categories. Edge weights w(e) are defined as the ratio of the number of GO terms (e.g. biological processes) co-occurring between two categories
#'   \mjseqn{BP_{u} \ n  BP_{v}} compared to the total number of GO terms available.
#'  A node weight \mjseqn{K_{w}(u)} is defined as the sum of the edge weights where the node u is a participant. Thus, the node weight represents how frequently
#'   GO terms are reported and expressed in a biological phenomenon.
#'
#'  \mjsdeqn{w(e) = \frac{\mid BP_{u} n {BP_{v}}\mid}{\mid BP\mid}} (1)
#'
#'  \mjsdeqn{K_{w} = \sum_{{v} \epsilon {V}}{w(u,v)}} (2)
#'
#'
#'  GO option:
#'
#'  The nodes \mjseqn{{V}} represent GO terms and the edges \mjseqn{{E}'} represent categories where a pair of GO terms co-occur. More specifically,
#'  two GO terms are connected by an edge \mjseqn{{e}'=({u},{v}')}. the edge weight \mjseqn{{w}'({e}')} corresponds to the number of categories co-occurring
#'   the GO terms \mjseqn{{u}} and \mjseqn{{v}'},compared with the total number of GO terms (Equation 3). A node weight \mjseqn{{K}'_w({u}')} is defined,in this case the weight
#'  represents the importance of a GO term (more frequent co-occurring).(Please be patient, it requires a long time to finish).
#'
#'  \mjsdeqn{{w}'({e}')=\frac{\mid{Cu}'\cap  {Cv}'\mid}{\mid BP \mid}} (3)
#'
#'  \mjsdeqn{{K}'_w({u}')=\sum_{{v}'\epsilon {V}'}{{w}'({u}',{v}')}} (4)
#'
#' @param df A data frame with the results of a functional enrichment analysis for
#'  a species with an extra column "feature" with the features to be compared
#' @param GOterm_field This is a string with the column name of the GO terms (e.g: "Functional.Category")
#' @param option  (values: "GO" or "Categories"). This option allows create either a graph where nodes are GO terms and edges are features or alternatively
#'  a graph where nodes are features and edges are GO terms (default value="Categories")
#' @param numCores numeric, Number of cores to use for the process (default value numCores=2). For the example below, only one core will be used
#' @param saveGraph logical, if \code{TRUE} the function will allow save the graph in graphml format
#' @param outdir This parameter will allow save the graph file in a folder described here (e.g: "D:").This parameter only
#' works when saveGraph=TRUE
#' @param filename The name of the graph filename to be saved in the outdir detailed by the user.This parameter only
#' works when saveGraph=TRUE
#'
#' @return This function will return a list with two slots: edges and nodes.
#'
#'  (Categories):
#'   Edges list columns:
#'   \tabular{rr}{
#'   Column \tab Description \cr
#'   SOURCE and TARGET \tab The source and target categories (Nodes in the edge) \cr
#'   FEATURES_N \tab The number of GO terms between the categories \cr
#'   WEIGHT \tab Edge weight \cr
#'   FEATURES \tab GO terms available for both nodes \cr
#'   }
#'
#'   Node list columns:
#'   \tabular{rr}{
#'   Column \tab Description \cr
#'   feature\tab Category name \cr
#'   GO_count  \tab GO terms counts for the node \cr
#'   WEIGHT  \tab Node weight \cr
#'   }
#'
#'  (GO):
#'
#'   Edges list columns:
#'   \tabular{rr}{
#'   Column \tab Description \cr
#'   SOURCE and TARGET \tab The source and target GO terms (Nodes in the edge) \cr
#'   FEATURE \tab The number of Categories where both GO Terms were found \cr
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
#' @examples
#'
#' #Loading example datasets
#' data(H_sapiens_compress)
#'
#' GOterm_field <- "Functional_Category"
#'
#' #Running function
#' x <- graphGOspecies(df=H_sapiens_compress,
#'                      GOterm_field=GOterm_field,
#'                      option = "Categories",
#'                      numCores=1,
#'                      saveGraph=FALSE,
#'                      outdir = NULL,
#'                      filename=NULL)
#'
#' @importFrom utils combn
#' @importFrom igraph graph_from_data_frame write.graph
#' @importFrom parallel makeCluster parLapplyLB stopCluster detectCores
#' @importFrom stats aggregate
#' @importFrom stringr str_count
#' @import mathjaxr
#' @export

graphGOspecies <- function(df, GOterm_field, option = "Categories", numCores=2,
                           saveGraph = FALSE,outdir = NULL,filename=NULL) {
  x_det <- NULL
  x_det <- parallel::detectCores()

  df <- df[,c("feature",GOterm_field)]

  if (is.null(option) | !option %in% c("Categories","GO")) {
    stop("Please use a valid option")
  }

  if (isTRUE(saveGraph) & is.null(outdir)) {
    stop("Please add a valid path to save your graph")
  }

  if (isFALSE(saveGraph) & !is.null(outdir)) {
    stop("Please select saveGraph=TRUE to save your graph")
  }

  if (isTRUE(saveGraph) & !is.null(outdir) & is.null(filename)) {
    message(paste("Your filename will be saved as GO.graphml or CAT.graphml
            dependent of the option selected in the path:",outdir))
    message("please provide a filename if you want to avoid this behavior")
  }


  if (numCores > x_det) {
    stop("Number of cores exceed the maximum allowed by the machine,
         use a coherent number of cores such as four")
  }

  message("Obtaining features and GO terms")

  features_list <- unique(df[, "feature"])
  GO_list <- unique(df[, GOterm_field])


  if (isTRUE(option == "GO")) {
    message("Using Categories as edges and GO terms as nodes (This will take a long time)")





    cl <- parallel::makeCluster(numCores)
    parallel::clusterExport(cl, varlist=c("features_list","df","GOterm_field"),envir=environment())

    option1 <- parallel::parLapplyLB(cl,
                                    X = seq_len(length(features_list)),
                                    fun = function (i){
                                      x <- df[, GOterm_field][which(df$feature == features_list[[i]])]
                                      if (length(x) > 1) {
                                        x <- data.frame(t(utils::combn(x, 2)), GO = features_list[[i]])
                                      }
                                    })
    parallel::stopCluster(cl)

    option1 <- option1[!sapply(option1,is.null)]
    option1 <- do.call(rbind, option1)
    colnames(option1) <- c("SOURCE", "TARGET", "FEATURE")
    option1 <- option1[option1$FEATURE %in%features_list,]

    rm(cl)

    message("Extracting edge weights")

    #opt <- unique(option1[,c("SOURCE","TARGET")])

    x <-  aggregate(FEATURE ~ paste(SOURCE," @ ",TARGET)  , data = option1, length)
    x_names <- strsplit(x[,1],"@",fixed = FALSE)
    x$SOURCE <- trimws(lapply(x_names, `[[`, 1))
    x$TARGET <- trimws(lapply(x_names, `[[`, 2))
    x <- x[,c("SOURCE","TARGET","FEATURE")]
    x$WEIGHT <-x$FEATURE/length(GO_list)

    res <- x
    rm(x,option1)

    #Extracting node weights
    message("Extracting node weights")

    cl <- parallel::makeCluster(numCores)
    parallel::clusterExport(cl, varlist=c("GO_list","res"),envir=environment())
    x_att <- parallel::parLapplyLB(cl,
                                 X = seq_len(length(GO_list)),
                                 fun = function (i){
                                   x <-res[res$SOURCE %in%  trimws(GO_list[[i]]) | res$TARGET %in%  trimws(GO_list[[i]]),]

                                   x <- data.frame(GO = trimws(GO_list[[i]]),
                                                   GO_WEIGHT = sum(x$WEIGHT,na.rm = T)

                                   )

                                 })
    parallel::stopCluster(cl)

    x_att <- do.call(rbind,x_att)
    colnames(x_att) <- c("GO","GO_WEIGHT")

    #Saving in a list object
    message("Saving results in a list object")
    res <- list(nodes=x_att,edges=res)

    rm(x_att)

  } else if(option=="Categories") {


    ##option 2

    message("Using GO terms as edges and categories as nodes")

    cl <- parallel::makeCluster(numCores)

    parallel::clusterExport(cl, varlist=c("GOterm_field","GO_list","df"),envir=environment())

    option2 <- parallel::parLapplyLB(cl,
                                    X = seq_len(length(GO_list)),
                                    fun = function (i){
                                      x <- df$feature[which(df[, GOterm_field] == GO_list[[i]])]
                                      if (length(x) > 1) {

                                        x <- data.frame(t(utils::combn(x, 2)), GO = GO_list[[i]])

                                      }

                                    })

    parallel::stopCluster(cl)

    option2 <- option2[!sapply(option2,is.null)]
    option2 <- do.call(rbind, option2)
    colnames(option2) <- c("SOURCE", "TARGET", "GO")
    option2 <- option2[option2$GO %in% GO_list,]
    rm(cl)
    message("Extracting edge weights")

    opt <- unique(option2[,c("SOURCE","TARGET")])

    cl <- parallel::makeCluster(numCores)
    parallel::clusterExport(cl, varlist=c("option2","GO_list","opt"),envir=environment())

    res <- parallel::parLapplyLB(cl,
                               X = seq_len(nrow(opt)),
                               fun = function (i){
                                 x_feat <- paste(option2$GO[option2$SOURCE %in% opt[i,1] &
                                                              option2$TARGET %in% opt[i,2]],
                                                 collapse = ";")
                                 x_feat_count <- (stringr::str_count(x_feat,";")+1)
                                 x <- data.frame(SOURCE = opt$SOURCE[[i]],
                                                 TARGET = opt$TARGET[[i]],
                                                 FEATURES_N = x_feat_count,
                                                 WEIGHT = round(x_feat_count/length(GO_list),3),
                                                 FEATURES = x_feat
                                 )

                               })
    parallel::stopCluster(cl)
    res <- do.call(rbind,res)
    rm(option2)


    #Extracting node weights
    message("Extracting node weights")
    x_att <- stats::aggregate(list(df[,GOterm_field]),
                              list(df$feature),function(i){length(unique(i))})
    colnames(x_att) <- c("feature","GO_count")

    ##extract node weights based on edge weights
    #Extracting node weights

    cl <- parallel::makeCluster(numCores)
    parallel::clusterExport(cl, varlist=c("features_list","res"),envir=environment())
    x_att2 <- parallel::parLapplyLB(cl,
                                   X = seq_len(length(features_list)),
                                   fun = function (i){
                                     x <-res[res$SOURCE %in%  trimws(features_list[[i]]) | res$TARGET %in%  trimws(features_list[[i]]),]

                                     x <- data.frame(feature  = trimws(features_list[[i]]),
                                                     WEIGHT=sum(x$WEIGHT,na.rm = T)

                                     )

                                   })
    parallel::stopCluster(cl)

    x_att2 <- do.call(rbind,x_att2)
    x_att <- merge(x_att,x_att2,by = "feature",all.y = T)
    rm(x_att2)
    #Saving in a list object
    message("Saving results in a list object")
    res <- list(nodes=x_att,edges=res)

  }

  if (isTRUE(saveGraph)) {
    if (isTRUE(option == "GO")) {
      if(is.null(filename)){

        message(paste("saving ",paste0(outdir, "/", "GO.graphml")))
        x1 <- igraph::graph_from_data_frame(res$edges, directed = FALSE,vertices =  res$nodes)
        igraph::write.graph(x1,
                          file = paste0(outdir, "/", "GO.graphml"),
                          format = "graphml")
      } else {

        message(paste("saving ",paste0(outdir, "/",filename,".graphml")))
        x1 <- igraph::graph_from_data_frame(res$edges, directed = FALSE,vertices =  res$nodes)
        igraph::write.graph(x1,
                            file = paste0(outdir, "/",filename,".graphml"),
                            format = "graphml")

      }
    } else if(option =="Categories") {
      if(is.null(filename)){

        x2 <- igraph::graph_from_data_frame(res$edges, directed = FALSE,vertices = res$nodes)
        message(paste("saving ",paste0(outdir, "/", "CAT.graphml")))

        igraph::write.graph(x2,
                            file = paste0(outdir, "/", "CAT.graphml"),
                            format = "graphml")
      } else {

        x2 <- igraph::graph_from_data_frame(res$edges, directed = FALSE,vertices = res$nodes)
        message(paste("saving ",paste0(outdir, "/",filename,".graphml")))

        igraph::write.graph(x2,
                            file = paste0(outdir, "/",filename,".graphml"),
                            format = "graphml")

      }
    }
  } else {
    message("Saving as edges and node lists")
  }


  return(res)
}
