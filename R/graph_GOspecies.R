#' @title network representation in format graphml for the results of functional
#'  enrichment analysis for one species
#' @name graphGOspecies
#' @description graphGOspecies is a function to create undirected graphs using two options:
#' 1.) Nodes are GO terms such as biological processes and the edges are features. Edges weights are calculated but this function
#' requires a long time to finish.
#' 2.) Nodes are features, the edges are the number of GO terms such as biological processes and they are represented
#' such as weights and collapsed in a column called FEATURES
#' @param df A data frame with the results of a functional enrichment analysis for
#'  a species with an extra column "feature" with the features to be compared
#' @param GOterm_field This is a string with the column name of the GO terms (e.g; "Functional.Category")
#' @param option  (values: 1 or 2). This option allows create either a graph where nodes are GO terms and edges are features or
#' a graph where nodes are features and edges are GO terms (default value=2).
#' @param numCores numeric, Number of cores to use for the process (default value numCores=2)
#' @param saveGraph logical, if \code{TRUE} the function will allow save the graph in graphml format
#' @param outdir This parameter will allow save the graph file in a folder described here (e.g: "D:").This parameter only
#' works when saveGraph=TRUE
#' @return This function will return a table representing an edge list
#' @examples
#'
#' #Loading example datasets
#' data(H_sapiens_compress)
#' #Defining the column with the GO terms to be compared
#' GOterm_field <- "Functional_Category"
#' #Running function
#' x <- graphGOspecies(df=H_sapiens_compress,
#'                      GOterm_field=GOterm_field,
#'                      option = 2,
#'                      numCores=2,
#'                      saveGraph=FALSE,
#'                      outdir = NULL)
#' #Displaying results
#' head(x)
#' @importFrom utils combn
#' @importFrom igraph graph_from_data_frame write.graph
#' @importFrom parallel makeCluster parLapply stopCluster detectCores
#' @importFrom stats aggregate
#' @importFrom stringr str_count
#' @export

graphGOspecies <- function(df, GOterm_field, option = 2, numCores=2,saveGraph = FALSE,outdir = NULL) {
  x_det <- NULL
  x_det <- parallel::detectCores()
  if (is.null(option) | option >3) {
    stop("Please use a valid option")
  }

  if (isTRUE(saveGraph) & is.null(outdir)) {
    stop("Please add a valid path to save your graph")
  }

  if (isFALSE(saveGraph) & !is.null(outdir)) {
    stop("Please select saveGraph=TRUE to save your graph")
  }

  if (numCores > x_det) {
    stop("Number of cores exceed the maximum allowed by the machine, use a coherent number of cores such as four")
  }

  message("Obtaining features and GO terms")

  features_list <- unique(df[, "feature"])
  GO_list <- unique(df[, GOterm_field])


  if (isTRUE(option == 1)) {
    message("Using features as edges (This will take a long time)")




    #FILTERING BY MORE THAN ONE CATEGORY
    # cl <- parallel::makeCluster(numCores)
    # parallel::clusterExport(cl, varlist=c("GO_list","df","GOterm_field"),envir=environment())
    #
    #
    #
    # GO_x <- parallel::parLapply (cl,
    #                                 X = seq_len(length(GO_list)),
    #                                 fun = function (i){
    #                                   x <- df$feature[which(df[, GOterm_field]==GO_list[[i]])]
    #                                   data.frame(GO=GO_list[[i]],FREQ=length(x))
    #                                 })
    # parallel::stopCluster(cl)
    # GO_x <- do.call(rbind,GO_x)
    # GO_list <- GO_x$GO[which(GO_x$FREQ>1)]
    # df <- df[(df[,GOterm_field] %in% GO_list),]
    #
    cl <- parallel::makeCluster(numCores)
    parallel::clusterExport(cl, varlist=c("features_list","df","GOterm_field"),envir=environment())

    option1 <- parallel::parLapply (cl,
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


  } else {


    ##option 2

    message("Using GO terms as edges")
    cl <- parallel::makeCluster(numCores)

    parallel::clusterExport(cl, varlist=c("GOterm_field","GO_list","df"),envir=environment())

    option2 <- parallel::parLapply (cl,
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

    res <- parallel::parLapply(cl,
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
  }

  if (isTRUE(saveGraph)) {
    if (isTRUE(option == 1)) {

      cl <- parallel::makeCluster(numCores)
      parallel::clusterExport(cl, varlist=c("GO_list","res"),envir=environment())
      x_att <- parallel::parLapply(cl,
                                   X = seq_len(length(GO_list)),
                                   fun = function (i){
                                     x <-res[res$SOURCE %in%  trimws(GO_list[[i]]) | res$TARGET %in%  trimws(GO_list[[i]]),]

                                     x <- data.frame(GO = trimws(GO_list[[i]]),
                                                     GO_WEIGHT = sum(x$WEIGHT)

                                     )

                                   })
      parallel::stopCluster(cl)

      x_att <- do.call(rbind,x_att)

      #x_att <- stats::aggregate(list(res$WEIGHT),list(unique(c(res$SOURCE,res$TARGET))),sum)
      colnames(x_att) <- c("GO","GO_WEIGHT")
      message(paste("saving ",paste0(outdir, "/", "option1.graphml")))
      x1 <- igraph::graph_from_data_frame(res, directed = FALSE,vertices = x_att)
      igraph::write.graph(x1,
                          file = paste0(outdir, "/", "option1.graphml"),
                          format = "graphml")
    } else {
      x_att <- stats::aggregate(list(df[,GOterm_field]),list(df$feature),function(i){length(unique(i))})
      colnames(x_att) <- c("feature","GO_count")
      x2 <- igraph::graph_from_data_frame(res, directed = FALSE,vertices = x_att)
      message(paste("saving ",paste0(outdir, "/", "option2.graphml")))

      igraph::write.graph(x2,
                          file = paste0(outdir, "/", "option2.graphml"),
                          format = "graphml")
    }
  }

  return(res)
}
