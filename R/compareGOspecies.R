#' @title Visual representation for the results of functional
#'  enrichment analysis to compare two species and a series of categories
#'
#' @name compareGOspecies
#' @description compareGOspecies function provides a simple workflow to compare results
#'  of functional enrichment analysis for two species.
#'
#'  To use this function you will need two matrices with a column which, represents the features to be compared (e.g.feature).
#'  This function will extract the unique GO terms for two matrices and it will generate a presence-absence matrix
#'  where rows will represent a combination of categories and species (e.g H.sapiens AID) and columns will represent
#'  the GO terms analyzed. Further, this function  will calculate Jaccard distances and it will provide as outputs a list with four slots:
#'   1.) A principal coordinates analysis (PCoA)
#'   2.) The Jaccard distance matrix
#'   3.) A list of shared GO terms between species
#'   4.) Finally, a list of the unique GO terms and the belonging to the respective species.
#'
#' @note Do not use "-" in the feature column. This will lead to wrong results!
#'
#' @param df1 A data frame with the results of a functional enrichment analysis for the species 1
#'  with an extra column "feature" with the features to be compared
#' @param df2 A data frame with the results of a functional enrichment analysis for the species 2
#'  with an extra column "feature" with the features to be compared
#' @param GOterm_field This is a string with the column name of the GO terms (e.g; "Functional_Category")
#' @param species1 This is a string with the species name for species 1 (e.g; "H. sapiens")
#' @param species2 This is a string with the species name for species 2 (e.g; "A. thaliana")
#' @param skipPCoA This is a boolean to indicate if the PCoA graphics can be skipped
#' @param paired_lists This is a boolean to indicate if both species have same comparable categories (gene lists).
#'  If the paired_lists is FALSE the counts will be done only for species and categories will be kept in
#'  the outcomes. Please use carefully when paired_lists = FALSE.
#' @return This function will return a list with four slots: graphics, distance shared_GO_list, and unique_GO_list
#' @examples
#'
#' #Loading example datasets
#' data(H_sapiens_compress)
#' data(A_thaliana_compress)
#' #Defining the column with the GO terms to be compared
#' GOterm_field <- "Functional_Category"
#' #Defining the species names
#' species1 <- "H. sapiens"
#' species2 <- "A. thaliana"
#'
#' #Running function
#' x <- compareGOspecies(df1=H_sapiens_compress,
#'                       df2=A_thaliana_compress,
#'                       GOterm_field=GOterm_field,
#'                       species1=species1,
#'                       species2=species2,
#'                       skipPCoA=FALSE,
#'                       paired_lists=TRUE)
#'
#' \dontrun{
#' #Displaying PCoA results
#'  x$graphics
#' # Checking shared GO terms between species
#'  print(tapply(x$shared_GO_list$feature,x$shared_GO_list$feature,length))
#'  }
#'
#' @importFrom vegan vegdist
#' @importFrom ape pcoa
#' @importFrom ggplot2 ggplot geom_polygon geom_point aes
#' @importFrom ggrepel geom_text_repel
#' @importFrom grDevices chull
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export

compareGOspecies <-
  function(df1,
           df2,
           GOterm_field,
           species1,
           species2,
           skipPCoA=FALSE,
           paired_lists=TRUE) {
    if (is.null(df1) | is.null(df2)) {
      stop("One input dataframe is absent, please add it and try again")
    }
    if (is.null(species1) | is.null(species2)) {
      stop("Missing Species name , please add it and try again")
    }

    if (species1 == species2) {
      stop("Same species name to compare, please fix")
    }
    if (is.null(paired_lists)) {
      stop("Use FALSE or TRUE")
    }
    if (is.null(skipPCoA)) {
      stop("Use FALSE or TRUE")
    }


    df1 <- df1[,c("feature",GOterm_field)]
    df2 <- df2[,c("feature",GOterm_field)]

    txtProgressBar <- NULL
    chull <- NULL
    Axis.1 <- NULL
    Axis.2 <- NULL
    species <- NULL
    grp <- NULL

    comb1 <- paste(species1, "-", unique(df1[, "feature"]))
    comb1_feat <- unique(df1[, "feature"])
    comb2 <- paste(species2, "-", unique(df2[, "feature"]))
    comb2_feat <- unique(df2[, "feature"])
    comb_all_feat <- c(comb1, comb2)

    GO_terms_unique <-
      unique(c(df1[, GOterm_field], df2[, GOterm_field]))
    mat_for_dist <-
      as.data.frame(matrix(
        nrow = length(comb1) + length(comb2),
        ncol = length(GO_terms_unique)
      ))
    row.names(mat_for_dist) <- comb_all_feat
    colnames(mat_for_dist) <- GO_terms_unique
    comb_all_feat <- strsplit(comb_all_feat, " - ")
    comb_all_feat <-
      lapply(seq_len(length(comb_all_feat)), function(i) {
        x <-
          data.frame(species = comb_all_feat[[i]][1], feature = comb_all_feat[[i]][2])
      })

    comb_all_feat <- do.call(rbind, comb_all_feat)
    #row.names(comb_all_feat)
    df1$species <- species1
    df2$species <- species2
    df_total <- rbind(df1, df2)



    message("Extracting data to calculate Jaccard distances")
    pb <-
      utils::txtProgressBar(min = 0,
                            max = nrow(mat_for_dist),
                            style = 3)

    for (i in seq_len(nrow(mat_for_dist))) {
      #message(paste(round((i/nrow(mat_for_dist)*100),2),"%"))
      utils::setTxtProgressBar(pb, i)

      for (j in seq_len(ncol(mat_for_dist))) {
        if (nrow(df_total[which(
          df_total$feature == comb_all_feat[i, 2] &
          df_total$species == comb_all_feat[i, 1] &
          df_total[, GOterm_field] == GO_terms_unique[[j]]
        ), ]) > 0) {
          mat_for_dist[i, j] <- 1
        } else {
          mat_for_dist[i, j] <- 0
        }

      }
      rm(j)
    }
    rm(i)
    close(pb)

    message("Calculating Jaccard distances")
    jacc_dist <-
      vegan::vegdist(mat_for_dist, method = "jaccard", na.rm = TRUE)

    if(sum(jacc_dist,na.rm = T)==0){
      pC_flag <- FALSE
      warning("No similarities found in the Jaccard distances. Jaccard distances matrix will be empty. Skipping PCoA")
    } else {
      pC_flag <- TRUE
      message("Valid Jaccard distances calculated")
    }

    if(isTRUE(pC_flag)){
      if(isFALSE(skipPCoA)){
        vare.mds <- ape::pcoa(jacc_dist)
        vare.mds2 <- as.data.frame(vare.mds$vectors[, 1:2])
        vare.mds2$species <-
          unlist(lapply(strsplit(row.names(vare.mds2), "-"), `[[`, 1))
        vare.mds2$species <- trimws(vare.mds2$species)

        grp.a <-
          vare.mds2[vare.mds2$species == species1,][chull(vare.mds2[vare.mds2$species ==
                                                                      species1, c("Axis.1", "Axis.2")]),]  # hull values for grp A
        grp.b <-
          vare.mds2[vare.mds2$species == species2,][chull(vare.mds2[vare.mds2$species ==
                                                                      species2, c("Axis.1", "Axis.2")]),]  # hull values for grp A
        hull.data <- rbind(grp.a, grp.b)  #combine grp.a and grp.b

        species.scores <-
          as.data.frame(vare.mds2)  #Using the scores function from vegan to extract the species scores and convert to a data.frame
        species.scores$grp <-
          unlist(lapply(strsplit(row.names(species.scores), "-"), `[[`, 2))

        message("Calculating PCoA")

        ov_plot <-
          ggplot2::ggplot() +
          ggplot2::geom_polygon(
            data = hull.data,
            ggplot2::aes(
              x = Axis.1,
              y = Axis.2,
              fill = species,
              group = species
            ),
            alpha = 0.30
          ) + # add the convex hulls
          #scale_fill_manual(values=c(species1 = "green", species2= "blue")) +
          ggplot2::geom_point(
            data = vare.mds2,
            ggplot2::aes(x = Axis.1, y = Axis.2, colour = species),
            size = 8,
            alpha = 0.5
          ) + # add the point markers
          ggrepel::geom_text_repel(
            data = species.scores,
            ggplot2::aes(x = Axis.1, y = Axis.2, label = grp),
            alpha = 0.5,
            size = 5,
            show.legend = FALSE,
            colour = 'black',
            na.rm = TRUE
            # force = 30, segment.colour = NA
          )
      } else {
        message("Skipping PCoA")
        ov_plot <- NULL
      }
    } else {
      ov_plot <- NULL
    }
    comb_feat_3 <- unique(comb_all_feat$feature)

    unique_GO_list <- list()
    shared_GO_list <- list()



    if (isFALSE(paired_lists)) {

      message("Extracting shared and unique GO terms (Assuming different categories for both species)")
      message("Unique and share lists must be used carefully!")
      #unpaired lists
      count_mat2  <-
        as.data.frame(matrix(ncol = 3, nrow = ncol(mat_for_dist)))
      colnames(count_mat2) <- c(species1, species2, "shared")
      row.names(count_mat2) <- colnames(mat_for_dist)
      count_mat2[, 1] <-
        colSums(mat_for_dist[comb_all_feat$species == species1, ])
      count_mat2[, 2] <-
        colSums(mat_for_dist[comb_all_feat$species == species2, ])

      count_mat2[, 3] <- colSums(mat_for_dist)

      ####getting counts
      for (i in seq_len(nrow(count_mat2))) {
        if (count_mat2[i, 1] > 0 & count_mat2[i, 2] == 0) {
          count_mat2[i, 1] = 1
          count_mat2[i, 3] = 0
        } else if (count_mat2[i, 1] == 0 & count_mat2[i, 2] > 0) {
          count_mat2[i, 2] = 1
          count_mat2[i, 3] = 0
        } else if (count_mat2[i, 1] > 0 & count_mat2[i, 2] > 0) {
          count_mat2[i, 1] = 1
          count_mat2[i, 2] = 1
          count_mat2[i, 3] = 1
        }
      }
      rm(i)

      #getting uniques for unpaired lists
      un_sp1 <- row.names(count_mat2)[which(count_mat2[, 1] == 1)]
      un_sp2 <- row.names(count_mat2)[which(count_mat2[, 2] == 1)]
      shared_sps <- row.names(count_mat2)[which(count_mat2[, 3] == 1)]



      un_sp1 <- lapply(seq_len(length(un_sp1)), function(k) {
        x_1 <- data.frame(feature = df1$feature[which(df1[, GOterm_field] == un_sp1[k])],
                          GO = un_sp1[k],
                          species = species1)
        return(x_1)
      })
      un_sp1 <- do.call(rbind, un_sp1)
      ############################################

      un_sp2 <- lapply(seq_len(length(un_sp2)), function(k) {
        x_1 <- data.frame(feature = df2$feature[which(df2[, GOterm_field] == un_sp2[k])],
                          GO = un_sp2[k],
                          species = species2)
        return(x_1)
      })
      un_sp2 <- do.call(rbind, un_sp2)
      ############################################

      shared_sps <- lapply(seq_len(length(shared_sps)), function(k) {
        x_1 <- data.frame(feature = df_total$feature[which(df_total[, GOterm_field] ==
                                                             shared_sps[k])],
                          GO = shared_sps[k],
                          species = "shared")
        return(x_1)
      })
      shared_sps <- do.call(rbind, shared_sps)
      ############################################
      unique_GO_list <- do.call(rbind, list(un_sp1, un_sp2))

      shared_GO_list <- shared_sps


      ############################################

    } else {

      message("Extracting shared and unique GO terms (Assuming same categories for both species)")

      pb <-
        utils::txtProgressBar(min = 0,
                              max = length(comb_feat_3),
                              style = 3)

      #paired lists
      for (i in seq_len(length(comb_feat_3))) {
        utils::setTxtProgressBar(pb, i)

        x <-
          as.numeric(row.names(comb_all_feat[comb_all_feat$feature %in% comb_feat_3[[i]], ]))
        x <- mat_for_dist[x, ]
        x_col <-
          data.frame(feature = comb_feat_3[i],
                     GO = GO_terms_unique,
                     status = colSums(x))
        x_un_col <- x_col[which(x_col$status < 2 & x_col$status > 0), ]

        x_un_list_i <- list()
        for (j in seq_len(nrow(x_un_col))) {
          x_st <- x[, which(colnames(x) == x_un_col$GO[[j]]), ]
          if (x_st[1] == 1) {
            x_go_un <-
              data.frame(feature = comb_feat_3[[i]],
                         GO = x_un_col$GO[[j]],
                         species = species1)
          } else if (x_st[2] == 1) {
            x_go_un <-
              data.frame(feature = comb_feat_3[[i]],
                         GO = x_un_col$GO[[j]],
                         species = species2)
          }

          x_un_list_i[[j]] <- x_go_un

        }
        rm(j)

        x_un_list_i <- do.call(rbind, x_un_list_i)
        unique_GO_list[[i]] <- x_un_list_i

        x_col <- x_col[which(x_col$status == 2), ]
        x_col$status <- NULL
        shared_GO_list[[i]] <- x_col
      }
      rm(i)

      close(pb)

      shared_GO_list <- do.call(rbind, shared_GO_list)
      unique_GO_list <- do.call(rbind, unique_GO_list)

      if (nrow(shared_GO_list) == 0) {
        shared_GO_list <- data.frame(feature = NA,
                                     GO = NA,
                                     species = "Shared")
      } else {
        shared_GO_list$species <- "Shared"

      }
      row.names(shared_GO_list) <- NULL
      row.names(unique_GO_list) <- NULL
    }

    if(nrow(shared_GO_list)==1){
      if(is.na(shared_GO_list$feature) & is.na(shared_GO_list$GO)){
        shared_GO_list <- NULL
        warning("NO SHARED GO TERMS!")
      }
    }

    if(nrow(unique_GO_list)==1){
      if(is.na(unique_GO_list$feature) & is.na(unique_GO_list$GO)){
        unique_GO_list <- NULL
        warning("NO UNIQUE GO TERMS, ALL GO TERMS ARE SHARED!")
      }
    }
    ov_plot_list <-
      list(
        graphics = ov_plot,
        distance = jacc_dist,
        shared_GO_list = shared_GO_list,
        unique_GO_list = unique_GO_list
      )

    if(is.null(ov_plot_list$graphics) & is.null(ov_plot_list$shared_GO_list) & is.null(ov_plot_list$unique_GO_list)){
      warning("graph_two_GOspecies function will not be usable (No shared GO terms were found).
              Please use evaluateCAT_species or evaluateGO_species instead!")
    }

    return(ov_plot_list)
  }
