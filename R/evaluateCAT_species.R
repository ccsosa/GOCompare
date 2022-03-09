#' @title Comprehensive comparison between species using categories and Pearson's
#'  Chi-squared Tests
#' @name evaluateCAT_species
#' @description evaluateGO_species provides a simple function  to compare results
#'  of functional enrichment analysis for two species through the use of proportion tests or Pearson's
#'  Chi-squared Tests and a False discovery rate correction
#'
#' @param df1 A data frame with the results of a functional enrichment analysis
#'  for the species 1 with an extra column "feature" with the features to be compared
#' @param df2 A data frame with the results of a functional enrichment analysis
#'  for the species 2 with an extra column "feature" with
#'  the features to be compared
#' @param GOterm_field This is a string with the column name of the GO terms
#'  (e.g; "Functional_Category")
#' @param species1 This is a string with the species name for the species 1 (e.g; "H. sapiens")
#' @param species2 This is a string with the species name for the species 2 (e.g; "A. thaliana")
#' @param test This is a string with the hypothesis test to be performed. Two options are provided,
#' "prop" and "chi-squared" (default value="prop")
#' @return This function will return a data.frame with the following fields:
#' \tabular{lcc}{
#' CAT \tab Category  \cr
#' pvalue \tab p-value obtained through the use of Pearson's Chi-squared Test\cr
#' FDR \tab Multiple comparison correction for the p-value column \cr
#' }
#' @examples
#'
#' #Loading example datasets
#' data(H_sapiens)
#' data(A_thaliana)
#' #Defining the column with the GO terms to be compared
#' GOterm_field <- "Functional_Category"
#' #Defining the species names
#' species1 <- "H. sapiens"
#' species2 <- "A. thaliana"
#' #Running function
#' x <- evaluateCAT_species(df1= H_sapiens,
#'                        df2=A_thaliana,
#'                        species1=species1,
#'                        species2=species2,
#'                        GOterm_field=GOterm_field,
#'                          test="prop")
#' print(x)
#' @importFrom stats chisq.test p.adjust
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export

evaluateCAT_species <-
  function(df1, df2, species1, species2, GOterm_field,test="prop") {


    if (!test %in% c("prop","chi-squared")) {
      stop("Incorrect test option chosen, please use 'prop' or 'chi-squared'")
    }

    df1$species <- NA
    df2$species <- NA
    df1$species <- species1
    df2$species <- species2
    join_db <- rbind(df1, df2)
    GO_list <- unique(join_db[, GOterm_field])
    unique_sp <- c(species1, species2)
    unique_features <- unique(join_db$feature)

    if (length(GO_list) < 15) {
      warning("Number of GO terms between two species is less than 15, be careful")
    }





    if(test=="prop"){

      message(paste("Applying proportion tests to ", length(unique_features), " groups"))

      agg1 <- join_db[which(join_db$species==species1),]
      agg2 <- join_db[which(join_db$species==species2),]


      prop_db1 <- lapply(seq_len(length(unique_features)), function(i){

        x_CAT <- data.frame(CAT=unique_features[[i]],
                            species1 = length(agg1[,GOterm_field][which(agg1[,"feature"]==unique_features[i])]),
                            species2 = length(agg2[,GOterm_field][which(agg2[,"feature"]==unique_features[i])])
        )
        return(x_CAT)
      })

      prop_db1 <- do.call(rbind,prop_db1)

      ns_int <- apply(prop_db1[,-1],2,sum)

      p_values_int <- array(NA,dim(prop_db1[,-1])[1])

      pb <-
        utils::txtProgressBar(min = 0,
                              max = length(unique_features),
                              style = 3)

      for(i in seq_len(nrow(prop_db1))){
        utils::setTxtProgressBar(pb, i)
        xs <-as.numeric(prop_db1[i,-1])
        test_p <- stats::prop.test(x=xs,n =ns_int)
        p_values_int[i] <- test_p$p.value
      };rm(i)

      rm(test_p,xs)

      close(pb)

      message("Applying FDR correction")

      x_df_p_int <- data.frame(CAT = prop_db1[,1],
                               pvalue = p_values_int,
                               FDR = stats::p.adjust(p_values_int,method = "fdr"))

      x_df_p_int <- x_df_p_int[order(x_df_p_int$pvalue, decreasing = FALSE), ]

      return(x_df_p_int)

    } else if(test=="chi-squared"){

      message(paste("Applying chi squared tests to ", length(unique_features), " groups"))
      pb <-
        utils::txtProgressBar(min = 0,
                              max = length(unique_features),
                              style = 3)

      chisq_db1 <- lapply(seq_len(length(unique_features)), function(i) {
        utils::setTxtProgressBar(pb, i)

        x <-
          join_db[which(join_db[, "feature"] == unique_features[[i]]),]
        xx <- matrix(nrow = 2, ncol = 2)

        xx[1, 1] <- nrow(x[which(x$species == species1),])
        xx[2, 1] <- nrow(x[which(x$species == species2),])
        xx[1, 2] <- length(GO_list) - xx[1, 1]
        xx[2, 2] <- length(GO_list) - xx[2, 1]


        x <- suppressWarnings(stats::chisq.test(
          as.table(xx),
          correct = FALSE,
          simulate.p.value = TRUE,
          B = 1000
        ))
        df <- data.frame(CAT = unique_features[[i]], pvalue = x$p.value)
        return(df)
      })

      close(pb)

      chisq_db1 <- do.call(rbind, chisq_db1)
      chisq_db1 <- chisq_db1[which(!is.na(chisq_db1$pvalue)), ]

      message("Applying FDR correction")
      chisq_db1$FDR <- NA
      chisq_db1$FDR <-
        stats::p.adjust(chisq_db1$pvalue, method = "fdr")
      chisq_db1 <-
        chisq_db1[order(chisq_db1$pvalue, decreasing = FALSE), ]

      return(chisq_db1)
    }
  }

