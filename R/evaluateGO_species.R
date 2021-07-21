#' @title Comprehensive comparison between species using GO terms and Pearson's Chi-squared Tests
#' @name evaluateGO_species
#' @description evaluateGO_species provides a simple workflow to compare results of functional enrichment analysis
#'  for two species through the use Pearson's Chi-squared Tests and a False discovery rate correction
#'
#' @param df1 A data frame with the results of a functional enrichment analysis for the species 1 with an extra column "feature" with
#'  the features to be compared
#' @param df2 A data frame with the results of a functional enrichment analysis for the species 2 with an extra column "feature" with
#'  the features to be compared
#' @param GOterm_field This is a string with the column name of the GO terms (e.g; "Functional.Category")
#' @param species1 This is a string with the species name for the species 1 (e.g; "H. sapiens")
#' @param species2 This is a string with the species name for the species 2 (e.g; "A. thaliana")
#' @return This function will return a data.frame with the following fields:
#' \tabular{lcc}{
#' GO \tab GO term analyzed \cr
#' pvalue \tab p-value obtained through the use of Pearson's Chi-squared Test\cr
#' FDR \tab Multiple comparison correction for the pvalue column \cr
#' }
#' @examples
#'
#' #Loading example datasets
#' data(H_sapiens)
#' data(A_thaliana)
#' #Defining the column with the GO terms to be compared
#' GOterm_field <- "Functional.Category"
#' #Defining the species names
#' species1 <- "H. sapiens"
#' species2 <- "A. thaliana"
#' #Running function
#' x <- evaluateGO_species(df1,df2, species1,species2,GOterm_field)
#' print(x)
#' @importFrom stats chisq.test p.adjust
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export

evaluateGO_species <-
  function(df1, df2, species1, species2, GOterm_field) {
    df1$species <- species1
    df2$species <- species2
    join_db <- rbind(df1, df2)
    GO_list <- unique(join_db[, GOterm_field])
    unique_sp <- c(species1, species2)
    unique_features <- unique(join_db$feature)

    if (length(unique_features) < 15) {
      warning("Number of features between two species is less than 15, be careful")
    }

    # x_tab <- table(join_db$species,join_db[,GOterm_field])
    # x_fish <- fisher.test(x_tab,simulate.p.value = TRUE,B=1000,)
    message(paste("Applying chi squared tests to ", length(GO_list), " GO terms"))
    pb <-
      utils::txtProgressBar(min = 0,
                            max = length(GO_list),
                            style = 3)

    chisq_db1 <- lapply(1:length(GO_list), function(i) {
      utils::setTxtProgressBar(pb, i)

      x <-
        join_db$species[which(join_db[, GOterm_field] == GO_list[[i]])]
      xx <- matrix(nrow = 2, ncol = 2)
      xx[1, 1] <- length(x[which(x == species1)])
      xx[2, 1] <- length(x[which(x == species2)])
      xx[1, 2] <- length(unique_features) - xx[1, 1]
      xx[2, 2] <- length(unique_features) - xx[2, 1]


      x <- suppressWarnings(stats::chisq.test(
        as.table(xx),
        correct = FALSE,
        simulate.p.value = TRUE,
        B = 1000
      ))
      df <- data.frame(GO = GO_list[[i]], pvalue = x$p.value)
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
