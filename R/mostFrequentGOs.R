#' @title Most frequent GO terms among groups for a data.frame
#' @name mostFrequentGOs
#' @description Provides an easy way to get the frequency of GO terms
#'  such as biological processes for a dataframe and a series of features.
#' @param df A data frame with the results of a functional enrichment analysis
#'  for a species with an extra column "feature" with the features to be compared
#' @param GOterm_field This is a string with the column name of the GO terms (e.g; "Functional.Category")
#' @return This function will return a table with the frequency of GO terms among a series of features
#' @examples
#'
#' #Loading example datasets
#' data(H_sapiens)
#' #Defining the column with the GO terms to be compared
#' GOterm_field <- "Functional_Category"
#' #Running function
#' x <- mostFrequentGOs(df=H_sapiens, GOterm_field=GOterm_field)
#' #Displaying results
#' head(x)
#' @export

mostFrequentGOs <- function(df, GOterm_field) {
   features_list <- unique(df[, "feature"])
   GO_list <- unique(df[, GOterm_field])
   x_freq <- lapply(seq_len(length(GO_list)), function(i) {
      x <- df$feature[which(df[, GOterm_field]==GO_list[[i]])]
      x <- data.frame(GO=GO_list[[i]],freq=length(x),features=paste(x,collapse = ";"))
      return(x)
   })
   x_freq <- do.call(rbind,x_freq)
   x_freq <- x_freq[order(x_freq$freq, decreasing = TRUE), ]
   return(x_freq)
}
