# GOCompare R package v1.0.2

## Description

## Installation
GOCompare can be installed as follows
```r
#CRAN
install.packages("GOCompare")
#Alternative: GitHub
library(devtools)
remotes::install_github("ccsosa/GOCompare")
```
A full list of libraries needed for the package is included below.

**Dependencies:** `R (>= 4.0.0)`

**Imports:** `base, utils, methods, stats, grDevices, ape, vegan, ggplot2, ggrepel, igraph, parallel, stringr`

**Suggests:** `testthat`


## Usage

This R package provides five functions to provide a simple workflow to compare results of functional enrichment analysis:

- Functions: `mostFrequentGOs. graphGOspecies, evaluateGO_species` are designed to provide analysis for one species.

- Functions: `compareGOspecies graph_two_GOspecies ` allow compare two species GO terms list belonging to the features needed by the user.

- Finally, a set of four datasets for test are provided: `A_thaliana, A_thaliana_compress, H_sapiens, H_sapiens_compress, comparison_example`


### Data inputs
**Functional enrichment analyses results**

As main inputs, you will need two dataframes with the results of functional enrichment analysis from ypur favorite resource such as:
BinGO, AmiGO, ShinnyGO or TopGO.
Each file  must have this structure:

- A `data.frame` of results from a functional enrichment analysis with a column with the GO terms to be analyzed and a 
column with the category to be compared

- Depending of the function you will need to specify the species name: species1 = "H_sapiens" and species2 = "A_thaliana" 

-  A field with the column name where your GO terms to analyzed are present must be provided (e.g: GOterm_field <- "Functional_Category")


Functional_Category | feature
------------ | -------------
Response to stress | AID
Defense response | AID
Regulation of cell size  | AID
Defense response | AIM
Response to external biotic stimulus  | DCE



### Workflow
An example of how to use this package is provided below:

```r

require(gprofiler2)

url_file = "https://raw.githubusercontent.com/ccsosa/R_Examples/master/Hallmarks_of_Cancer_AT.csv"
x <- read.csv(url_file)
x[,1] <- NULL
CH <- c("AID","AIM","DCE","ERI","EGS","GIM","IA","RCD","SPS","TPI")


x_Hsap <- lapply(seq_len(length(CH)), function(i){
  x_unique <- unique(na.omit(x[,i]))
  x_unique <- x_unique[which(x_unique!="")]
  x_unique <- as.list(x_unique)
  return(x_unique)
})

names(x_Hsap) <- CH

  x_s <-  gprofiler2::gost(query = x_Hsap,
                           organism = "hsapiens", ordered_query = FALSE,
                           multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                           measure_underrepresentation = FALSE, evcodes = FALSE,
                           user_threshold = 0.05, correction_method = "g_SCS",
                           domain_scope = "annotated", custom_bg = NULL,
                           numeric_ns = "", sources = "GO:BP", as_short_link = FALSE)

colnames(x_s$result)[1] <- "feature"

#Check number of enriched terms per category
tapply(x_s$result$feature,x_s$result$feature,length)

#Running function to get graph of a list of features and GO terms

x <- graphGOspecies(df=x_s$result,
                      GOterm_field=GOterm_field,
                      option = "Categories",
                      numCores=1,
                      saveGraph=FALSE,
                      outdir = NULL,
                      filename=NULL)

# visualize nodes 
View(x$nodes)

#Get nodes with values greater than 95%
perc <- x$nodes[which(x$nodes$GO_WEIGHT > quantile(x$nodes$GO_WEIGHT,probs = 0.95)),]
# visualize nodes filtered
View(perc)
```

## Authors
Main:Chrystian C. Sosa, Diana Carolina Clavijo-Buriticá, Mauricio Quimbaya

Other contributors: Maria Victoria Diaz, Camila Riccio Rengifo, Arlen James Mosquera, Andrés Álvarez

## References

Clavijo-Buriticá, Sosa, C.C., Mosquera, A.J. Álvarez, A., Medina, J. Quimbaya, M.A. A systematic comparison of the molecular machinery associated with Cancer-Hallmarks between plants and humans reveals Arabidopsis thaliana as a useful model to understand specific carcinogenic events (to be submitted, Target journal: Plos Biology)
 
 
## License
GNU GENERAL PUBLIC LICENSE Version 3
