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

- Functions: `compareGOspecies graph_two_GOspecies evaluateCAT_species evaluateGO_species ` allow compare two species GO terms list belonging to the features needed by the user.

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

require(gprofiler2);require(stringr)

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

#Using as background the unique genes for the ten CH.
GOterm_field <- "term_name"
x_s <-  gprofiler2::gost(query = x_Hsap,
                         organism = "hsapiens", ordered_query = FALSE,
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                         measure_underrepresentation = FALSE, evcodes = FALSE,
                         user_threshold = 0.05, correction_method = "g_SCS",
                         domain_scope = "annotated", custom_bg = unique(unlist(x_Hsap)),
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
perc <- x$nodes[which(x$nodes$WEIGHT > quantile(x$nodes$WEIGHT,probs = 0.95)),]
# visualize nodes filtered
View(perc)



#########

#Running function to get graph of a list of GO terms  and categories

x_GO <- graphGOspecies(df=x_s$result,
                    GOterm_field=GOterm_field,
                    option = "GO",
                    numCores=1,
                    saveGraph=FALSE,
                    outdir = NULL,
                    filename=NULL)

# visualize nodes 
View(x_GO$nodes)

#Get GO terms nodes with values greater than 95%
perc_GO <- x_GO$nodes[which(x_GO$nodes$GO_WEIGHT > quantile(x_GO$nodes$GO_WEIGHT,probs = 0.95)),]

# visualize GO terms nodes filtered
View(perc_GO)


########################################################################################################
#two species comparison assuming they are the same genes in Drosophila melanogaster

GOterm_field <- "term_name"
x_s2 <-  gprofiler2::gost(query = x_Hsap,
                         organism = "dmelanogaster", ordered_query = FALSE,
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                         measure_underrepresentation = FALSE, evcodes = FALSE,
                         user_threshold = 0.05, correction_method = "g_SCS",
                         domain_scope = "annotated", custom_bg = unique(unlist(x_Hsap)),
                         numeric_ns = "", sources = "GO:BP", as_short_link = FALSE)

colnames(x_s2$result)[1] <- "feature"

#preparing input for compare two species
x_input <- GOCompare::compareGOspecies(x_s$result,x_s2$result,GOterm_field,species1 = "H. sapiens",species2 = "D. melanogaster",paired_lists = T)


#Comparing species results

comp_species_graph <- GOCompare::graph_two_GOspecies(x_input,species1  = "H. sapiens",species2 = "D. melanogaster",option = "Categories")

#View nodes order by combined weight (RCD category has more frequent GO terms co-occurring)
View(comp_species_graph$nodes[order(comp_species_graph$nodes$COMBINED_WEIGHT,decreasing = T),])

comp_species_graph_GO <- GOCompare::graph_two_GOspecies(x_input,species1  = "H. sapiens",species2 = "D. melanogaster",option = "GO")
#Get GO terms nodes with values greater than 95%
perc_GO_two <- comp_species_graph_GO$nodes[which(comp_species_graph_GO$nodes$GO_WEIGHT > quantile(comp_species_graph_GO$nodes$GO_WEIGHT,probs = 0.95)),]

# visualize GO terms nodes filtered and ordered (more frequent GO terms in both species and categories)

View(perc_GO_two[order(perc_GO_two$GO_WEIGHT,decreasing = T),])


#evaluating if there are different in proportions of GO terms for each category 
x_CAT <- GOCompare::evaluateCAT_species(x_s$result,x_s2$result,species1  = "H. sapiens",species2 = "D. melanogaster",GOterm_field = "term_name",test = "prop")
x_CAT <- x_CAT[which(x_CAT$FDR<=0.05),]
#View Categories with FDR <0.05 (GIM)
View(x_CAT)

#evaluating if there are different in proportions of categories for GO terms
x_GO <- GOCompare::evaluateGO_species(x_s$result,x_s2$result,species1  = "H. sapiens",species2 = "D. melanogaster",GOterm_field = "term_name",test = "prop")
x_GO <- x_GO[which(x_GO$FDR<=0.05),]
#View Categories with FDR <0.05 (No significant results in proportions)
View(x_GO)


```

## Authors
Main:Chrystian C. Sosa, Diana Carolina Clavijo-Buriticá, Mauricio Quimbaya, Victor Hugo García-Merchán

Other contributors: Nicolas Lopéz-Rozo, Camila Riccio Rengifo, David Arango Londoño, Maria Victoria Diaz

## References

Sosa, Chrystian Camilo and Clavijo-Buriticá, Diana Carolina and García-Merchán, Victor Hugo and López-Rozo, Nicolas and Riccio-Rengifo, Camila and Diaz, Maria Victoria and Londoño, David Arango and Quimbaya, Mauricio Alberto, GOCompare: An R Package to Compare Functional Enrichment Analysis between Two Species. Available at SSRN: https://ssrn.com/abstract=4201186 or http://dx.doi.org/10.2139/ssrn.4201186 
 
## License
GNU GENERAL PUBLIC LICENSE Version 3
