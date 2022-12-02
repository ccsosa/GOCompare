# GOCompare R package v1.0.2

## Description:

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

This R package provides six functions to provide a simple workflow to compare results of functional enrichment analysis:

- Functions: `mostFrequentGOs. graphGOspecies` are designed to provide analysis for one species.

- Functions: `compareGOspecies graph_two_GOspecies evaluateCAT_species evaluateGO_species ` allow compare two species GO terms list belonging to the categories needed by the user.

- Finally, a set of four datasets for test are provided in the package: `A_thaliana, A_thaliana_compress, H_sapiens, H_sapiens_compress, comparison_example`


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
An example of how to use this package is provided below (We will compare cancer genes for 10 hallmarks associated to carcinogenesis process and their possible orthologues in D. melanogaster as an example. Please use an accurate program or platform to obtain orthologue genes. In this case we provide an example using the `gorth` function of gprofiler2 R package.

```r
require(gprofiler2);require(stringr);require(GOCompare)

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


orth_genes <- gprofiler2::gorth(query=unique(unlist(x_Hsap)),source_organism = "hsapiens",target_organism = "dmelanogaster")

#assigning genes

x_Dmap <- list()
for(i in 1:length(x_Hsap)){
  
  D_list <- list()
  for(j in 1:length(x_Hsap[[i]])){
    x_orth <- orth_genes[orth_genes$input==x_Hsap[[i]][j],]
    if(nrow(x_orth)>0){
      D_list[[j]] <- data.frame(orth=x_orth$ortholog_name)
    } else {
      D_list[[j]] <- NULL
    }
    rm(x_orth)
  };rm(j)

  D_list <- unique(do.call(rbind,D_list))
  D_list <- D_list[which(!is.null(D_list))]
 x_Dmap[[i]] <- D_list
 rm(D_list)
};rm(i)

names(x_Dmap) <- CH


GOterm_field <- "term_name"
x_s2 <-  gprofiler2::gost(query = x_Dmap,
                          organism = "dmelanogaster", ordered_query = FALSE,
                          multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                          measure_underrepresentation = FALSE, evcodes = FALSE,
                          user_threshold = 0.05, correction_method = "g_SCS",
                          domain_scope = "annotated", custom_bg = unique(unlist(x_Dmap)),
                          numeric_ns = "", sources = "GO:BP", as_short_link = FALSE)

colnames(x_s2$result)[1] <- "feature"

#preparing input for compare two species
x_input <- GOCompare::compareGOspecies(x_s$result,x_s2$result,GOterm_field,species1 = "H. sapiens",species2 = "D. melanogaster",paired_lists = T)

#try to test similarities using clustering

plot(hclust(x_input$distance,method = "ward.D"))

#Comparing species results

comp_species_graph <- GOCompare::graph_two_GOspecies(x_input,species1  = "H. sapiens",species2 = "D. melanogaster",option = "Categories")

#View nodes order by combined weight (SPS and GIM categories have more frequent GO terms co-occurring)
View(comp_species_graph$nodes[order(comp_species_graph$nodes$COMBINED_WEIGHT,decreasing = T),])

comp_species_graph_GO <- GOCompare::graph_two_GOspecies(x_input,species1  = "H. sapiens",species2 = "D. melanogaster",option = "GO")
#Get GO terms nodes with values greater than 95%
perc_GO_two <- comp_species_graph_GO$nodes[which(comp_species_graph_GO$nodes$GO_WEIGHT > quantile(comp_species_graph_GO$nodes$GO_WEIGHT,probs = 0.95)),]

# visualize GO terms nodes filtered and ordered (more frequent GO terms in both species and categories)

View(perc_GO_two[order(perc_GO_two$GO_WEIGHT,decreasing = T),])


#evaluating if there are different in proportions of GO terms for each category 
x_CAT <- GOCompare::evaluateCAT_species(x_s$result,x_s2$result,species1  = "H. sapiens",species2 = "D. melanogaster",GOterm_field = "term_name",test = "prop")
x_CAT <- x_CAT[which(x_CAT$FDR<=0.05),]
#View Categories with FDR <0.05 (RCD,SPS,GIM, AIM,ERI,DCE)

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

Sosa, Chrystian C., Diana Carolina Clavijo-Buriticá, Victor Hugo García-Merchán, Nicolas López-Rozo, Camila Riccio-Rengifo, Maria Victoria Diaz, David Arango Londoño, y Mauricio Alberto Quimbaya. «GOCompare: An R Package to Compare Functional Enrichment Analysis between Two Species». Genomics 115, n.º 1 (January 2023): 110528. https://doi.org/10.1016/j.ygeno.2022.110528.

## License
GNU GENERAL PUBLIC LICENSE Version 3
