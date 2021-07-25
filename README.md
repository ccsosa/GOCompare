# GOCompare R package

## Description

## Installation
GOCompare can be installed as follows
```r
#CRAN
#install.packages("GOCompare") #TO BE SUBMITTED TO CRAN
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
The recommended workflow is as follows:


## Authors
Main:Chrystian C. Sosa, Diana Carolina Clavijo-Buriticá, Mauricio Quimbaya

Other contributors: Maria Victoria Diaz, Camila Riccio Rengifo, Arlen James Mosquera, Andrés Álvarez

## References

Clavijo-Buriticá, Sosa, C.C., Mosquera, A.J. Álvarez, A., Medina, J. Quimbaya, M.A. A systematic comparison of the molecular machinery associated with Cancer-Hallmarks between plants and humans reveals Arabidopsis thaliana as a useful model to understand specific carcinogenic events (to be submitted, Target journal: Plos Biology)
 
 
## License
GNU GENERAL PUBLIC LICENSE Version 3
