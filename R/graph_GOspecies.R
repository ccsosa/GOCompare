#' @title Comprehensive visual and graph comparison between two species and a series of categories
#' @name graphGOspecies
#' @description The GRSex process provides a geographic measurement of the proportion of a species’ range
#'  that can be considered to be conserved in ex situ repositories. The GRSex uses buffers (default 50 km radius)
#'  created around each G coordinate point to estimate geographic areas already well collected within the distribution
#'  models of each taxon, and then calculates the proportion of the distribution model covered by these buffers.
#' @param Occurrence_data A data frame object with the species name, geographical coordinates,
#'  and type of records (G or H) for a given species
#' @param Species_list A vector of characters with the species names to calculate the GRSex metrics.
#' @param Raster_list A list of rasters representing the species distribution models for the species list provided
#'  in \var{Species_list}. The order of rasters in this list must match the same order as \var{Species_list}.
#' @param Buffer_distance Geographical distance used to create circular buffers around germplasm.
#'  Default: 50000 (50 km) around germplasm accessions (CA50)
#' @param Gap_Map logical, if \code{TRUE} the function will calculate gap maps for each species analyzed and
#'  will return a list with two slots GRSex and gap_maps. If any value is provided, the function will assume that
#'  Gap_Map = TRUE
#' @return This function returns a data frame with two columns:
#'
#' \tabular{lcc}{
#' species \tab Species name \cr
#' GRSex \tab GRSex value calculated\cr
#' }
#'
#' @examples
#' ##Obtaining occurrences from example
#' data(CucurbitaData)
#' Cucurbita_splist <- unique(CucurbitaData$species)
#' ## Obtaining rasterList object. ##
#' data(CucurbitaRasters)
#' CucurbitaRasters <- raster::unstack(CucurbitaRasters)
#' #Running GRSex
#' GRSex_df <- GRSex(Species_list = Cucurbita_splist,
#'                     Occurrence_data = CucurbitaData,
#'                     Raster_list = CucurbitaRasters,
#'                     Buffer_distance = 50000,
#'                     Gap_Map = TRUE)
#'
#' @references
#' Ramirez-Villegas et al. (2010) PLOS ONE, 5(10), e13497. doi: 10.1371/journal.pone.0013497
#' Khoury et al. (2019) Ecological Indicators 98:420-429. doi: 10.1016/j.ecolind.2018.11.016
#'
#' @export
#' @importFrom sp coordinates proj4string SpatialPoints over CRS
#' @importFrom stats median
#' @importFrom fasterize fasterize
#' @importFrom raster overlay crop raster extent ncell projection


dir2 <- "D:/PROGRAMAS/Dropbox/shared/CANCER_HALLMARKS/GOCompare/GOCompare/data"

load(paste0(dir2,"/","A_thaliana.rda"))
load(paste0(dir2,"/","H_sapiens.rda"))


df1 <- H_sapiens
df2 <- A_thaliana
GOterm_field <- "Functional.Category"
species1 <- "H. sapiens"
species2 <- "A. thaliana"


compareGOspecies <- function(df1,df2,GOterm_field,species1,species2){

  comb1 <- paste(species1,"-",unique(df1[,"feature"]))
  comb1_feat <- unique(df1[,"feature"])
  comb2 <- paste(species2,"-",unique(df2[,"feature"]))
  comb2_feat <- unique(df2[,"feature"])
  comb_all_feat <- c(comb1,comb2)
  #comb_feat_3 <- unique(comb1_feat,comb2_feat)

  GO_terms_unique <- unique(df1[,GOterm_field],df2[,GOterm_field])
  mat_for_dist <- as.data.frame(matrix(nrow=length(comb1)+length(comb2),ncol = length(GO_terms_unique)))
  row.names(mat_for_dist) <- comb_all_feat
  #colnames(mat_for_dist) <- GO_terms_unique
  comb_all_feat <- strsplit(comb_all_feat," - ")
  comb_all_feat <- lapply(1:length(comb_all_feat),function(i){
    x <- data.frame(species=comb_all_feat[[i]][1],feature=comb_all_feat[[i]][2])
  })

  comb_all_feat <- do.call(rbind,comb_all_feat)
  row.names(comb_all_feat)
  df1$species <-species1
  df2$species <-species2
  df_total <- rbind(df1,df2)


  pb <- txtProgressBar(min = 0, max = nrow(mat_for_dist), style = 3)

  for(i in 1:nrow(mat_for_dist)){
    #message(paste(round((i/nrow(mat_for_dist)*100),2),"%"))
    setTxtProgressBar(pb, i)

    for(j in 1:ncol(mat_for_dist)){
      if(nrow(df_total[which(df_total$feature == comb_all_feat[i,2] &
                             df_total$species == comb_all_feat[i,1] &
                             df_total[,GOterm_field]==GO_terms_unique[[j]]),])>0){
        mat_for_dist[i,j] <- 1
      } else {
        mat_for_dist[i,j] <- 0
      }

    };rm(j)
  };rm(i)
  close(pb)

  jacc_dist <- vegan::vegdist(mat_for_dist,method = "jaccard",na.rm = T)

  vare.mds <- ape::pcoa(jacc_dist)
  vare.mds2 <- as.data.frame(vare.mds$vectors[,1:2])
  vare.mds2$species <- unlist(lapply(strsplit(row.names(vare.mds2),"-"), `[[`, 1))
  vare.mds2$species <- trimws(vare.mds2$species)

  grp.a <- vare.mds2[vare.mds2$species == species1, ][chull(vare.mds2[vare.mds2$species ==
                                                                        species1, c("Axis.1", "Axis.2")]), ]  # hull values for grp A
  grp.b <- vare.mds2[vare.mds2$species == species2, ][chull(vare.mds2[vare.mds2$species ==
                                                                        species2, c("Axis.1", "Axis.2")]), ]  # hull values for grp A
  hull.data <- rbind(grp.a, grp.b)  #combine grp.a and grp.b

  species.scores <- as.data.frame(vare.mds2)  #Using the scores function from vegan to extract the species scores and convert to a data.frame
  species.scores$grp <- unlist(lapply(strsplit(row.names(species.scores),"-"), `[[`, 2))

  ov_plot <-
    ggplot() +
    geom_polygon(data=hull.data,aes(x=Axis.1,y=Axis.2,fill=species,group=species),alpha=0.30) + # add the convex hulls
    #scale_fill_manual(values=c(species1 = "green", species2= "blue")) +
    geom_point(data=vare.mds2,aes(x=Axis.1,y=Axis.2,colour=species),size=8,alpha=0.5) + # add the point markers
    geom_text_repel(data=species.scores,aes(x=Axis.1,y=Axis.2,label=grp),
                    alpha=0.5,size=5, show.legend = FALSE,colour= 'black',na.rm=T
                    # force = 30, segment.colour = NA
    )

  comb_feat_3 <-unique(comb_all_feat$feature)

  shared_GO_list <- list()
  for(i in 1:length(comb_feat_3)){

    x <- as.numeric(row.names(comb_all_feat[comb_all_feat$feature %in% comb_feat_3[[1]],]))
    x <-mat_for_dist[x,]
    x_col <- data.frame(feature=comb_feat_3[i],GO=GO_terms_unique,status=colSums(x))
    x_col <- x_col[which(x_col$status==2),]
    x_col$status <- NULL
    shared_GO_list[[i]] <- x_col
  };rm(i)

  shared_GO_list <- do.call(rbind,shared_GO_list)
  row.names(shared_GO_list) <- NULL
  ov_plot_list <- list(graphic= ov_plot, distance = jacc_dist,shared_GO_list=shared_GO_list)

  return(ov_plot_list)
}


x <- compareGOspecies(df1,df2,GOterm_field,species1,species2)
