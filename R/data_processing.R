
## PREPARE DATA FOR GLOBAL MODELS
## Note: For gdal to work on Mac in R, see: https://github.com/smwindecker/gdaltools


## Set work environmet
## ------------------------------------
library(sp)
library(raster)
library(rgdal)
devtools::install_github('skiptoniam/sense')
library(sense) 
library(tools)
library(data.table)



## Biodiversity data
## ------------------------------------
## Occurrence records for Tyto alba
## GBIF.org (17 September 2018) GBIF Occurrence Download https://doi.org/10.15468/dl.7kootk
raw_dat <- fread("data/gbif_tytoalba/0012585-180824113759888.csv", quote = "")


## GBIF backbone taxonomy
## GBIF Secretariat (2017). GBIF Backbone Taxonomy. Checklist dataset https://doi.org/10.15468/39omei accessed via GBIF.org on 2019-08-26.
backbone <-fread("data/Taxon.tsv")


## Clean GBIF data
source("R/filter.gbifraw.R")
dat <- filter.gbifraw(raw_dat, backbone, subset.gbifnubtaxonomy.byclass = "Aves",
                      output_folder = "data/processed/", domain.mask = global_mask)
  ## view log file
  ## explore dat

## Remove unwanted columns from data
dat <- dat[,.(species, decimallatitude, decimallongitude)]

rm(raw_dat, backbone)
gc()


## Covariate data
## ------------------------------------
## WorldClim: http://www.worldclim.org/ (10 min resolution)
bio_current <- list.files("data/wc10", pattern = '.bil$', full.names = TRUE)
bio_future <- list.files("data/cmip5/10m", full.names = TRUE)
bio_rcp26 <- bio_future[grepl(paste0("(?=.*26)(?=.*tif)"), bio_future, perl = TRUE)]
bio_rcp85 <- bio_future[grepl(paste0("(?=.*85)(?=.*tif)"), bio_future, perl = TRUE)]
bioclim <- c(bio_current, bio_rcp26, bio_rcp85)


## Elevation: https://webmap.ornl.gov/ogc/dataset.jsp?dg_id=10008_1
srtm <- "data/srtm.tif"


## Soil: https://daac.ornl.gov/SOILS/guides/igbp-surfaces.html
  ## bulkdens - bulk density (weight of soil in a given volume)
  ## pawc - profile available water capacity (amount of water that can be stored in a soil profile)
  ## soilcarb - soil carbon density 
  ## totaln - total nitrogen density
soil <- list.files("data/soil", pattern = "*.dat", full.names = T)


## Landuse: https://earthexplorer.usgs.gov/
landuse <- "data/landuse.tif"


## Create Mask from WorldClim layer (10 minn resolution)
global_mask <- raster("data/wc10/bio1.bil")
global_mask[which(!is.na(global_mask[]))] <- 1


## Define extent
e <- c(-180,180,-60,90)


## Define resolution
reso <- res(global_mask)


## Specify files for sense::gdal_crop function
file_in <- c(bioclim, srtm, soil, landuse)
file_out <- paste0("data/processed/", paste0(tools::file_path_sans_ext(basename(file_in)), "_treated.", tools::file_ext(file_in)))


## Crop layers as per mask
## see: https://github.com/skiptoniam/sense/blob/master/R/gdal_raster_functions.R
mapply(gdal_crop, inpath = file_in, outpath = file_out, MoreArgs = list(extent=e,res=reso)) 


## SRTM variables
elevation <- raster("data/processed/srtm_treated.tif")
aspect <- terrain(elevation, opt = "aspect")
slope <- terrain(elevation, opt = 'slope')
roughness <- terrain(elevation, opt = "roughness")


## Landuse reclassification
## Original land use classes in GLCC data: https://www.usgs.gov/media/images/global-land-cover-characteristics-data-base-version-20
##    # 0	Water
##    # 1	Evergreen Needle leaf Forest
##    # 2	Evergreen Broadleaf Forest
##    # 3	Deciduous Needle leaf Forest
##    # 4	Deciduous Broadleaf Forest
##    # 5	Mixed Forests
##    # 6	Closed Shrublands
##    # 7	Open Shrublands
##    # 8	Woody Savannas
##    # 9	Savannas
##    # 10	Grasslands
##    # 11	Permanent Wetland
##    # 12	Croplands
##    # 13	Urban and Built-Up
##    # 14	Cropland/Natural Vegetation Mosaic
##    # 15	Snow and Ice
##    # 16	Barren or Sparsely Vegetated
landuse <- raster("data/processed/landuse_treated.tif")
landuse[landuse[]%in%c(12,14)] <- 100
landuse[landuse[]%in%c(6,8,7,9,10)] <- 200
landuse[landuse[]%in%c(1,2,3,4,5)] <- 300
landuse[landuse[]%in%c(13)] <- 400
landuse[landuse[]%in%c(0, 11,15,16,17)] <- 500
landuse <- landuse/100 - 1 ## these steps ensure that substituted alues do not overlap wth existing values
  ## check
  unique(values(landuse))

  
## Stack, mask and save covariates as .rds
covariates_all <- stack(setdiff(file_out, file_out[grep("srtm|landuse", file_out)]), 
                        elevation, slope, roughness, aspect, landuse) 
  ## load and stack covariates files, except files for srtm and landuse
  ##  which were updated and therefore these rasters are loaded from memory.
crs(covariates_all) <- crs(global_mask)
# covariates_all <- mask(covariates_all, global_mask)
names(covariates_all) <- sub("_treated","", names(covariates_all))
saveRDS(covariates_all, file = "data/processed/covariates_all.rds")
  
  ## Check to see same number of NAS in covariates as in global mask
  summary(global_mask)
  summary(covariates_all)

  
## Check for correlations and subset covariates: for modelling 
covariates_all <- readRDS(paste0(data_gsdms, "/covariates_all.rds"))
cov_keep <- c("bio1", "bio4","bio12", "bio15", "bulkdens","pawc","soilcarb","totaln",
                          "srtm","slope","roughness","aspect","landuse")
covariates <- covariates_all[[which(names(covariates_all) %in% cov_keep)]]
cov_values <- getValues(covariates)
corr_val <- cor(cov_values, use = 'complete.obs', method = 'pearson') 
rm(cov_values)


## Visulaisation
  library(corrplot)
  library(ggplot2)
  library(mnormt); library(psych)
  library(reshape); library(GGally)
corrplot::corrplot(corr_val, type = "upper", method = "number")
psych::pairs.panels(corr_val, scale = TRUE)


## Remove highly correlated covariates (> 0.8)
'%!in%' <- function(x,y)!('%in%'(x,y))
covariates <- covariates[[which(names(covariates) %!in% c("bio4", "totaln", "roughness"))]]  
saveRDS(covariates, file = "data/processed/covariates.rds")


## Check for correlations and subset covariates: for prediction
rm(cov_keep)
cov_keep <- names(covariates_all)[grep('26|85', names(covariates_all))]
cov_keep <- cov_keep[grep('701$|704$|7012$|7015$', cov_keep)]
cov_keep <- c(names(covariates), cov_keep)
covariates_predict <- covariates_all[[which(names(covariates_all) %in% cov_keep)]]
covariates_predict <- covariates_predict[[which(names(covariates_predict) %!in% c("bc26bi704" ,"bc85bi704","totaln", "roughness"))]]
saveRDS(covariates_predict, file = "data/processed/covariates_predict.rds")

