#function to make dummies from raster factor maps
dummyRaster <- function(rast){
  rast <- as.factor(rast)
  result <- list()
  for(i in 1:length(levels(rast)[[1]][[1]])){
    result[[i]] <- rast == levels(rast)[[1]][[1]][i]
    names(result[[i]]) <- paste0(names(rast), 
                                 levels(rast)[[1]][[1]][i])
  }
  return(stack(result))
}
#libraries
library(GSIF)
library(raster)
library(OGC)
library(rgeos)
library(sp)
library(caret)
library(FactoMineR)
library(kknn)

#topographic parameters
#s <- stack(list.files(pattern='sdat'))
#s <- aggregate(s, 15, mean)
#NA2mean <- function(x) replace(x, is.na(x), median(x, na.rm = TRUE))
#x=as(x,'SpatialPixelsDataFrame')
#x=as(s,'SpatialPixelsDataFrame')
#x@data[] <- lapply(x@data, NA2mean)
#oblique coordinates and topo
#x<- as(stack(S, projectRaster(stack(ogcs), S[[1]], method='ngb')), 'SpatialPixelsDataFrame')
#only oblique coordinates(n 6)
#x <- ogcs

#http://193.43.36.20/map?entryId=baa463d0-88fd-11da-a88f-000d939bc5d8
#ecological zones map FAO
eco_zones <- shapefile('/home/mguevara/Downloads/eco_zone/eco_zone.shp')
#define projection
proj4string(eco_zones) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#factorize the column with the ecological areas
eco_zones@data$GEZ_TERM  <- as.factor(eco_zones@data$GEZ_TERM)
#read topographic predictors 15km grids
x <- readRDS('/home/mguevara/Downloads/topographic_predictors_15km_grids.rds')
#copy them in a raster stack
S <-stack(x)
#rasterize ecological map 
eco.r <- rasterize(x = eco_zones, y = S[[1]], field = "GEZ_TERM")
#create dummy variables
eco.r_dummy <- dummyRaster(eco.r)
#get the names
names(eco.r_dummy) <- levels(eco_zones$GEZ_TERM)
#to spatial pixels
eco.r_dummy <- as(eco.r_dummy, 'SpatialPixelsDataFrame')
#ogcs <- makeOGC(eco.r_dummy[[1]], 6)
#create oblique coordinates
ogcs <- makeOGC(eco.r, 6)
#transform to spatial pixels
ogcs <- as(ogcs, 'SpatialPixelsDataFrame')
#combine oblique coordinates and ecological dummies
eco_coords <- cbind(eco.r_dummy, ogcs)
#load available water capacity from the regridded HWSD
#https://daac.ornl.gov/SOILS/guides/HWSD.html
awc <- raster('/home/mguevara/Downloads/AWC_CLASS.nc4')
#use a 15 km grid integrity
awc <- projectRaster(awc, S[[1]], method='ngb')
#combine topography, awc, ecological zones, oblique coordinates 
x<- as(stack(S, awc, projectRaster(stack(eco_coords), S[[1]], method='ngb')), 'SpatialPixelsDataFrame')
#fit a PCA to reduce dimensions
fit <- PCA(x@data,  graph=FALSE)
#get the PCAs
x@data <- cbind(x@data, data.frame(fit$ind))
#save prediction domain (covariate file)
saveRDS(x, file='covariates_dem_topo_awc_hwsd_eco_type_oblique_xy_PCA.rds')
#end
