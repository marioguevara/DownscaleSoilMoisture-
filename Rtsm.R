
rm(list=ls())
#kknn tuning parameters
kmax=25
kvc=10
repeats=3
samps=1000
#user provided information
year_of_prediction <- 2018
path_to_target <- '/home/mguevara/Downloads/annual/'
path_to_covariates <- '/home/mguevara/Downloads/ESSDD_covariates_v2.rds'
external_covariates_by_country <- 'RUS'
external_covariates_by_DEM <- 0
covariate_space <- 0
polygon_of_interest <- 0
country_of_interest <-0
eco__climate_area <- 0
get_bioclim_vars = TRUE
predictor_variables <- c(1:44)
# Available variables to predict soil moisture:
# [1] "aspect"                           "carea"                            "chnl_base"                        "chnl_dist"
# [5] "convergence"                      "hcurv"                            "land"                             "lsfactor"
# [9] "rsp"                              "shade"                            "sinks"                            "slope"
# [13] "vall_depth"                       "vcurv"                            "wetness"                          "available.water.storage.capacity"
# [17] "Boreal.coniferous.forest"         "Boreal.mountain.system"           "Boreal.tundra.woodland"           "No.data"
# [21] "Polar"                            "Subtropical.desert"               "Subtropical.dry.forest"           "Subtropical.humid.forest"
# [25] "Subtropical.mountain.system"      "Subtropical.steppe"               "Temperate.continental.forest"     "Temperate.desert"
# [29] "Temperate.mountain.system"        "Temperate.oceanic.forest"         "Temperate.steppe"                 "Tropical.desert"
# [33] "Tropical.dry.forest"              "Tropical.moist.deciduous.forest"  "Tropical.mountain.system"         "Tropical.rainforest"
# [37] "Tropical.shrubland"               "Water"                            "pi0.00"                           "pi0.17"
# [41] "pi0.33"                           "pi0.50"                           "pi0.67"                           "pi0.83"
# [45] "coord.Dim.1"                      "coord.Dim.2"                      "coord.Dim.3"                      "coord.Dim.4"
# [49] "coord.Dim.5"                      "cos2.Dim.1"                       "cos2.Dim.2"                       "cos2.Dim.3"
# [53] "cos2.Dim.4"                       "cos2.Dim.5"                       "contrib.Dim.1"                    "contrib.Dim.2"
# [57] "contrib.Dim.3"                    "contrib.Dim.4"                    "contrib.Dim.5"                    "dist"

topo_sm_15km_1991_2018_kknn<- function(year_of_prediction,  path_to_target, path_to_covariates, predictor_variables) {

library(raster)
if(class(external_covariates_by_country)=="character"){
library(raster)
print('downloading_elevation_data_from_SRTM')
elevation <- getData('alt', country=external_covariates_by_country)
}

if(class(external_covariates_by_DEM)=="RasterLayer"){
library(raster)
elevation <-external_covariates_by_DEM
}

print('verifying_data_structure')
if(class(elevation)=="list"){
print(paste0('this_object_list_of_rasters_including_', length(elevation), '_using_', length(elevation[[1]]), 'kb_for_',
external_covariates_by_country))
print(paste0('using_only_the_first_element'))
print(paste0('aggregating_spatial_resolution_by_a_mean_factor_of_5_on_this_big_file'))
print(paste0('download_the_other_elements_using_raster::getData'))
print(paste0('provide_them_independenlty_as_external_covariates_by_DEM'))
elevation <- elevation[[1]]
elevation <- aggregate(elevation, 5, mean)
print(class(elevation))
}

print('calculating_slope_aspect_and_TPI_from_elevation')
slp_asp <- terrain(elevation, opt=c('slope', 'aspect'), unit='degrees')
shade <- hillShade(slp_asp[[1]], slp_asp[[2]], 40, 270)
names(shade) <- 'shade'
# TPI for different neighborhood size:
tpiw <- function(x, w=5) {
	m <- matrix(1/(w^2-1), nc=w, nr=w)
	m[ceiling(0.5 * length(m))] <- 0
	f <- focal(x, m)
	x - f
}
tpi <- tpiw(elevation, w=5)
names(tpi) <- 'tpi'
covariate_space <- as(stack(elevation, slp_asp, tpi, shade), 'SpatialPixelsDataFrame')
NA2median <- function(x) replace(x, is.na(x), median(x, na.rm = TRUE))
covariate_space@data[] <- lapply(covariate_space@data, NA2median)
predictor_variables <- names(covariate_space)


if(get_bioclim_vars==TRUE){
print('downloading_bioclimatic_variables_from_worldclim')
bio <- raster::getData('worldclim', var='bio',  res=10)
bio <- crop(bio, elevation)
bio <- projectRaster(bio, elevation)
bio <- mask(bio, elevation)
print('reducing_dimensions_of_bioclimatic_data')
library(RStoolbox)
bio <- rasterPCA(bio)
 print('triming_covariates_to_study_area')
covariate_space  <- trim(stack(stack(covariate_space), stack(bio$map[[1:3]])))
#library(spatialEco)
#x <- sp.na.omit(x)
covariate_space  <- as(covariate_space , 'SpatialPixelsDataFrame')
NA2median <- function(x) replace(x, is.na(x), median(x, na.rm = TRUE))
covariate_space@data[] <- lapply(covariate_space@data, NA2median)
predictor_variables <- names(covariate_space)
}

lis <- list.files(path_to_target, full.names=TRUE)
library(raster)
library(kknn)
library(spatialEco)
library(rgdal)
#load covariates

x1 <- covariate_space

if(class(x1)=="numeric"){
x <- readRDS(path_to_covariates)
x1 <- x[predictor_variables]
}

if(is.character(eco__climate_area)==TRUE){
x1 <- x[grep(eco__climate_area, names(x))]
x1 <- raster(x1)
x1[x1==0] <- NA	
x1 <- x[as(x1, 'SpatialPixelsDataFrame'),]
 }
 
if(is.character(country_of_interest)==TRUE){
 print('downloading_country_limit')
 country_of_interest <- getData('GADM', country=country_of_interest, level=0)

 x1 <- x1[predictor_variables]
 print('triming_covariates_to_study_area')
 x1 <- as(trim(stack(x[country_of_interest,])), 'SpatialPixelsDataFrame')
 }
 
if(is.numeric(polygon_of_interest)==FALSE){
 
 x1 <- x1[predictor_variables]
 print('triming_covariates_to_study_area')
 x1 <- as(trim(stack(x[polygon_of_interest,])), 'SpatialPixelsDataFrame')
 }

training_image <- raster(lis[grep(year_of_prediction, lis)])
#define proyection
#proj4string(training_image) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

training_image <-projectRaster(training_image , raster(x1), method='ngb')

#convert to data frame
df=raster::as.data.frame(training_image, xy=T)
#remove NAs
df=na.omit(df)
#define coordinates
coordinates(df)=~x+y
#define proyection
proj4string(df) <- CRS(projection(training_image))
df <- spTransform(df, CRS(projection(x1)))
#proj4string(x1) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
ov=over(df, x1)
#to data frame
d=as.data.frame(df)
#combine extraction with data
X=cbind(d,  ov)
#copy data frame (only columns of interes)
df=data.frame(y=X[,1], X[,4:dim(X)[2]])
#only oblique coordinates
#remove NAs again
df=na.omit(df)
#df=data.frame(y=X[,3], X[,4:9])
#regression matrix
print('n_train , n_cols')
print(dim(df))

preds <- stack()
results1 <- data.frame()

for (j in 1:repeats){	

df1 <- df[sample(nrow(df), samps), ]

##you must select the best parameters by tunning them with CV,  the parameter K and the parameter kernel
#tune the model
#df <- df[1:100,]
 knnTuning <- train.kknn(y~., data=df1, kmax = kmax, distance = 2,
                         kernel = c("rectangular", "triangular", "epanechnikov","gaussian", "rank", "optimal"),
                         ykernel = NULL, scale = TRUE,kcv=kvc)
# #extract best paramters
 n<-which(knnTuning$best.parameters$kernel==c("rectangular", "triangular", "epanechnikov","gaussian", "rank", "optimal"))
# #save best paramters
 mejoresresultados <- data.matrix(unlist(knnTuning$fitted.values[[(kmax*(n-1))+knnTuning$best.parameters$k]]))
# #calculate RMSE
 rmse <- sqrt(knnTuning$MEAN.SQU[knnTuning$best.parameters$k,n])
# #calculate correlation obs pred
 (cd <- cor(df1[,1], mejoresresultados))
# #run the best model (with the best parameters) and make predictions (maps)
# #x <- x[vars][1:100,]

 mejorKKNN <- kknn(y~.,train=df1,test=x1,kernel=unlist(knnTuning$best.parameters[1]),
                   scale=TRUE,k=as.numeric(knnTuning$best.parameters[2]))
# #get the fitted values for all the area (x)
 x1$kknn=mejorKKNN$fitted.values
sm_raw_values <- raster(x1['kknn'])
sm_scaled_values <- raster.transformation(sm_raw_values, trans = "stretch", smin = min(df$y), smax = max(df$y))
 results<-data.frame(year=numeric(), cor=numeric(),rmse=numeric(), n=numeric(),kernel=as.character(),stringsAsFactors=FALSE,k=numeric())
 k=as.numeric(knnTuning$best.parameters[2])
 #store in the results data frame (month), replace by year if using yearly averages
 results[1,1]<-year_of_prediction
 #correlation obs pred
 results[1,2]<-cd
#root mean squared error
 results[1,3]<-rmse
 #data available for that month/year
 results[1,4]<-dim(df1)[1]
 #best kernel
 results[1,5]<-unlist(knnTuning$best.parameters[1])
#k value
results[1,6]<-k
#print results for model i
results1 <- rbind(results1, results)
preds <- stack(preds, sm_scaled_values)
print(paste0('repeat_' , j, '_of_', repeats, '_realizations'))

}
print(paste0('calculating_mean_value_of_', repeats, '_realizations'))
sm_predicted <- calc(preds, mean)
print(paste0('calculating_sdev_value_of_', repeats, '_realizations'))
sm_predicted_sd <- calc(preds, sd)

print(year_of_prediction)
print(cd)
print(rmse)
training_image <- crop(training_image, sm_predicted)
training_image <- mask(projectRaster(training_image, sm_predicted, method = 'ngb'), sm_scaled_values)
s <- stack(sm_predicted, sm_predicted_sd, training_image)
names(s) <- c('predicted_soil_moisture', 'prediction_variance', 'original_satellite_soil_moisture')
output <-list(s, results1, x1)
return(output)
}

out <- topo_sm_15km_1991_2018_kknn(year_of_prediction,  path_to_target, path_to_covariates, predictor_variables)
plot(raster(out[[3]]['shade']), col=grey(0:100/100), main='', cex.axis=1.5, legend=FALSE)
plot(out[[1]][[1]], col=rainbow(7, alpha=0.35), main='', cex.axis=1.5, add=TRUE)
#library(mapview)
#mapview::mapview(out[[1]][[1]])
#mapa <- out[[1]][[1]]
#library(plotKML)
#plotKML(mapa,  colour_scale = rev(SAGA_pal[[1]]))


