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

#start 2pm anual 2017 2018
#libraries
library(rgdal)
library(raster)
library(kknn)
#list to paths and names of ESA raster files
(lis <- list.files('/home/mguevara/Downloads/ESACCI-SM-COMBINED_2017-2018-v045_metrics/', full.names=TRUE))
#e.g., select only 1 year in a monthly basis
#here selecting the mean of 2017 and 2018
lis <- lis[c(27, 40)]
#load covariates
x <- readRDS('covariates_dem_topo_awc_hwsd_eco_type_oblique_xy_PCA.rds')
#empty data frame for storing results
results<-data.frame(year=numeric(), cor=numeric(),rmse=numeric(), n=numeric(),kernel=as.character(),stringsAsFactors=FALSE,k=numeric())
#initialize the loop
for (i in 1:length(lis)){
#ras <- raster('/home/mguevara/Downloads/ESACCI-SM-COMBINED_2017-2018-v045_metrics/World_sm_mean_2017_version_45.tif')
#read raster i
ras <- raster(lis[i])
#define proyection
#proj4string(ras) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
#convert to data frame
df=as.data.frame(ras, xy=T)
#remove NAs
df=na.omit(df)
#define coordinates
coordinates(df)=~x+y
#define proyection
proj4string(df) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
#extract values to points, only the first 5 PCAs
#change here if needed (names(df))
ov=over(df, x[45:49])
#to data frame
d=as.data.frame(df)
#combine extraction with data
X=cbind(d,  ov)
#remove NAs again
X=na.omit(X)
#copy data frame (only columns of interes)
df=data.frame(y=X[,3], X[,4:8])
#only oblique coordinates
#df=data.frame(y=X[,3], X[,4:9])
names(df)
print(dim(df))
##you must select the best parameters by tunning them with CV,  the parameter K and the parameter kernel
#test a max value of 25 for k
kmax=25
#tune the model
knnTuning <- train.kknn(y~., data=df, kmax = kmax, distance = 2,
kernel = c("rectangular", "triangular", "epanechnikov","gaussian", "rank", "optimal"),
           ykernel = NULL, scale = TRUE,kcv=10)
#extract best paramters
n<-which(knnTuning$best.parameters$kernel==c("rectangular", "triangular", "epanechnikov","gaussian", "rank", "optimal"))
#save best paramters
mejoresresultados <- data.matrix(unlist(knnTuning$fitted.values[[(kmax*(n-1))+knnTuning$best.parameters$k]]))
#calculate RMSE
rmse <- sqrt(knnTuning$MEAN.SQU[knnTuning$best.parameters$k,n])
#calculate correlation obs pred
(cd <- cor(df[,1], mejoresresultados))
#run the best model (with the best parameters) and make predictions (maps)
mejorKKNN <- kknn(y~.,train=df,test=x,kernel=unlist(knnTuning$best.parameters[1]),
                        scale=TRUE,k=as.numeric(knnTuning$best.parameters[2]))
#get the fitted values for all the area (x)
x$kknn=mejorKKNN$fitted.values
#coordinates(kknn)=~x+y
#gridded(kknn) = TRUE
#kknn=as(kknn,'SpatialPixelsDataFrame')
r <- raster(x['kknn'])
#save the prediction raster (check year/month )
writeRaster(r, file=paste0('predicted_month_', 2016+i, '_SM_15km_essd_ESA_CCI_env_corr_V2.tif'), overwrite=TRUE) ### 
#save month, replace by year if using yearly averages
year=2016+i
#best k value
k=as.numeric(knnTuning$best.parameters[2])
#store in the results data frame (month), replace by year if using yearly averages
results[i,1]<-year
#correlation obs pred
results[i,2]<-cd
#root mean squared error
results[i,3]<-rmse
#data available for that month/year
results[i,4]<-dim(df)[1]
#best kernel
results[i,5]<-unlist(knnTuning$best.parameters[1])
#k value
results[i,6]<-k
#print results for model i
print(i)
print(cd)
print(rmse)
	}
#modeling results 
print(results)
#save csv with modeling results
write.csv(results, file='ESSD_sm_prediction_kknn_2017_2018_yearly.csv' )
#end



#tiles system
e <- extent(s)
e <- as(e, 'SpatialPolygons')
tiles <- getSpatialTiles(e, 10, 10,  overlap.percent = 0.5)
proj4string(tiles) <- CRS(projection(ras))
tiles <- tiles[eco_zones,]
#build a regression matrix
ras <- cbind(ras@data, ras@coords, extract(s, ras))
ras$eco_zone <- as.factor(ras$eco_zone)
lev <- levels(ras$eco_zone)
#Tropical rain forest is 19
#registerDoMC(6)
#for(i in 1: length(lev)){
#train <- ras[ras$eco_zone == lev[i],]
train <- ras
train <- na.omit(train)
train <- train[train$World_sm_mean_2017_version_45>0,]
y <- train[,1]
x <- train[5:19]
df <- data.frame(y=y, x)


kmax=30
#find op	timal k and kernel type via 10 fold cross validation
knnTuning <- train.kknn(y~., data=df, kmax = kmax, distance = 2,
kernel = c("rectangular", "triangular", "epanechnikov","gaussian", "rank", "optimal"),
           ykernel = NULL, scale = TRUE, kcv=10)

#identify best kernel
n<-which(knnTuning$best.parameters$kernel==c("rectangular", "triangular", "epanechnikov","gaussian", "rank", "optimal"))
#store best parameters
mejoresresultados <- data.matrix(unlist(knnTuning$fitted.values[[(kmax*(n-1))+knnTuning$best.parameters$k]]))
#lowest rmse
rmse <- sqrt(knnTuning$MEAN.SQU[knnTuning$best.parameters$k,n])
#highest correlation
cd <- cor(y, mejoresresultados)
#best results

covs <- crop(s, eco_zones[eco_zones$GEZ_TERM==lev[19],])

e <- crop(covs, drawExtent())

e <- as(e, 'SpatialPixelsDataFrame')

mejorKKNN <- kknn(y~.,train=df, test=e@data,kernel=unlist(knnTuning$best.parameters[1]),
                        scale=TRUE,k=as.numeric(knnTuning$best.parameters[2]))

x$kknn=mejorKKNN$fitted.values
#make a raster of the prediction
r <- raster(x['kknn'])

#modeling year (1991-2016)
year=1990+i
#prepare best kernel
as.character(unlist(knnTuning$best.parameters[1]))
#prepare best k
k=as.numeric(knnTuning$best.parameters[2])
#store results for year i in the data frame 
results[i,1]<-year
results[i,2]<-cd
results[i,3]<-rmse
results[i,4]<-dim(z)[1]
results[i,5]<-unlist(knnTuning$best.parameters[1])
results[i,6]<-k
#print accuracy results on screen 
print(results)


# define the control using a random forest selection function
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# run the RFE algorithm
registerDoMC(6)
results <- rfe(x, y, sizes=c(1:8), rfeControl=control)
# summarize the results
#14:37

y <- train[,1]
x <- train[5:19]






e <- extent(s)
e <- as(e, 'SpatialPolygons')
tiles <- getSpatialTiles(e, 10, 10,  overlap.percent = 0.5)
tiles <- tiles[eco_zones,]

