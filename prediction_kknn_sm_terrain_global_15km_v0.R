 qlogin -l exclusive=1
vpkg_require r/3.3.3
#r-search GSIF
#r-search caret
vpkg_require r-cran/20170524
vpkg_require r-gdal/20170524
cd /home/work/spac/mguevara/SOIL_MOISTURE_ESA/

R

rm(list=ls())
###time counter
CronometroON<- function(){
      tic<-proc.time()[3]
      assign(".tic", tic, envir=baseenv())
      invisible(tic)
      }

 CronometroOFF<- function(){
      tic <- get(".tic", envir=baseenv())
      toc<-proc.time()[3]-tic
      hrs<-as.integer(toc/3600)
      minu<- as.integer(((toc/3600)-hrs)*60)
      seg<- ((((((toc/3600)-hrs)*60)))-minu)*60
      time<-paste(as.character(hrs),"hrs ",as.character(minu),"min ",as.character(round(seg,digit=2)),"seg",sep="")
      return(time)
      }
CronometroON()
###generate data bases
#setwd("/home/1579/sample_soil_moisture/input_soil_moisture_stats")### the final outputs will be stored here

library(raster)

#COVARIATES
#x=stack('~/sample_soil_moisture/SAGAoutput.tif')### set the path to the SAGA output

lis1 <- '/home/work/spac/anita/sample/sample_soil_moisture/SAGAIO-361587'
lis2 <-list.files(lis1, full.name=TRUE, pattern ='sdat')
s <- stack(lis2)
#extent to CONUS
#e <- readRDS('/home/work/spac/mguevara/ext1.rds')
#x <- crop (s, e)

s <- aggregate(s, 15, mean)

NA2mean <- function(x) replace(x, is.na(x), median(x, na.rm = TRUE))


#x=as(x,'SpatialPixelsDataFrame')
x=as(s,'SpatialPixelsDataFrame')

x@data[] <- lapply(x@data, NA2mean)


###start a loop
###z will be each year of soil moisture data for the US (from 1978 top 2013)
##save the results of 10 fold cross validation (correlation obs vs mod [cd] and RMSE as well as N) in a data frame named results
library(kknn)
results<-data.frame(year=numeric(), cor=numeric(),rmse=numeric(), n=numeric(),kernel=as.character(),stringsAsFactors=FALSE,k=numeric())

##Start a loop for take all years at once
dirs <- list.dirs()[-1]

i=1
#for (i in 1:26){

r <- stack(paste0(dirs[i], '/',list.files(dirs[i], pattern='nc')), varname='sm')
#r <- crop (r, e)
r <- calc(r, median,  na.rm=TRUE)

proj4string(r) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
df=as.data.frame(r, xy=T)

df=na.omit(df)
coordinates(df)=~x+y
proj4string(df) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

ov=over(df, x)

d=as.data.frame(df)
y=cbind(d,  ov)
y=na.omit(y)

z=data.frame(z=y[,3],y[,1:2],  y[,4:18])
print(dim(z))
##you must select the best parameters by tunning them with CV,  the parameter K and the parameter kernel
##

kmax=50

knnTuning <- train.kknn(z~., data=z, kmax = kmax, distance = 2,
kernel = c("rectangular", "triangular", "epanechnikov","gaussian", "rank", "optimal"),
           ykernel = NULL, scale = TRUE,kcv=10)

n<-which(knnTuning$best.parameters$kernel==c("rectangular", "triangular", "epanechnikov","gaussian", "rank", "optimal"))
mejoresresultados <- data.matrix(unlist(knnTuning$fitted.values[[(kmax*(n-1))+knnTuning$best.parameters$k]]))
rmse <- sqrt(knnTuning$MEAN.SQU[knnTuning$best.parameters$k,n])
cd <- cor(z[,1], mejoresresultados)

year=1990+i
as.character(unlist(knnTuning$best.parameters[1]))
k=as.numeric(knnTuning$best.parameters[2])

results[i,1]<-year
results[i,2]<-cd
results[i,3]<-rmse
results[i,4]<-dim(z)[1]
results[i,5]<-unlist(knnTuning$best.parameters[1])
results[i,6]<-k

##build a model with best parameters, predict to all world and make a map of it
##
mejorKKNN <- kknn(z~.,train=z,test=x,kernel=unlist(knnTuning$best.parameters[1]),
                        scale=TRUE,k=as.numeric(knnTuning$best.parameters[2]))

x$kknn=mejorKKNN$fitted.values

print(CronometroOFF())

#coordinates(kknn)=~x+y
#gridded(kknn) = TRUE
#kknn=as(kknn,'SpatialPixelsDataFrame')
r <- raster(x['kknn'])

#plot(r)
#proj4string(r) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
#writeRaster(r, file=paste0('predicted-', 1990+i, '-esa-sm-topo-CONUS.tif'), overwrite=TRUE) ### change the name of the output files if necessary
writeRaster(r, file=paste0('predicted-', 1990+i, '-esa-sm-topo-GLOBAL.tif'), overwrite=TRUE) ### 
##see cross validation results
##save cross validation results
#(results=results[3:38,])
#write.table(results, file='resultsCONUS.csv', sep=',', col.names=T, row.names=F)
##
print(paste0('predicted-', 1990+i, '-esa-sm-topo-CONUS.tif'))
print(cd)
print(rmse)
}
CronometroOFF()
##END





beginCluster()
     
    r <- stack(list.files(), varname='sm')
     
     f1 <- function(x) calc(x, mean)
     r <- clusterR(r, f1)


proj4string(r) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
df=as.data.frame(r, xy=T)

df=na.omit(df)
coordinates(df)=~x+y
proj4string(df) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

x=stack('/home/1579/sample_soil_moisture/SAGAoutput.tif')### set the path to the SAGA output

