
#time counter function
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

##define prediction coordinates of Delaware
#coords <- c(-76.5, -74.5 , 37.75,  40.25  )
#xmin       : -76.5 
#xmax       : -74.5 
#ymin       : 37.75 
#ymax       : 40.25
#start time 1979
yearmo_init <- 1978
#model results output name
res_out <- 'DE_sm_prediction_kknn_1979_2018_year_basis.csv'
#list to paths and names of ESA raster files
lis <-  list.files('ESA_CCI_SM_v45_Delaware', full.names=TRUE, pattern='tif')
#problem definition (all covariates, monthly predictions year 2017-2018)
probdef <-'delaware_ALL_COVS_1km_'
#libraries
library(rgdal)
library(raster)
library(kknn)
#load covariates
S <- stack(list.files('/home/mguevara/Downloads//SAGAIO-361587/', full.names=TRUE, pattern='sdat'))
refRaster <-extent(raster("/home/mguevara/Downloads//ESA_CCI_SM_v45_Delaware/ESA_CCI_SM_v45_Delaware_20170210.tif"))
S <- crop(S, refRaster)
x <- as(S, 'SpatialPixelsDataFrame' )
NA2mean <- function(x) replace(x, is.na(x), median(x, na.rm = TRUE))
x@data[] <- lapply(x@data, NA2mean)
library(OGC)
ogcs <- makeOGC(raster(x), 6)
x <- cbind(x, as(ogcs, 'SpatialPixelsDataFrame'))


#function to use a median value for NAs, at the borders 
#NA2median <- function(x) replace(x, is.na(x), median(x, na.rm = TRUE))
#replace with median
#x@data[] <- lapply(x@data, NA2median)
#empty data frame for storing results
#all <- stack(lis[-1])

####empty vector for date and year

	date <- character()
	year <- character()


####THESE NUMBERS ARE SENSITIVE TO THE NUMBER OF CHARACTERS IN THE FILES PATH
####SEE EXAMPLE FIRST (function noquote):
#> s <- noquote(strsplit(as.character("/home/mguevara/Downloads//ESA_CCI_SM_v45_Delaware/ESA_CCI_SM_v45_Delaware_20170210.tif"), NULL)[[1]])
#> s
#[1] / h o m e / m g u e v a r a / D o w n l o a d s / / E S A _ C C I _ S M _ v
#[39] 4 5 _ D e l a w a r e / E S A _ C C I _ S M _ v 4 5 _ D e l a w a r e _ 2 0
#[77] 1 7 0 2 1 0 . t i f

	for (i in 1:length(lis[-1])){
	s <- noquote(strsplit(as.character(lis[i+1]), NULL)[[1]])
	date[i] <-paste0((s)[[49]], (s)[[50]], (s)[[51]], (s)[[52]], '-', (s)[[53]], (s)[[54]],'-', (s)[[55]], (s)[[56]])
	year[i]  <- paste0((s)[[49]], (s)[[50]], (s)[[51]], (s)[[52]])
	}

library(package=lubridate)
#Set Weeks number. Date already of class `Date`
Week <- week( as.Date(date))



results1<-data.frame(yearmo=numeric(), cor=numeric(),rmse=numeric(), n=numeric(),kernel=as.character(),stringsAsFactors=FALSE,k=numeric())
predictions1 <- stack()
SDS <- stack()
lev <- levels(as.factor(year))
f1 <- function(x) calc(x, mean, na.rm=TRUE)
f2 <- function(x) calc(x, sd, na.rm=TRUE)

beginCluster()
#i=2018
for (i in 2:length(lev)){
#for (i in 2:3){

#res <- calc(stack(lis[grep(lev[1], lis)]), mean)
ras <- clusterR(stack(lis[grep(lev[1], lis)]), f1)
ras <-projectRaster(ras, raster(x), method='ngb')
df=as.data.frame(ras, xy=T)
#define proyection
#proj4string(ras) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
#convert to data frame
#remove NAs
df=na.omit(df)
#define coordinates
coordinates(df)=~x+y
#define proyection
proj4string(df) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
#extract values to points, only the first 5 PCAs
#change here if needed (names(df))
#ov=over(df, x[c(1:15, 39:44)])
ov=over(df, x)
#to data frame
d=as.data.frame(df)
#combine extraction with data
X=cbind(d,  ov)
#copy data frame (only columns of interes)
df=data.frame(y=X[,3], X[,4:dim(X)[2]])
#only oblique coordinates
#remove NAs again
df=na.omit(df)
#df=data.frame(y=X[,3], X[,4:9])
names(df)
print(dim(df))
predictions <- stack()

results<-data.frame(yearmo=numeric(), cor=numeric(),rmse=numeric(), n=numeric(),kernel=as.character(),stringsAsFactors=FALSE,k=numeric())

for (j in 1:10){	



df1 <- df[sample(nrow(df), 1000), ]

##you must select the best parameters by tunning them with CV,  the parameter K and the parameter kernel
#test a max value of 25 for k
kmax=25
#tune the model
#df <- df[1:100,]
knnTuning <- train.kknn(y~., data=df1, kmax = kmax, distance = 2,
kernel = c("rectangular", "triangular", "epanechnikov","gaussian", "rank", "optimal"),
           ykernel = NULL, scale = TRUE,kcv=10)
#extract best paramters
n<-which(knnTuning$best.parameters$kernel==c("rectangular", "triangular", "epanechnikov","gaussian", "rank", "optimal"))
#save best paramters
mejoresresultados <- data.matrix(unlist(knnTuning$fitted.values[[(kmax*(n-1))+knnTuning$best.parameters$k]]))
#calculate RMSE
rmse <- sqrt(knnTuning$MEAN.SQU[knnTuning$best.parameters$k,n])
#calculate correlation obs pred
(cd <- cor(df1[,1], mejoresresultados))
#run the best model (with the best parameters) and make predictions (maps)
print(CronometroON())
#x <- x[1:100,]
mejorKKNN <- kknn(y~.,train=df1,test=x,kernel=unlist(knnTuning$best.parameters[1]),
                        scale=TRUE,k=as.numeric(knnTuning$best.parameters[2]))
print(CronometroOFF())
#get the fitted values for all the area (x)
x$kknn=mejorKKNN$fitted.values
#coordinates(kknn)=~x+y
#gridded(kknn) = TRUE
#kknn=as(kknn,'SpatialPixelsDataFrame')
r <- raster(x['kknn'])
#save the prediction raster (check year/month )
#writeRaster(r, file=#index
#out <- paste0(probdef, unlist(strsplit(lis[i], '//'))[2]), overwrite=TRUE) ### 
#save month, replace by year if using yearly averages
yearmo=yearmo_init+i
#best k value
k=as.numeric(knnTuning$best.parameters[2])
#store in the results data frame (month), replace by year if using yearly averages
results[j,1]<-yearmo
#correlation obs pred
results[j,2]<-cd
#root mean squared error
results[j,3]<-rmse
#data available for that month/year
results[j,4]<-dim(df)[1]
#best kernel
results[j,5]<-unlist(knnTuning$best.parameters[1])
#k value
results[j,6]<-k
#print results for model i
print(i)
print(cd)
print(rmse)
predictions <- stack(predictions, r)
}
results1 <- rbind(results1, results)

preds <- clusterR(predictions, f1)
sds <- clusterR(predictions, f2)
SDS<- stack(SDS, sds)

predictions1 <- stack(predictions1, preds)

}
endCluster()

#modeling results 
print(results)


####
####
####

####
####
####

####
####
####















names(all[[grep(i, names(all))]])



for( j  in 1:12){

S2 <- all[[grep(j, names(S))]]


}
print(i)
}



tm <- seq(as.Date('1978-11-01'), as.Date('2018-12-31'), 'days')
s <- setZ(all, tm, 'days')
library(zoo)
x <- zApply(all, by=as.yearqtr, fun=mean, name='quarters')




results<-data.frame(yearmo=numeric(), cor=numeric(),rmse=numeric(), n=numeric(),kernel=as.character(),stringsAsFactors=FALSE,k=numeric())


s <- stack()
for (i in 2:length(lis)){

char1 <- noquote(strsplit(as.character(lis[10]), NULL)[[1]])
paste0(s[53], s[54]) 


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
#ov=over(df, x[c(1:15, 39:44)])
ov=over(df, x[vars])
#to data frame
d=as.data.frame(df)
#combine extraction with data
X=cbind(d,  ov)
#copy data frame (only columns of interes)
df=data.frame(y=X[,3], X[,4:dim(X)[2]])
#only oblique coordinates
#remove NAs again
df=na.omit(df)
#df=data.frame(y=X[,3], X[,4:9])
names(df)
print(dim(df))
##you must select the best parameters by tunning them with CV,  the parameter K and the parameter kernel
#test a max value of 25 for k
kmax=25
#tune the model
#df <- df[1:100,]
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
print(CronometroON())
#x <- x[vars][1:100,]
mejorKKNN <- kknn(y~.,train=df,test=x,kernel=unlist(knnTuning$best.parameters[1]),
                        scale=TRUE,k=as.numeric(knnTuning$best.parameters[2]))
print(CronometroOFF())
#get the fitted values for all the area (x)
x$kknn=mejorKKNN$fitted.values
#coordinates(kknn)=~x+y
#gridded(kknn) = TRUE
#kknn=as(kknn,'SpatialPixelsDataFrame')
r <- raster(x['kknn'])
#save the prediction raster (check year/month )
writeRaster(r, file=#index
out <- paste0(probdef, unlist(strsplit(lis[i], '//'))[2])
, overwrite=TRUE) ### 
#save month, replace by year if using yearly averages
yearmo=yearmo_init+i
#best k value
k=as.numeric(knnTuning$best.parameters[2])
#store in the results data frame (month), replace by year if using yearly averages
results[i,1]<-yearmo
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
write.csv(results, file=res_out )
#total time 
CronometroOFF()
#end


