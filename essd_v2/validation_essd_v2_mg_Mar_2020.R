library(raster)
library(dplyr)
library(openair)
library(Metrics)

library(scales)

valset <- readRDS('ALL_ismn.rds')
valset[valset$sm<=0,] <- NA
valset <- na.omit(valset)
valset[valset$sm>=1,] <- NA
valset <- na.omit(valset)
summary(valset)

all <- list.files(pattern='tif')[2:27]
all <- stack(all)
library(spatialEco)
#all <- raster.transformation(all, trans = "stretch", smin = min(valset$sm), smax = max(valset$sm))
res <- data.frame()
for (i in 1991:2016){
val <- valset[grep(i,valset$date),]
val <-data.frame(val)
val$e <- extract(all[[i-1990]], val[c('x', 'y')])
res <- rbind(res, val)
}
all_covs <- na.omit(res)
#all_covs$e <- rescale(all_covs$e, c(min(all_covs$sm), max(all_covs$sm)))

topo <- stack(list.files('/home/mguevara/Downloads/predicted-2001-2016-esa-sm-topo-GLOBALmean/', full.names=TRUE
, pattern='tif'))
#topo[topo<0] <- NA
#topo[topo>1] <- NA
library(spatialEco)
#topo <- raster.transformation(topo, trans = "stretch", smin = min(valset$sm), smax = max(valset$sm))
res <- data.frame()
for (i in 1991:2016){
val <- valset[grep(i,valset$date),]
val <-data.frame(val)
val$e <- extract(topo[[i-1990]], val[c('x', 'y')])
res <- rbind(res, val)
}
topo <- na.omit(res)
topo$e[topo$e < 0] <- NA
topo <- na.omit(topo)
#topo$e <- rescale(topo$e, c(min(topo$sm), max(topo$sm)))

res <- data.frame()
esa <- stack(list.files('/home/mguevara/Downloads/annual/', full.names=TRUE))
library(spatialEco)
#esa <- raster.transformation(esa, trans = "stretch", smin = min(valset$sm), smax = max(valset$sm))
res <- data.frame()
for (i in 1991:2016){
val <- valset[grep(i,valset$date),]
val <-data.frame(val)
val$e <- extract(esa[[i-1990]], val[c('x', 'y')])
res <- rbind(res, val)
}

esa <- na.omit(res)
#esa$e <- rescale(esa$e, c(min(esa$sm), max(esa$sm)))

all_covs$yr <-  as.numeric(substr(all_covs$date, 1, 4))
ALL=all_covs %>% group_by(yr) %>% summarise_each(list(mean = mean))
ALL$e <- rescale(ALL$e_mean, c(min(ALL$sm_mean), max(ALL$sm_mean)))
cor(ALL$sm_mean, ALL$e)
rmse(ALL$sm_mean, ALL$e)

ISMN <- data.frame(date=as.Date(paste0(ALL$yr, '/01/01/')), sm = ALL$sm_mean)
TheilSen(ISMN, 'sm', slope.percent = TRUE, alpha = 0.01, ylim=c(0.1, 0.4),  dec.place=3, deseason=TRUE, shade = "white", lab.cex=1.5, data.col='black', text.col='black',  scales=list(tck=c(1,0), x=list(cex=1.5), y=list(cex=1.5)))

scatterPlot(ALL, x = "sm_mean", y = "e",  mod.line=TRUE, smooth=TRUE)
ALL <- data.frame(date=as.Date(paste0(ALL$yr, '/01/01/')), sm = ALL$e)
TheilSen(ALL, 'sm', slope.percent = TRUE, alpha = 0.01, ylim=c(0.1, 0.4),  dec.place=3, deseason=TRUE, shade = "white", lab.cex=1.5, data.col='black', text.col='black',  scales=list(tck=c(1,0), x=list(cex=1.5), y=list(cex=1.5)))

topo$yr <-  as.numeric(substr(topo$date, 1, 4))
TOPO=topo %>% group_by(yr) %>% summarise_each(list(mean = mean))
TOPO$e <- rescale(TOPO$e_mean, c(min(TOPO$sm_mean), max(TOPO$sm_mean)))
 cor(TOPO$sm_mean, TOPO$e)
 rmse(TOPO$sm_mean, TOPO$e)
scatterPlot(TOPO, x = "sm_mean", y = "e",  mod.line=TRUE, linear=TRUE)


TOPO <- data.frame(date=as.Date(paste0(TOPO$yr, '/01/01/')), sm = TOPO$e)
TheilSen(TOPO, 'sm', slope.percent = TRUE, alpha = 0.01, ylim=c(0.1, 0.4),  dec.place=3, deseason=TRUE, shade = "white", lab.cex=1.5, data.col='black', text.col='black',  scales=list(tck=c(1,0), x=list(cex=1.5), y=list(cex=1.5)))

esa$yr <-  as.numeric(substr(esa$date, 1, 4))
ESA=esa %>% group_by(yr) %>% summarise_each(list(mean = mean))
ESA$e <- rescale(ESA$e_mean, c(min(ESA$sm_mean), max(ESA$sm_mean)))
cor(ESA$sm_mean, ESA$e)
rmse(ESA$sm_mean, ESA$e)
scatterPlot(ESA, x = "sm_mean", y = "e",  mod.line=TRUE, linear=TRUE)
ESA <- data.frame(date=as.Date(paste0(ESA$yr, '/01/01/')), sm = ESA$e)
TheilSen(ESA, 'sm', slope.percent = TRUE, alpha = 0.01, ylim=c(0.1, 0.4),  dec.place=3, deseason=TRUE, shade = "white", lab.cex=1.5, data.col='black', text.col='black',  scales=list(tck=c(1,0), x=list(cex=1.5), y=list(cex=1.5)))










ben <- read.csv('/home/mguevara/Downloads/SRDB_V4_1578/data/srdb-data-V4.csv')
prec <- ben[c(2,13:14, 34, 21)]
res <- data.frame()
#ras <- stack(list.files('/home/mguevara/Downloads/annual/', full.names=TRUE))
#ras<- stack(list.files('/home/mguevara/Downloads/predicted-2001-2016-esa-sm-topo-GLOBALmean/', full.names=TRUE, pattern='tif'))
ras <- stack(list.files(pattern='tif')[4:31])
res <- data.frame()
for (i in 1991:2016){
val <- prec[grep(i,prec$Entry_date),]
val <-data.frame(val)
val$e <- extract(ras[[i-1990]], val[c('Longitude', 'Latitude')])
res <- rbind(res, val)
}

all <- res
all$e[all$e < 0] <- NA
all <- na.omit(all)
all$sm <- log1p(all$e)
cor(all$MAP, all$e)
tro <- all[all$Biome=='Tropical',]
bor <- all[all$Biome=='Boreal',]
cor(bor$MAP, bor$e)
cor(tro$MAP, tro$e)




 scatterPlot(all, x = "MAP", y = "sm", smooth = TRUE)

