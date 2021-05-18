#########################################################
## CODE TO AGGREGATE WORLDPOP DATA TO MA CENSUS TRACTS ##
#########################################################

library(sf)
library(raster)
library(sp)
library(ggplot2)
library(USAboundaries)
library(tigris)
options(tigris_use_cache = FALSE)

yrs<-c(2008:2010)
ages<-c(0,1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80)

## extract MA CT shapefile ##
ma_shp<-tracts(state = 'MA',year=2010)
ma_shp<-spTransform(ma_shp, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

poplist<-list()
for (i in yrs){
  ctpop<-ma_shp@data
  for (j in c('m','f')){
    for (k in ages){
      
      ## read in population raster ##
      ma07<-raster(paste0('wp_',i,'_',j,'_',k,'.tif'))
      
      ## convert to spatial points data frame ##
      ma07point <- rasterToPoints(ma07,spatial = T)
      
      ## sum of population at all points within each MA CT ##
      sum_overct<-over(ma_shp,ma07point,fn=sum)
      names(sum_overct)<-paste0('wp_',j,'_',k)
      
      ctpop<-cbind(ctpop,sum_overct)
      
    }
  }
  poplist<-c(poplist,list(ctpop))
}

names(poplist)<-yrs

save(poplist,file='worldpop_ctagg_strat.RData')
