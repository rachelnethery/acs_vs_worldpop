####################################################
## THIS CODE PROCESSES RAW WORLDPOP DATA FROM URL ##
## NEEDS TO BE RUN ON CANNON IN PARALLEL          ##
####################################################

## at command line run the following for interactive R job
## module load gcc/7.1.0-fasrc01 R/3.3.3-fasrc01 udunits/2.2.26-fasrc01 gdal/2.3.0-fasrc01 proj/5.0.1-fasrc01 geos/3.6.2-fasrc01
## export R_LIBS_USER=$HOME/apps/R:$R_LIBS_USER
## srun -p test --pty --mem 10000 -t 0-02:00 /bin/bash
## R --quiet

## read command line arguments ##
args<-commandArgs(TRUE)
for (i in 1:length(args)) { eval (parse (text = args[[i]] )) }

library(raster)
library(rgdal)

## years of data ##
yrs<-as.character(c(2008:2010))
ages<-as.character(c(0,1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80))
sexes<-c('m','f')
kk<-expand.grid(yrs,sexes,ages,stringsAsFactors = F)
kk_i<-kk[m,]
print(kk_i)

setwd('/n/holyscratch01/dominici_lab/rachel/waller/worldpop/')

if (file.exists(paste0('wp_',kk_i[1],'_',kk_i[2],'_',kk_i[3],'.tif'))){
}else{

## MA lat/long boundaries ##
xmin<-(-73.508142)
ymin<-41.237964
xmax<-(-69.928393)
ymax<-42.886589

cropbox<-extent(xmin,xmax,ymin,ymax)

## read in raster ##
rname<-paste0('https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020/',
              kk_i[1],'/USA/usa_',kk_i[2],'_',kk_i[3],'_',kk_i[1],'.tif')
download.file(url=rname,destfile=paste0('/n/holyscratch01/dominici_lab/rachel/waller/worldpop/raw',
                                        kk_i[2],'_',kk_i[3],'_',kk_i[1],'.tif'))
#download.file(url=rname,destfile=paste0('raw',kk_i[2],'_',kk_i[3],'_',kk_i[1],'.tif'))

wp<-raster(paste0('/n/holyscratch01/dominici_lab/rachel/waller/worldpop/raw',
                  kk_i[2],'_',kk_i[3],'_',kk_i[1],'.tif'))

## view attributes ##
#wp

## crop to MA boundaries ##
ma_pop<-crop(wp,cropbox)

## export MA raster ##
writeRaster(ma_pop, paste0('wp_',kk_i[1],'_',kk_i[2],'_',kk_i[3],'.tif'), overwrite=TRUE)

## delete big raster file ##
file.remove(paste0('/n/holyscratch01/dominici_lab/rachel/waller/worldpop/raw',
                   kk_i[2],'_',kk_i[3],'_',kk_i[1],'.tif'))

}