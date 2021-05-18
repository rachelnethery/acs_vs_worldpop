## pseudo-simulations using 2010 ACS data ##
## Generate data with "true" census pop and use ACS pop or worldpop estimate in models ##
## created by Rachel Nethery ##
## date: 3/1/20 ##

## read command line arguments ##
args<-commandArgs(TRUE)
for (i in 1:length(args)) { eval (parse (text = args[[i]] )) }

setwd(wd)

library(maptools)
library(spdep)
library(MASS)
library(msm)
library(tigris)
library(CARBayes)
library(sp)
library(reshape2)

#reps<-100
#nburn<-10000
#nsamp<-30000

if (yr==2008){
  ind<-1
} else if (yr==2009){
  ind<-2
} else{
  ind<-3
}

####################################################
## 1. population size and covariate info from ACS ##
####################################################

## load the 2010 data to get %black to be used as covariate in analyses for all years ##
x_2010<-read.csv('pmrates_strat_2010.csv',stringsAsFactors = F)

## extract number black per CT ##
n_black<-dcast(subset(x_2010,race=='black' | race=='total',select=c('GEOID','race','acs_1yr_pop')),GEOID~race)
n_black$black_pct<-n_black$black/n_black$total

## load the pre-prepped MA CT-level ACS/census/dp data, subset to only total population (not race-stratified) ##
adat<-subset(read.csv(paste0('pmrates_strat_',yr,'.csv'),stringsAsFactors = F),
                   race=='total')

## merge in black count/percent ##
adat<-merge(adat,n_black,by='GEOID')

## remove CTs with 0 pop ##
adat<-subset(adat,acs_1yr_pop>0 & ce_1yr_pop>0 & wp_1yr_pop>0)

## remove missing covariates ##
adat<-adat[complete.cases(adat[,c('black_pct','acs_pov_pct')]),]

## center/scale variables ##
adat$acs_pov_pct<-as.numeric(scale(adat$acs_pov_pct/100))
adat$black_pct<-as.numeric(scale(adat$black_pct))

adat<-adat[order(adat$GEOID),]

X<-cbind(1,adat$black_pct,adat$acs_pov_pct)
colnames(X)<-c('intercept','bw','pov')
P1<-matrix(adat$acs_1yr_istexp)
P2<-matrix(adat$wp_1yr_istexp)
Ptrue<-matrix(adat$ce_1yr_istexp)

N<-nrow(X)
Nct<-length(unique(adat$GEOID))

###########################################
## 2. adjacency matrix for census tracts ##
###########################################

## extract shapefile ##
ma_shp<-as_Spatial(tracts(state = 'MA',year=2010))

## get adjacency matrix ##
foo = poly2nb(ma_shp, queen=TRUE, row.names=ma_shp@data$GEOID10)
W=nb2mat(foo,style='B')

##########################################################################
## 3. align ordering of adjacency matrix and covariates/population data ##
##########################################################################

## census tracts in W with matches in the data ##
tractrm<-which(rownames(W) %in% as.character(adat$GEOID))

W<-W[tractrm,tractrm]

W<-W[order(as.numeric(rownames(W))),order(as.numeric(rownames(W)))]

identical(unique(as.character(adat$GEOID)),rownames(W))

########################################
## 4. create simulated data structure ##
########################################

## spatial information ##
Dw<-diag(rowSums(W))
#rho_bds<-eigen(sqrt(solve(Dw))%*%W%*%sqrt(solve(Dw)))$values
#print(1/min(rho_bds))
#print(1/max(rho_bds))
rho<-.2

## lambda=scalar background rate ##
#lambda<-.03

## beta coefficients for covariates ##
beta<-c(0,0.1,0.2)

###################################
## 5. set up storage for results ##
###################################

car1<-matrix(NA,nrow=reps,ncol=length(beta)*3+6)
car2<-matrix(NA,nrow=reps,ncol=length(beta)*3+6)

########################
## 6. run simulations ##
########################

set.seed(simnum)

for (xx in 1:reps){
  
  ## U (spatial component) ##
  U<-matrix(mvrnorm(n=1,mu=rep(0,Nct),Sigma=solve(Dw-(rho*W))))
  U<-U[as.numeric(as.factor(adat$GEOID)),]
  
  ## V (unsturctured variance component) ##
  V<-matrix(rnorm(n=N,mean=0,sd=.5))
  
  ## generate outcome ##
  theta<-X%*%beta+U+V
  
  Y<-matrix(rpois(n=N,lambda=Ptrue*exp(theta)))
  
  ####################################
  ## FIT STANDARD SPATIAL CAR MODEL ##
  ####################################
  adat$Y<-Y
  adat$offs1<-log(P1)
  adat$offs2<-log(P2)
  car_offs1<-S.CARbym(Y~black_pct+acs_pov_pct+offset(offs1),family="poisson",data=adat,
                                W=W,burnin=nburn,n.sample=nsamp)
  
  car_offs2<-S.CARbym(Y~black_pct+acs_pov_pct+offset(offs2),family="poisson",data=adat,
                             W=W,burnin=nburn,n.sample=nsamp)
  
  car1[xx,]<-c(apply(car_offs1$samples$beta,2,mean),
               apply(car_offs1$samples$beta,2,quantile,.025),
               apply(car_offs1$samples$beta,2,quantile,.975),
               apply(car_offs1$samples$tau2,2,mean),
               apply(car_offs1$samples$tau2,2,quantile,.025),
               apply(car_offs1$samples$tau2,2,quantile,.975),
               apply(car_offs1$samples$sigma2,2,mean),
               apply(car_offs1$samples$sigma2,2,quantile,.025),
               apply(car_offs1$samples$sigma2,2,quantile,.975))
  
  car2[xx,]<-c(apply(car_offs2$samples$beta,2,mean),
               apply(car_offs2$samples$beta,2,quantile,.025),
               apply(car_offs2$samples$beta,2,quantile,.975),
               apply(car_offs2$samples$tau2,2,mean),
               apply(car_offs2$samples$tau2,2,quantile,.025),
               apply(car_offs2$samples$tau2,2,quantile,.975),
               apply(car_offs2$samples$sigma2,2,mean),
               apply(car_offs2$samples$sigma2,2,quantile,.025),
               apply(car_offs2$samples$sigma2,2,quantile,.975))
  
  
}

car1<-as.data.frame(car1)
car2<-as.data.frame(car2)

names(car1)<-c(paste0('beta',0:(length(beta)-1)),paste0('beta',0:(length(beta)-1),'_lo'),
               paste0('beta',0:(length(beta)-1),'_hi'),'tau2','tau2_lo','tau2_hi',
               'sigma2','sigma2_lo','sigma2_hi')

names(car2)<-c(paste0('beta',0:(length(beta)-1)),paste0('beta',0:(length(beta)-1),'_lo'),
               paste0('beta',0:(length(beta)-1),'_hi'),'tau2','tau2_lo','tau2_hi',
               'sigma2','sigma2_lo','sigma2_hi')

save(car1,car2,file=paste0('results/sim_total_',yr,'_',simnum,'.RData'))
