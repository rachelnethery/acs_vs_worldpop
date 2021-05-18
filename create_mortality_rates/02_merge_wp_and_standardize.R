##########################################################################
## merge the worldpop data for ages <65 with the census and ACS for <65 ##
## then add in premature mortality counts and age/sex standardize       ##
##########################################################################

library(maptools)
library(spdep)
library(MASS)
library(msm)
library(tigris)
library(CARBayes)
library(ggplot2)
library(sf)
library(gtools)
library(reshape2)
library(rstudioapi)

###########################################################
## 1. process worldpop and acs population size estimates ##
###########################################################

load('merged_denom_cov_data.RData')
## remove 07 data due to spatial misalignment ##
acs_list<-adat[2:4]
## remove unnecessary columns ##
acs_list<-lapply(acs_list,
                 function(x) subset(x,select=c("GEOID","race","agecat","sex","acs_pop","ce_pop","acs_pov_pct")))
## aggregate age group 18-19 with 15-17 ##
agg_teen<-function(x){
  acs_cov <- aggregate(acs_pov_pct~GEOID,data=x,
                       FUN = mean)
  ## make a coded dataset with new age categories ##
  age_convert<-data.frame(agecat=c("Age0-4","Age5-9","Age10-14","Age15-17",
                                  "Age18-19","Age20-24","Age25-29",
                                  "Age30-34","Age35-44","Age45-54","Age55-64"),
                          agenew=c("Age0-4","Age5-9","Age10-14","Age15-19",
                                            "Age15-19","Age20-24","Age25-29",
                                            "Age30-34","Age35-44","Age45-54","Age55-64"))
  temp<-merge(x,age_convert,by='agecat')
  temp2<-aggregate(cbind(acs_pop,ce_pop)~GEOID+race+agenew+sex,data=temp,FUN=sum,na.rm=T)
  names(temp2)[3]<-'agecat'
  temp3<-merge(temp2,acs_cov,by='GEOID')
  return(temp3)
}

acs_list<-lapply(acs_list, agg_teen)

load('worldpop_ctagg_strat.RData')
## remove unnecessary columns ##
wp_list<-lapply(poplist, function(x) x[,c(4,grep('wp_',names(x)))])
## reshape stratified worldpop data ##
rshp_wp<-function(x){
  temp<- reshape2::melt(x,id.vars='GEOID10',value.name = 'wp_pop')
  temp$variable<-as.character(temp$variable)
  temp$sex<-substr(temp$variable,start=4,stop=4)
  temp$agecat<-substr(temp$variable,start=6,stop=nchar(temp$variable))
  temp2<-subset(temp,agecat<65)
  
  ## make a coded dataset with new age categories ##
  age_convert<-data.frame(agecat=c("0","1","5","10","15","20","25","30","35","40","45","50","55","60"),
                          agenew=paste0('Age',c("0-4","0-4","5-9","10-14","15-19","20-24","25-29",
                            "30-34","35-44","35-44","45-54","45-54","55-64","55-64")))
  temp3<-merge(temp2,age_convert,by='agecat')
  temp4<-aggregate(wp_pop~GEOID10+sex+agenew,data=temp3,FUN=sum,na.rm=T)
  names(temp4)[c(1,3)]<-c('GEOID','agecat')
  
  return(temp4)
}
wp_list<-lapply(wp_list,rshp_wp)

## make age- gender- and race-specific worldpop CT population counts based on ACS CT demographic compositions ##
adat<-list()
for (i in 1:length(acs_list)){
  ## determine the proportion of each CT population that belongs to each race group from ACS ##
  acs_cov <- aggregate(acs_pov_pct~GEOID,data=acs_list[[i]],
                       FUN = mean)
  
  totrace<-subset(acs_list[[i]],race=='total',select=c('GEOID','agecat','sex','acs_pop'))
  names(totrace)[4]<-'tot_axs_gp'

  temp<-merge(acs_list[[i]],totrace,by=c('GEOID','agecat','sex'))
  temp$race_prop<-temp$acs_pop/temp$tot_axs_gp
  ## when denom is zero for a group, assign race proportions based on the state totals ##
  temp<-temp[order(temp$GEOID,temp$agecat,temp$sex,temp$race),]
  aggtab<-aggregate(acs_pop~race,data=acs_list[[i]],FUN=sum,na.rm=T)
  aggtab$acs_pop<-aggtab$acs_pop/aggtab$acs_pop[which(aggtab$race=='total')]
  names(aggtab)[2]<-'race_prop'
  temp$gp_id<-paste0(temp$GEOID,temp$agecat,temp$sex)
  for (j in unique(temp$gp_id)){
    if (sum(temp$tot_axs_gp[which(temp$gp_id==j)])==0){
      temp$race_prop[which(temp$gp_id==j)]<-aggtab$race_prop
    }
  }
  
  ## merge in the worldpop population estimates ##
  temp<-merge(temp,wp_list[[i]],by=c('GEOID','agecat','sex'))
  temp$wp_pop<-temp$wp_pop*temp$race_prop
  temp<-temp[,c('GEOID','race','agecat','sex','acs_pop','ce_pop','wp_pop')]
  
  ## merge covariate data back in ##
  temp<-merge(temp,acs_cov,by='GEOID')
  
  adat<-c(adat,list(temp))
}

names(adat)<-2008:2010


#############################################
## 2. merge in premature mortality counts  ##
#############################################

## make a function to process mortality data ##
process_mort<-function(mort){
  ## remove deaths that are not assigned to a CT or gender or race or ethnicity ##
  mort<-subset(mort,areakey10>0 & !(racecat=='unknown') & !(gencat=='unknow') & !(hispanic==9))
  
  ## subset to premature mortalities (<65) ##
  pmort<-subset(mort,age<65)
  
  ## add age categories ##
  agecat<-paste0('Age',c('0-4','5-9','10-14','15-19','20-24','25-29','30-34','35-44','45-54','55-64'))
  
  pmort$agecat<-cut(pmort$age,breaks=c(0,5,10,15,20,25,30,35,45,55,65),right = F)
  levels(pmort$agecat)<-agecat
  
  ## count deaths within CT, age, sex groups (for total across races) ##
  cpmort<-as.data.frame(table(pmort$areakey10,pmort$agecat,pmort$gencat))
  names(cpmort)<-c('GEOID','agecat','sex','ndeaths')
  cpmort$race<-'total'
  
  ## count deaths within CT, age, sex groups for hispanic ethnicity ##
  pmort_eth<-subset(pmort,hispanic %in% 1:7 | race==15)
  
  cpmort_eth<-as.data.frame(table(pmort_eth$areakey10,pmort_eth$agecat,pmort_eth$gencat))
  names(cpmort_eth)<-c('GEOID','agecat','sex','ndeaths')
  cpmort_eth$race<-'hispanic'
  
  ## subset to only races non-hispanic white, black, asian/pacislander, american indian ##
  pmort$race2<-NA
  pmort$race2[which(pmort$hispanic==0 & pmort$race==1)]<-'white'
  pmort$race2[which(pmort$race==2)]<-'black'
  pmort$race2[which(pmort$race==3)]<-'am_indian'
  pmort$race2[which(pmort$race %in% c(4, 5, 7, 8, 9, 10, 12))]<-'asian'
  pmort$race2[which(pmort$race %in% c(6,11))]<-'pac_islander'
  
  pmort_sub<-subset(pmort,!is.na(race2))
  
  ## count deaths within CT, race, age, sex groups ##
  cpmort_sub<-as.data.frame(table(pmort_sub$areakey10,pmort_sub$race2,pmort_sub$agecat,pmort_sub$gencat))
  names(cpmort_sub)<-c('GEOID','race','agecat','sex','ndeaths')
  
  adat_pmort<-rbind(cpmort,cpmort_sub,cpmort_eth)
  levels(adat_pmort$sex)<-c('f','m')
  
  return(adat_pmort)
}

## read/process 2010 mortality data ##
mort1yr<-read.csv('mort_geo_10.csv',
                  stringsAsFactors = F,header=T)

adat_mort1yr<-process_mort(mort=mort1yr)

adat<-lapply(adat,function(x) merge(x,adat_mort1yr,by=c('GEOID','race','agecat','sex'),all.x=T))

## missing ndeaths ==> 0 deaths ##
rr<-function(x){
  x$ndeaths[which(is.na(x$ndeaths))]<-0
  return(x)
}
adat<-lapply(adat,FUN = rr)

save(adat,file='merged_pm_ce_acs_wp_strat.RData')

#######################################################
## 3. age standardization by each population dataset ##
#######################################################

adat<-lapply(adat,function(x) subset(x,race %in% c('total','white','black','hispanic')))

## get age-stratified mortality counts for standardization ##
mort_ref<-aggregate(ndeaths~agecat+sex,data=subset(adat[[1]],race=='total'),FUN=sum,na.rm=T)
names(mort_ref)[3]<-'deaths_agexsex'

agestand<-function(dat,popnm,incnm,mort_ref,d,nyrs){
  datnew<-dat[,c('GEOID','race','agecat','sex',popnm,incnm)]
  datnew[[popnm]]<-datnew[[popnm]]*nyrs
  dattot<-datnew[which(datnew$race=='total'),]
  
  ## get pop counts in each age/sex group for standardization ##
  pop_ref<-aggregate(dattot[[popnm]],by=list(dattot$agecat,dattot$sex),FUN=sum,na.rm=T)
  names(pop_ref)<-c('agecat','sex','pop_agexsex')
  
  ## get reference mortality rates ##
  mort_ref<-merge(mort_ref,pop_ref,by=c('agecat','sex'))
  mort_ref$mr<-mort_ref$deaths_agexsex/mort_ref$pop_agexsex
  mort_ref$pop_prop<-mort_ref$pop_agexsex/sum(mort_ref$pop_agexsex,na.rm=T)
  
  ## make zero populations missing (otherwise we get infinite rates) ##
  datnew[[popnm]][which(datnew[[popnm]]==0)]<-NA

  ## data to compute crude rates by CT and race ##
  crudedat<-aggregate(list(datnew[[incnm]],datnew[[popnm]]),by=list(datnew$GEOID,datnew$race),FUN=sum,na.rm=T)
  names(crudedat)<-c('GEOID','race','ndeaths','pop')
  ## compute crude rates ##
  crudedat$pmrate_cr<-crudedat$ndeaths/crudedat$pop
  
  ## data to do direct age standardization ##
  dstdat<-merge(datnew,mort_ref,by=c('agecat','sex'))
  dstdat$pmrate_dst<-(dstdat[[incnm]]/dstdat[[popnm]])*dstdat[['pop_prop']]
  dstdat<-aggregate(pmrate_dst~GEOID+race,data=dstdat,FUN=sum,na.rm=T)
  
  ## data to do indirect age standardization ##
  ## also add both observed counts and expected (separately) to the data ##
  istdat<-merge(datnew,mort_ref,by=c('agecat','sex'))
  istdat$istexp<-(istdat[[popnm]]*istdat[['mr']])
  istdat<-aggregate(list(istdat$ndeaths,istdat$istexp),by=list(istdat$GEOID,istdat$race),FUN=sum,na.rm=T)
  names(istdat)<-c('GEOID','race','ndeaths','istexp')
  istdat$pmrate_ist<-istdat$ndeaths/istdat$istexp

  ## merge crude and standardized rates ##
  m1<-merge(crudedat,dstdat,by=c('GEOID','race'),all.x=T)
  m2<-merge(m1,istdat[,-3],by=c('GEOID','race'),all.x=T)
  
  ## make NA rates in places with 0 population ##
  m2[which(m2$pop==0),5:ncol(m2)]<-NA
  
  names(m2)[3:ncol(m2)]<-c(paste0('ndeaths_',nyrs,'yr'),paste0(d,'_',nyrs,'yr_',names(m2)[4:ncol(m2)]))

  return(m2)
}

apply_agestand<-function(x){
  foo1<-agestand(dat=x,popnm='acs_pop',incnm='ndeaths',mort_ref=mort_ref,d='acs',nyrs=1)
  foo2<-agestand(dat=x,popnm='wp_pop',incnm='ndeaths',mort_ref=mort_ref,d='wp',nyrs=1)
  foo3<-agestand(dat=x,popnm='ce_pop',incnm='ndeaths',mort_ref=mort_ref,d='ce',nyrs=1)
  foo4<-merge(foo1,foo2[,-3],by=c('GEOID','race'), all=T)
  foo5<-merge(foo4,foo3[,-3],by=c('GEOID','race'), all=T)
  return(foo5)
}

pop_as<-lapply(adat, apply_agestand)

#########################################
## 4. merge back in the covariate data ##
#########################################

covdat<-adat[[1]][!duplicated(adat[[1]][,c('GEOID','race')]),c('GEOID','race','acs_pov_pct')]

ss<-function(x){
  foo1<-merge(x,covdat,by=c('GEOID','race'))
  ## sort by geoid ##
  foo2<-foo1[order(foo1$GEOID),]
  
  return(foo2)
}

adat<-lapply(pop_as,ss)

save(adat,file='merged_pmrates_ce_acs_wp_strat.RData')
