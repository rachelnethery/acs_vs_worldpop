################################################################################
## extract ACS and census denominator data for ages <65 only using tidycensus ##
## then format and merge together for analyses                                ##
################################################################################

library(tidycensus)
library(tigris)
library(sp)
library(stringr)
library(haven)
library(reshape2)
library(rstudioapi)

#################
## 1. ACS DATA ##
#################

## variable names ##
#v12 <- load_variables(2012, "acs5")
acs_list<-list()
for (yr in 2009:2012){
  ##extract using tidycensus ##
  ma_acs<- get_acs(geography = "tract",
                   variables = c(## total pop <65
                     paste0('B01001_',str_pad(as.character(c(3:19,27:43)),width=3,side='left',pad='0')),
                     ## non-hispanic white <65
                     paste0('B01001H_',str_pad(as.character(c(3:13,18:28)),width=3,side='left',pad='0')),
                     ## black <65
                     paste0('B01001B_',str_pad(as.character(c(3:13,18:28)),width=3,side='left',pad='0')),
                     ## american indian/alaska native <65
                     paste0('B01001C_',str_pad(as.character(c(3:13,18:28)),width=3,side='left',pad='0')),
                     ## asian <65
                     paste0('B01001D_',str_pad(as.character(c(3:13,18:28)),width=3,side='left',pad='0')),
                     ## native hawaiian/pacific islander
                     paste0('B01001E_',str_pad(as.character(c(3:13,18:28)),width=3,side='left',pad='0')),
                     ## hispanic
                     paste0('B01001I_',str_pad(as.character(c(3:13,18:28)),width=3,side='left',pad='0'))
                   ),
                   state = "MA",year=yr)
  
  ## organize ##
  ma_acs<-as.data.frame(ma_acs)
  
  ## racevar is a variable that tells us what racegroup the row represents ##
  ma_acs$racevar<-substr(ma_acs$variable, 1, 7)
  ## vnum is a variable that tells us the sex and age group the row represents ##
  ma_acs$vnum<-as.numeric(substr(ma_acs$variable,nchar(ma_acs$variable)-2,nchar(ma_acs$variable)))
  
  ## create a dataset with proper race and sex variables and an age group variable with consistent groupings (called rxsxa) based on the racevar and vnum variables ##
  ## age ranges ##
  agecat<-paste0('Age',c('0-4','5-9','10-14','15-17','18-19','20-24','25-29','30-34','35-44','45-54','55-64'))
  
  rxsxa<-ma_acs[!duplicated(ma_acs[,c('racevar','vnum')]),c('racevar','vnum')]
  rxsxa<-rxsxa[order(rxsxa$racevar,rxsxa$vnum),]
  
  realign_fineage<-c(1:5,6,6,6,7,8,9,9,10,10,11,11,11)
  rxsxa$age<-c(rep(realign_fineage,2),rep(1:length(agecat),2*(length(unique(rxsxa$racevar))-1)))
  rxsxa$sex<-c(rep(c('m','f'),each=length(realign_fineage)),rep(rep(c('m','f'),each=length(agecat)),length(unique(rxsxa$racevar))-1))
  rxsxa<-merge(rxsxa,data.frame('age'=1:length(agecat),'agecat'=agecat),by='age')
  rxsxa<-merge(rxsxa,data.frame(
    'race'=c('total','white','black','am_indian','asian','pac_islander','hispanic'),
    'racevar'=paste0('B01001',c('_','H','B','C','D','E','I'))),by='racevar')
  
  ## merge this new dataset into ma_acs to add proper race, age, sex variables ##
  ma_acs<-merge(ma_acs,rxsxa,by=c('racevar','vnum'))
  
  ## need to now aggregate population size and MOE where applicable due to the inconsistent age groupings ##
  temp1<-aggregate(estimate~GEOID+race+agecat+sex,data=ma_acs,
                   FUN=sum,na.rm=T)
  names(temp1)[5]<-'acs_pop'
  
  temp2<-aggregate(moe~GEOID+race+agecat+sex,data=ma_acs,
                   FUN = function(x) 1.96*(sqrt(sum((x/1.645)^2,na.rm=T))))
  names(temp2)[5]<-'acs_moe'
  
  ma_acs<-merge(temp1,temp2,by=c('GEOID','race','agecat','sex'))
  
  ma_acs<-ma_acs[order(ma_acs$GEOID),]
  
  ## add in the poverty and ICE measures from Pam computed from ACS ##
  ice<-read_sas(data_file = 'icemeasures_acs_0812_12_12_18.sas7bdat')
  ice<-as.data.frame(ice)
  ice<-ice[,c('GEOid2','perc_belowpov','ICEracewb','ICEinc','ICEwbinc')]
  names(ice)<-c('GEOID','acs_pov_pct','acs_ice_racewb','acs_ice_inc','acs_ice_raceinc')
  
  ma_acs<-merge(ma_acs,ice,by='GEOID')
  
  acs_list<-c(acs_list,list(ma_acs))
  
}

save(acs_list,file='acs_data.RData')
rm(list=ls())

####################
## 2. CENSUS DATA ##
####################

## census 2010 data ##
## variable names ##
#v10 <- load_variables(2010, "sf1")
racemg<-data.frame('race'=c('total','white','black','am_indian','asian','pac_islander','hispanic'),'racevar'=paste0('P012',c('0','I','B','C','D','E','H')))
nums65<-str_pad(c(3:19,27:43),width=3,side='left',pad='0')
vnames<-c(paste0('P012',nums65),paste0(rep(racemg$racevar[2:nrow(racemg)],each=length(nums65)),nums65))

##extract using tidycensus ##
ma_ce<-get_decennial(geography = "tract",
                     variables=vnames,
                     state = "MA",year=2010,sumfile='sf1')

## organize ##
ma_ce<-as.data.frame(ma_ce)

## racevar is a variable that tells us what racegroup the row represents ##
ma_ce$racevar<-substr(ma_ce$variable,1,5)
## vnum is a variable that tells us the sex and age group the row represents ##
ma_ce$vnum<-as.numeric(substr(ma_ce$variable,nchar(ma_ce$variable)-2,nchar(ma_ce$variable)))

## create a dataset with proper race and sex variables and an age group variable with consistent groupings (called rxsxa) based on the racevar and vnum variables ##
## age ranges ##
agecat<-paste0('Age',c('0-4','5-9','10-14','15-17','18-19','20-24','25-29','30-34','35-44','45-54','55-64'))

rxsxa<-ma_ce[!duplicated(ma_ce[,c('racevar','vnum')]),c('racevar','vnum')]
rxsxa<-rxsxa[order(rxsxa$racevar,rxsxa$vnum),]

realign_fineage<-c(1:5,6,6,6,7,8,9,9,10,10,11,11,11)
rxsxa$age<-rep(realign_fineage,2*(length(unique(rxsxa$racevar))))
rxsxa$sex<-rep(rep(c('m','f'),each=length(realign_fineage)),length(unique(rxsxa$racevar)))
rxsxa<-merge(rxsxa,data.frame('age'=1:length(agecat),'agecat'=agecat),by='age')
rxsxa<-merge(rxsxa,racemg,by='racevar')

ma_ce<-merge(ma_ce,rxsxa,by=c('racevar','vnum'))

## need to now aggregate population size and MOE where applicable due to the inconsistent age groupings ##
ma_ce<-aggregate(value~GEOID+race+agecat+sex,data=ma_ce,
                 FUN=sum,na.rm=T)
names(ma_ce)[5]<-'ce_pop'

## merge in the ice for race from the census (P003002-P003003)/P001001 ##
ice<-get_decennial(geography = "tract",
                   variables=c('P001001','P003002','P003003'),
                   state = "MA",year=2010,sumfile='sf1',output = 'wide')
ice<-as.data.frame(ice)
ice<-data.frame('GEOID'=ice$GEOID,'ce_ice_racewb'=(ice$P003002-ice$P003003)/ice$P001001)

ma_ce<-merge(ma_ce,ice,by='GEOID')

ma_ce<-ma_ce[order(ma_ce$GEOID),]

save(ma_ce,file='ce_data.RData')

rm(list=ls())

##################################
## 3. merge ACS and census data ##
##################################

## read each processed dataset ##
load('acs_data.RData')

load('ce_data.RData')

## merge them by fips code, race, age, and sex ##
adat<-lapply(acs_list,FUN=function(x) merge(x,ma_ce,by=c('GEOID','race','agecat','sex')))

save(adat,file='merged_denom_cov_data.RData')
