#####################################################
## merge worldpop denominators with census and ACS ##
#####################################################

library(reshape2)
library(data.table)

yrs<-2008:2010

## load in both the acs/census combined data and the worldpop data ##
load('merged_denom_cov_data_fullpop.RData')
acs_list<-adat[2:4]
load('worldpop_ctagg_strat.RData')

## further process worldpop data ##
## remove unnecessary columns ##
wp_list<-lapply(poplist, function(x) x[,c(4,grep('wp_',names(x)))])
## reshape stratified worldpop data ##
rshp_wp<-function(x){
  temp<- reshape2::melt(x,id.vars='GEOID10',value.name = 'wp_pop')
  temp$variable<-as.character(temp$variable)
  temp$sex<-substr(temp$variable,start=4,stop=4)
  temp$agecat<-substr(temp$variable,start=6,stop=nchar(temp$variable))
  
  ## make a coded dataset with new age categories ##
  age_convert<-data.frame(agecat=c("0","1","5","10","15","20","25","30","35","40","45","50","55","60","65","70","75","80"),
                          agenew=paste0('Age',c("0-4","0-4","5-9","10-14","15-19","20-24","25-29",
                                                "30-34","35-44","35-44","45-54","45-54","55-64","55-64",'65-74','65-74','75+','75+')))

  temp3<-merge(temp,age_convert,by='agecat')
  temp4<-aggregate(wp_pop~GEOID10+sex+agenew,data=temp3,FUN=sum,na.rm=T)
  names(temp4)[c(1,3)]<-c('GEOID','agecat')
  
  return(temp4)
}
wp_list<-lapply(wp_list,rshp_wp)

## make race-specific worldpop CT population counts based on ACS CT demographic compositions ##
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

rs_pop<-adat

names(rs_pop)<-yrs

save(rs_pop,file='merged_ce_acs_wp_fullpop.RData')