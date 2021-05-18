####################################################################
## visualize trends over time (2008-2010) in ACS and WP estimates ##
####################################################################

library(data.table)
library(ggplot2)
library(GGally)
library(cowplot)
library(reshape2)
library(plyr)
library(RColorBrewer)
library(tidycensus)
library(ggthemr)

ggthemr('flat')

load('merged_ce_acs_wp_fullpop.RData')

adat<-rs_pop

## aggregate across age and sex groups ##
rs_pop<-lapply(rs_pop,function(x) aggregate(cbind(acs_pop,ce_pop,wp_pop)~GEOID+race,data=x,FUN=sum,na.rm=T))

rs_pop<-lapply(rs_pop,FUN=setDT)

nyrs<-length(rs_pop)

data(fips_codes)
setDT(fips_codes)

################################################################
## create a function that formats data to facilitate plotting ##
################################################################

create_dat<-function(rs_pop,vnm,pnm,racegp,brks,labs){
  
  ## data structure for scatterplot matrices ##
  for_spm<-merge(merge(rs_pop[[1]][race==racegp,.(GEOID,est08=get(vnm))],
                       rs_pop[[2]][race==racegp,.(GEOID,est09=get(vnm))],by='GEOID'),
                 rs_pop[[3]][race==racegp,.(GEOID,est10=get(vnm),census10=ce_pop)],by='GEOID')
  
  for_spm<-for_spm[census10>0]
  
  for (i in 2:4){
    for (j in (i+1):5){
      for_spm[,(paste0('diff_',substr(names(for_spm)[j],1,1),
                       substr(names(for_spm)[j],nchar(names(for_spm)[j])-1,nchar(names(for_spm)[j])),
                       substr(names(for_spm)[i],1,1),
                       substr(names(for_spm)[i],nchar(names(for_spm)[i])-1,nchar(names(for_spm)[i])))):=get(names(for_spm)[j])-get(names(for_spm)[i])]
    }
  }
  
  ## data structure for histograms ##
  for_hists<-melt(for_spm[,!(est08:census10)],id.vars = 'GEOID')
  setDT(for_hists)
  for_hists[,variable:=as.character(variable)
            ][,base_yr:=substr(variable,start=nchar(variable)-2,stop=nchar(variable))
              ][,comp_yr:=substr(variable,start=nchar(variable)-5,stop=nchar(variable)-3)
                ][,':='(base_yr=revalue(base_yr,c(e08=paste(pnm,'08'),e09=paste(pnm,'09'),e10=paste(pnm,'10'))),comp_yr=revalue(comp_yr,c(e09=paste(pnm,'09'),e10=paste(pnm,'10'),c10='Census 10')))]
  
  ## data structure for spaghetti plots ##
  for_tmp<-melt(for_spm[,.(GEOID,est08,est09,est10)],id.vars='GEOID')
  setDT(for_tmp)
  for_tmp[,fips:=substr(GEOID,1,5)
          ][,GEOID:=factor(GEOID,levels=unique(for_spm[order(census10),GEOID]))
            ][,variable:=revalue(variable,c(est08='2008',est09='2009',est10='2010'))]
  
  ma_fips<-fips_codes[state=='MA']
  ma_fips<-ma_fips[,fips:=paste0(state_code,county_code)
                   ][,.(county,fips)]
  
  for_tmp<-merge(for_tmp,ma_fips,by='fips')
  for_tmp<-merge(for_tmp,for_spm[,.(GEOID,bline=census10)],by='GEOID')
  for_tmp[,bline_cat:=cut(bline,breaks=c(0,brks,max(bline)+1),labels=labs,include.lowest=T,right=F)]
  
  
  return(list('for_spm'=for_spm,'for_hists'=for_hists,'for_lng'=for_tmp))
  
}

########################################################
## make plots of ACS/census estimates and time trends ##
########################################################

acs<-create_dat(rs_pop,'acs_pop','ACS',racegp = 'total',
                brks=c(2000,4000,6000,8000),labs=c('<2000','[2000,4000)',
                                                   '[4000,6000)','[6000,8000)','8000+'))

## scatterplot matrices of ACS population size estimates across years 2008-2010 ##
pdf('figures/acs_time_diffs_spm.pdf',width=10,height=6)
ggpairs(acs$for_spm,columns = 2:5)
dev.off()

## histograms of differences in ACS population size estimates across years ##
pdf('figures/acs_time_diffs_hists.pdf',width=6,height=8)
  ggplot(acs$for_hists, aes(x=abs(value),color=comp_yr))+
  geom_freqpoly(binwidth = 50,size=1.5)+
  xlim(0,1000)+
  facet_wrap(base_yr~.,ncol=1,scales = 'free_y')+
  theme_bw()+
  theme(legend.title=element_blank(),text=element_text(size=16,colour='black'),
        axis.title.y=element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())+
    scale_color_brewer(palette="Set1")+
  xlab('Absolute Difference')
dev.off()

mm<-sample(unique(acs$for_lng[,GEOID]),size=length(unique(acs$for_lng[,GEOID])),replace=F)
acs$for_lng[,GEOID:=factor(GEOID,levels=mm)]

## acs spaghetti plots for total pop (figure 3) ##
pdf('figures/acs_longitudinal_plots.pdf',width=12)
ggplot(acs$for_lng, aes(x=variable, y=value, group=GEOID,color=bline_cat)) +
  geom_line(size=.5)+
  facet_wrap(~county,ncol=5,scales='free_y')+
  scale_color_brewer(palette="Spectral")+
  #theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(color="'10 Census Size",y='Population size',x='ACS Center Year')

tot_long<-ggplot(acs$for_lng, aes(x=variable, y=value/bline, group=GEOID,color=bline_cat)) +
  geom_line(size=.4)+
  facet_wrap(~county,ncol=5,scales='free_y')+
  scale_color_brewer(palette="Spectral")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.position = c(0.92, 0.05)) +
  labs(color="'10 Census Size",y="Population size as a proportion of Census '10",x='ACS Center Year')
print(tot_long)
dev.off()

## acs spaghetti plots for black only ##
acs_blk<-create_dat(rs_pop,'acs_pop','ACS',racegp = 'black',
                    brks=c(50,100,500,1000),labs=c('<50','[50,100)',
                                                   '[100,500)','[500,1000)','1000+'))

mm<-sample(unique(acs_blk$for_lng[,GEOID]),size=length(unique(acs_blk$for_lng[,GEOID])),replace=F)
acs_blk$for_lng[,GEOID:=factor(GEOID,levels=mm)]

pdf('figures/acs_longitudinal_plots_blk.pdf',width=12)
ggplot(acs_blk$for_lng, aes(x=variable, y=value, group=GEOID,color=bline_cat)) +
  geom_line(size=.5)+
  facet_wrap(~county,ncol=5,scales='free_y')+
  scale_color_brewer(palette="Spectral")+
  #theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(color="'10 Census Size",y='Population size',x='ACS Center Year')

blk_long<-ggplot(acs_blk$for_lng, aes(x=variable, y=value/bline, group=GEOID,color=bline_cat)) +
  geom_line(size=.4)+
  facet_wrap(~county,ncol=5,scales='free_y')+
  scale_color_brewer(palette="Spectral")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.position = c(0.92, 0.05)) +
  labs(color="'10 Census Size",y="Population size as a proportion of Census '10",x='ACS Center Year')
print(blk_long)
dev.off()

pdf('figures/acs_longitudinal_plots_comb.pdf',width=9,height=10)
plot_grid(tot_long, blk_long, labels = c('A', 'B'),ncol=1)
dev.off()

## acs spaghetti plots for white only ##
acs_wht<-create_dat(rs_pop,'acs_pop','ACS',racegp = 'white',
                    brks=c(500,1000,2000,4000),labs=c('<500','[500,1000)',
                                                      '[1000,2000)','[2000,4000)','4000+'))

pdf('figures/acs_longitudinal_plots_wht.pdf',width=12)
ggplot(acs_wht$for_lng, aes(x=variable, y=value, group=GEOID,color=bline_cat)) +
  geom_line(size=.4)+
  facet_wrap(~county,ncol=5,scales='free_y')+
  scale_color_brewer(palette="Spectral")+
  #theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(color="'10 Census Size",y='Population size',x='ACS Center Year')

ggplot(acs_wht$for_lng, aes(x=variable, y=value/bline, group=GEOID,color=bline_cat)) +
  geom_line()+
  facet_wrap(~county,ncol=5,scales='free_y')+
  scale_color_brewer(palette="Spectral")+
  #theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.position = c(0.92, 0.05)) +
  labs(color="'10 Census Size",y="Population size as a proportion of Census '10",x='ACS Center Year')
dev.off()


#######################################################
## make plots of WP/census estimates and time trends ##
#######################################################

wp<-create_dat(rs_pop,'wp_pop','WP',racegp='total',
               brks=c(2000,4000,6000,8000),labs=c('<2000','[2000,4000)',
                                                  '[4000,6000)','[6000,8000)','8000+'))

## scatterplot matrices of WP population size estimates across years 2008-2010 ##
pdf('figures/wp_time_diffs_spm.pdf',width=10,height=6)
ggpairs(wp$for_spm,columns = 2:5)
dev.off()

## histograms of differences in WP population size estimates across years ##
pdf('figures/wp_time_diffs_hists.pdf',width=6,height=8)
ggplot(wp$for_hists, aes(x=value, color=comp_yr))+
  geom_density()+
  xlim(-2000,2000)+
  facet_wrap(base_yr~.,ncol=1,scales = 'free_y')+
  theme_bw()+
  theme(legend.title=element_blank(),text=element_text(size=16,colour='black'),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())+
  xlab('Difference')
dev.off()

## WP spaghetti plots ##
pdf('figures/wp_longitudinal_plots.pdf',width=12)
ggplot(wp$for_lng, aes(x=variable, y=value, group=GEOID,color=bline_cat)) +
  geom_line()+
  facet_wrap(~county,ncol=5,scales='free_y')+
  scale_color_brewer(palette="Spectral")+
  #theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(color="'10 Census Size",y='Population size',x='WP Year')

ggplot(wp$for_lng, aes(x=variable, y=value/bline, group=GEOID,color=bline_cat)) +
  geom_line(size=.4)+
  facet_wrap(~county,ncol=5,scales='free_y')+
  scale_color_brewer(palette="Spectral")+
  #theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.position = c(0.92, 0.05)) +
  labs(color="'10 Census Size",y="Population size as a proportion of Census '10",x='WP Year')
dev.off()

## plots of association between poverty and errors in ACS/WP by age group (figure 2) ##
kk<-subset(adat[[3]],race=='total')
ll<-aggregate(cbind(acs_pop,ce_pop,wp_pop)~GEOID+agecat,data=kk,FUN=sum,na.rm=T)
mm<-aggregate(acs_pov_pct~GEOID,data=kk,FUN=mean,na.rm=T)
nn<-merge(ll,mm,by='GEOID')
nn$agecat<-factor(nn$agecat,levels=paste0('Age',c('0-4','5-9','10-14','15-19','20-24','25-29','30-34','35-44','45-54','55-64','65-74','75+')))
#tapply(nn$ce_pop-nn$wp_pop,nn$agecat,median)

pdf('figures/age_strat_plots.pdf',width=12)
ggplot(nn,aes(x=ce_pop,y=wp_pop))+
  geom_point()+
  geom_smooth(method='lm')+
  facet_wrap(agecat~.)+
  ylab('WorldPop 2010')+
  xlab('Census 2010')+
  ylim(0,2000)+xlim(0,2000)

ggplot(nn,aes(x=acs_pov_pct,y=ce_pop-wp_pop))+
  geom_point()+
  geom_smooth(method='lm')+
  facet_wrap(agecat~.,scales = 'free_y')+
  ylab('Census-WorldPop')+
  xlab('ACS Percent Poverty')

ggplot(nn,aes(x=ce_pop,y=acs_pop))+
  geom_point()+
  geom_smooth(method='lm')+
  facet_wrap(agecat~.,scales = 'free_y')+
  ylab('ACS 2010')+
  xlab('Census 2010')+
  ylim(0,2000)+xlim(0,2000)

ggplot(nn,aes(x=acs_pov_pct,y=ce_pop-acs_pop))+
  geom_point()+
  geom_smooth(method='lm')+
  facet_wrap(agecat~.,scales = 'free_y')+
  ylab('Census-ACS')+
  xlab('ACS Percent Poverty')
dev.off()