################################################
## bias plots for race stratified simulations ##
################################################

yrs<-2008:2010
fls<-5

proc_sims<-function(simtype){

  simout<-NULL
  for (i in 1:length(yrs)){
    for (j in 1:fls){
      load(paste0('results/',simtype,'/sim_rs_',yrs[i],'_',j,'.RData'))
      simout<-rbind(simout,data.frame('simnm'=yrs[i],'pop_est'=rep(c('ACS','WorldPop'),each=nrow(car1)),rbind(car1,car2)))
    }
  }
  
  betaonly<-simout[,1:(ncov+2)]
  
  ggplotdat<-melt(data=betaonly,id.vars=c('simnm','pop_est'))
  
  betatrue<-data.frame('variable'=unique(ggplotdat$variable),
                       #'truth'=c(beta,(exp(0+0.1)-exp(0))*1000),
                       'truth'=beta,
                       'variablenew'=c('beta[0]: Intercept','beta[1]: I(Black)','beta[2]: PropPov'))
  
  ggplotdat<-merge(ggplotdat,betatrue,by='variable')
  return(ggplotdat)
}

## process crude results ##
dat_crude<-proc_sims('crude')

## process standardized results ##
dat_stand<-proc_sims('standardized')

save(dat_crude,dat_stand,file='results/compiled_rs.RData')

bp_crude <- ggplot(dat_crude, aes(x=factor(simnm), y=value, fill=pop_est)) +
  geom_boxplot()+
  geom_hline(data=dat_crude,aes(yintercept=truth))+
  labs(y = "Estimate")+
  facet_wrap(~ variablenew, ncol=3,scales='free_y',label='label_parsed')+
  scale_fill_brewer(palette="Paired") +
  theme(axis.title.x=element_blank(),text=element_text(size=16,colour='black'),legend.title =element_blank())

bp_stand <- ggplot(dat_stand, aes(x=factor(simnm), y=value, fill=pop_est)) +
  geom_boxplot()+
  geom_hline(data=dat_stand,aes(yintercept=truth))+
  labs(y = "Estimate")+
  facet_wrap(~ variablenew, ncol=3,scales='free_y',label='label_parsed')+
  scale_fill_brewer(palette="Paired") +
  theme(axis.title.x=element_blank(),text=element_text(size=16,colour='black'),legend.title =element_blank())

pdf('bias_plots_rs.pdf',width=10,height=6)
print(bp_crude)
print(bp_stand)
dev.off()