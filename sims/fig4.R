############################################################
## create bias boxplots for simulation results (figure 4) ##
############################################################

library(cowplot)
library(ggplot2)

load('results/compiled_total.RData')
total<-dat_stand
load('results/compiled_rs.RData')
rs<-dat_stand

bp_total_stand <- ggplot(total, aes(x=factor(simnm), y=value, fill=pop_est)) +
  geom_boxplot()+
  geom_hline(data=total,aes(yintercept=truth))+
  labs(y = "Estimate")+
  facet_wrap(~ variablenew, ncol=3,scales='free_y',label='label_parsed')+
  scale_fill_brewer(palette="Paired") +
  theme(axis.title.x=element_blank(),text=element_text(size=18,colour='black'),legend.title =element_blank())

bp_rs_stand <- ggplot(rs, aes(x=factor(simnm), y=value, fill=pop_est)) +
  geom_boxplot()+
  geom_hline(data=rs,aes(yintercept=truth))+
  labs(y = "Estimate")+
  facet_wrap(~ variablenew, ncol=3,scales='free_y',label='label_parsed')+
  scale_fill_brewer(palette="Paired") +
  theme(axis.title.x=element_blank(),text=element_text(size=18,colour='black'),legend.title =element_blank())



load('results/compiled_total.RData')
total<-dat_crude
load('results/compiled_rs.RData')
rs<-dat_crude

bp_total_crude <- ggplot(total, aes(x=factor(simnm), y=value, fill=pop_est)) +
  geom_boxplot()+
  geom_hline(data=total,aes(yintercept=truth))+
  labs(y = "Estimate")+
  facet_wrap(~ variablenew, ncol=3,scales='free_y',label='label_parsed')+
  scale_fill_brewer(palette="Paired") +
  theme(axis.title.x=element_blank(),text=element_text(size=18,colour='black'),legend.title =element_blank())

bp_rs_crude <- ggplot(rs, aes(x=factor(simnm), y=value, fill=pop_est)) +
  geom_boxplot()+
  geom_hline(data=rs,aes(yintercept=truth))+
  labs(y = "Estimate")+
  facet_wrap(~ variablenew, ncol=3,scales='free_y',label='label_parsed')+
  scale_fill_brewer(palette="Paired") +
  theme(axis.title.x=element_blank(),text=element_text(size=18,colour='black'),legend.title =element_blank())

pdf('bias_plots.pdf',width=12,height=8)
plot_grid(bp_total_stand, bp_rs_stand, labels = c('A', 'B'),ncol=1)
plot_grid(bp_total_crude, bp_rs_crude, labels = c('A', 'B'),ncol=1)
dev.off()
