##########################################################################################
## test for significant time trends in ACS population estimates for overlapping periods ##
##########################################################################################

library(data.table)
library(ggplot2)
library(GGally)
library(cowplot)
library(reshape2)
library(plyr)
library(RColorBrewer)
library(tidycensus)

ma_acs_08<- get_acs(geography = "tract",
                 variables = c('B01001H_001','B01001B_001'),
                 state = "MA",year=2010)

ma_acs_09<- get_acs(geography = "tract",
                    variables = c('B01001H_001','B01001B_001'),
                    state = "MA",year=2011)

ma_acs_10<- get_acs(geography = "tract",
                    variables = c('B01001H_001','B01001B_001'),
                    state = "MA",year=2012)

## compare over time ##
comp1<-merge(ma_acs_08,ma_acs_09,by=c('GEOID','variable'),suffixes=c('_08','_09'))
setDT(comp1)
comp1[,':='(se_08=moe_08/1.645,se_09=moe_09/1.645)]
comp1[variable=='B01001B_001',mean(abs((estimate_08-estimate_09)/(sqrt(1-(4/5))*sqrt(se_08^2+se_09^2)))>1.96)]
comp1[variable=='B01001H_001',mean(abs((estimate_08-estimate_09)/(sqrt(1-(4/5))*sqrt(se_08^2+se_09^2)))>1.96)]

comp2<-merge(ma_acs_09,ma_acs_10,by=c('GEOID','variable'),suffixes=c('_09','_10'))
setDT(comp2)
comp2[,':='(se_09=moe_09/1.645,se_10=moe_10/1.645)]
comp2[variable=='B01001B_001',mean(abs((estimate_09-estimate_10)/(sqrt(1-(4/5))*sqrt(se_09^2+se_10^2)))>1.96)]
comp2[variable=='B01001H_001',mean(abs((estimate_09-estimate_10)/(sqrt(1-(4/5))*sqrt(se_09^2+se_10^2)))>1.96)]

comp3<-merge(ma_acs_08,ma_acs_10,by=c('GEOID','variable'),suffixes=c('_08','_10'))
setDT(comp3)
comp3[,':='(se_08=moe_08/1.645,se_10=moe_10/1.645)]
comp3[variable=='B01001B_001',mean(abs((estimate_08-estimate_10)/(sqrt(1-(4/5))*sqrt(se_08^2+se_10^2)))>1.96)]
comp3[variable=='B01001H_001',mean(abs((estimate_08-estimate_10)/(sqrt(1-(4/5))*sqrt(se_08^2+se_10^2)))>1.96)]