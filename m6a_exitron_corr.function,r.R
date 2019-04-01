##########################2019-04-01###############################
##Aim: test the corrleation between m6a and exitron################
##Result: from the scatterplot,there is no obvious correlation#####
###################################################################
library(dplyr)
library(tidyr)
library(ggplot2)


##function for correlation between m6a and exitron##
m6a_exitron<-function(exitron1,exitron2,m6a3){
##read exitron data##
nilo1<-read.delim(exitron1,header=FALSE)
nilo2<-read.delim(exitron2,header=FALSE)

##subset df (group on genes and choose maximum pso value for each gene)##
nilo1gene<-nilo1 %>% 
  group_by(V7) %>% 
  filter(V11==max(V11)) %>%   ##filter row
  select(one_of(c("V7","V11"))) ##select column

nilo2gene<-nilo2 %>% 
  group_by(V7) %>% 
  filter(V11==max(V11))%>%
  select(one_of(c("V7","V11")))
nilogene<-full_join(nilo1gene,nilo2gene,by="V7")%>%rename(nilo1=V11.x,nilo2=V11.y)
nilogene[is.na(nilogene)]<-0
nilogene$exitron<-apply(nilogene[,2:3],1,mean)
nilogene<-rename(nilogene,name=V7)

##read m6a data##
##methylation info from nilo and parental##
nilom6a<-read.delim(m6a3,sep=",",skip=1,header =TRUE)
nilom6a<-nilom6a[,4:5]
nilom6a<-rename(nilom6a,m6a=score) %>% 
  group_by(name) %>% 
  filter(m6a==max(m6a))%>%
  distinct()

##full_join exitron and m6a data##
nilo_corr<-full_join(nilogene,nilom6a,by="name")
nilo_corr[is.na(nilo_corr)]<-0
nilo_corr$group<-rep(strsplit(m6a3,"[.]")[[1]][1],nrow(nilo_corr))
return(nilo_corr)
}

parental_corr<-m6a_exitron("parental1_IN_sorted.junction.exitron","parental2_IN_sorted.junction.exitron","Parental.m6a.csv")
nilo_corr<-m6a_exitron("nilo1_IN_sorted.junction.exitron","nilo2_IN_sorted.junction.exitron","resistant.m6a.csv")

##all corr data##
corrplot<-rbind(nilo_corr[,c(1,4:6)],parental_corr[,c(1,4:6)])
##plot##
sp <- ggplot(corrplot, aes(x=m6a, y=exitron)) + geom_point(shape=1) +facet_grid(. ~ group)
sp

