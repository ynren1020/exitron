########2019-01-25########
##nilo parental release###
##exitron statistics######

##packages##
library(dplyr)
library(tidyr)

##read data##
nilo1<-read.delim("nilo1_IN_sorted.junction.exitron",header=FALSE)
nilo2<-read.delim("nilo2_IN_sorted.junction.exitron",header=FALSE)
nilo<-union(nilo1,nilo2)

unionfunc<-function(file1,file2){
  temp1<-read.delim(file1,header=FALSE)
  temp2<-read.delim(file2,header=FALSE)
  temp<-dplyr::union(temp1,temp2)
  return(temp)
  }
  
parental<-unionfunc("parental1_IN_sorted.junction.exitron","parental2_IN_sorted.junction.exitron")
release<-unionfunc("release1_IN_sorted.junction.exitron","release2_IN_sorted.junction.exitron")
exitron<-dplyr::data_frame(nilo=nrow(nilo),parental=nrow(parental),release=nrow(release))
write.csv(exitron,"k562.exitron.csv",sep=",",col.names = TRUE,row.names = FALSE,quote=FALSE)

##select four columns,chr,start,end,strand, and check union for nilo,parental,and release##
##V11 is pso,V14 is total spliced read, each sample have the same value##
unionfuncsub<-function(file1,file2){
  temp1<-read.delim(file1,header=FALSE)
  temp1<-temp1[,c(1,2,3,6)]
  temp2<-read.delim(file2,header=FALSE)
  temp2<-temp2[,c(1,2,3,6)]
  temp<-dplyr::union(temp1,temp2)
  return(temp)
}

nilosub<-unionfuncsub("nilo1_IN_sorted.junction.exitron","nilo2_IN_sorted.junction.exitron")
parentalsub<-unionfuncsub("parental1_IN_sorted.junction.exitron","parental2_IN_sorted.junction.exitron")
releasesub<-unionfuncsub("release1_IN_sorted.junction.exitron","release2_IN_sorted.junction.exitron")

exitronsub<-dplyr::data_frame(nilo=nrow(nilosub),parental=nrow(parentalsub),release=nrow(releasesub))
write.csv(exitronsub,"k562.exitronsub.csv",sep=",",col.names = TRUE,row.names = FALSE,quote=FALSE)


##test if each gene's pso changed signicantly in nilo/parental##
##two replicates for each sample##
##maximum pso as the gene's pso##
parental1<-read.delim("parental1_IN_sorted.junction.exitron",header=FALSE)
parental2<-read.delim("parental2_IN_sorted.junction.exitron",header=FALSE)
release1<-read.delim("release1_IN_sorted.junction.exitron",header=FALSE)
release2<-read.delim("release2_IN_sorted.junction.exitron",header=FALSE)

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


parental1gene<-parental1 %>% 
  group_by(V7) %>% 
  filter(V11==max(V11))%>%
  select(one_of(c("V7","V11")))

parental2gene<-parental2 %>% 
  group_by(V7) %>% 
  filter(V11==max(V11))%>%
  select(one_of(c("V7","V11")))
parentalgene<-full_join(parental1gene,parental2gene,by="V7")%>%rename(parental1=V11.x,parental2=V11.y)

niloparentalgene<-full_join(nilogene,parentalgene,by="V7")
niloparentalgene<-as.data.frame(niloparentalgene)
niloparentalgene$nilo1<-as.numeric(as.character(niloparentalgene$nilo1))
niloparentalgene$nilo2<-as.numeric(as.character(niloparentalgene$nilo2))
niloparentalgene$parental1<-as.numeric(as.character(niloparentalgene$parental1))
niloparentalgene$parental2<-as.numeric(as.character(niloparentalgene$parental2))
##remove NA observations##
niloparentalgenesub<-na.omit(niloparentalgene)

niloparentalgenesub$pvals <- rep(NA, nrow(niloparentalgenesub))
for(i in 1:nrow(niloparentalgenesub)) niloparentalgenesub$pvals[i] <- t.test(niloparentalgenesub[i,c(2,3)],niloparentalgenesub[i,c(4,5)])$p.value
niloparentalgenesub$pvals.adj<-round(p.adjust(niloparentalgenesub$pvals, "BH"), 3)

##replace NA with 0##
niloparentalgene[is.na(niloparentalgene)]<-0
niloparentalgene$pvals <- rep(NA, nrow(niloparentalgene))
for(i in 1:nrow(niloparentalgene)) niloparentalgene$pvals[i] <- t.test(niloparentalgene[i,c(2,3)],niloparentalgene[i,c(4,5)])$p.value
niloparentalgene$pvals.adj<-round(p.adjust(niloparentalgene$pvals, "BH"), 3)

##linear regression nilo and parental pso data##
group<-factor(c(1,1,2,2))
design<-model.matrix(~0+group)
#fit<-NULL
pval<-NULL
for (i in 1:nrow(niloparentalgene)){
#fit[i]<-lm(t(as.matrix(niloparentalgene[i,2:5]))~design)
pval[i]<-anova(lm(t(as.matrix(niloparentalgene[i,2:5]))~design))$'Pr(>F)'[1]
}
niloparentalgene$pval.lm<-pval
niloparentalgene.sig<-niloparentalgene[niloparentalgene$pval.lm<0.1,]  ##20genes differential pso between nilo and parental
niloparentalgene.sig.nilo<-niloparentalgene.sig[niloparentalgene.sig$nilo1!=0|niloparentalgene.sig$nilo2!=0,]
niloparentalgene.sig.nilo<-rename(niloparentalgene.sig.nilo,name=V7)
niloparentalgene.sig.nilonoDEG<-left_join(niloparentalgene.sig.nilo,nilogeneidDEG,by="name")

niloparentalgene.sig.nilonoDEGsub<-niloparentalgene.sig.nilonoDEG[c(1,2,9),]
write.csv(niloparentalgene.sig.nilonoDEGsub,"niloparentalgene.sig.nilonoDEGsub.csv",sep=",",col.names = TRUE,row.names = FALSE,quote=FALSE)

##pso diff,no DEG, m6a#######FINAL STOP SITE,MAY USE THIS ONE####
niloparentalgene_all<-rename(niloparentalgene,name=V7)
niloparentalgene_all<-distinct(niloparentalgene_all)
niloparentalgene_all<-left_join(niloparentalgene_all,nilogeneidDEG,by="name")
niloparentalgene_all<-left_join(niloparentalgene_all,nilom6a,by="name")
niloparentalgene_all<-left_join(niloparentalgene_all,parentalm6a,by="name")
niloparentalgene_all<-niloparentalgene_all[!is.na(niloparentalgene_all$Length),]
#niloparentalgene_allsub<-niloparentalgene_all[(!is.na(niloparentalgene_all$nilom6a))&(!is.na(niloparentalgene_all$parentalm6a)),]
niloparentalgene_all$nilo<-(niloparentalgene_all$nilo1+niloparentalgene_all$nilo2)/2
niloparentalgene_all$parental<-(niloparentalgene_all$parental1+niloparentalgene_all$parental2)/2
niloparentalgene_all[is.na(niloparentalgene_all)]<-0
##nilo > parental and nilom6a > parentalm6a##
niloparentalgene_all_large<-filter(niloparentalgene_all,(nilo>parental)&(nilom6a>parentalm6a))
write.csv(niloparentalgene_all_large,"niloparentalgene_all_large.csv",sep=",",col.names = TRUE,row.names = FALSE,quote=FALSE)

##################

##nilo exitron##
niloparentalgene_nilo2<-niloparentalgene[niloparentalgene$nilo1!=0|niloparentalgene$nilo2!=0,]
##parental exitron##
niloparentalgene_parental<-niloparentalgene[niloparentalgene$nilo1==0&niloparentalgene$nilo2==0,]


##subset##
##3NAs##
keep<-rowSums(niloparentalgene[,2:5]==0)>2
niloparentalgene3NAs<-niloparentalgene[keep,]

##less than or equal 2NAs##
keep2<-rowSums(niloparentalgene[,2:5]==0)<=2
niloparentalgenesubs<-niloparentalgene[keep2,]

##methylation info from nilo and parental##
nilom6a<-read.delim("resistant.m6a.csv",sep=",",skip=1,header =TRUE)
nilom6a<-nilom6a[,4:5]
nilom6a<-rename(nilom6a,nilom6a=score) %>% 
  group_by(name) %>% 
  filter(nilom6a==max(nilom6a))%>%
  distinct()
  

parentalm6a<-read.delim("Parental.m6a.csv",sep=",",skip=1,header =TRUE)
parentalm6a<-parentalm6a[,4:5]
parentalm6a<-rename(parentalm6a,parentalm6a=score)%>% 
  group_by(name) %>% 
  filter(parentalm6a==max(parentalm6a))%>%
  distinct()

##nilo exitron;nilo m6a##
niloparentalgene_nilo_m6a<-niloparentalgene_nilo[niloparentalgene_nilo$V7%in%nilom6a$name,]
niloparentalgene_nilo_m6a<-rename(niloparentalgene_nilo_m6a,name=V7)
niloparentalgene_nilo_m6a<-left_join(niloparentalgene_nilo_m6a,nilom6a)
niloparentalgene_nilo_m6a<-left_join(niloparentalgene_nilo_m6a,parentalm6a)

##only exitron and m6a in nilo##
niloparentalgene_nilo_m6a_onlynilo<-niloparentalgene_nilo_m6a[is.na(niloparentalgene_nilo_m6a$parentalm6a),]  #15
##"TCIRG1" "TESPA1" "ZZEF1"  "ITGA2B" "KCNC3"  "HOXD8"  "CCRL2"  "HIVEP2" "ZER1"   "MLLT10" "CD44"   "DICER1" "BZRAP1" "POTEF"  "KANK1" 
##only exitron in nilo,m6a score higher in nilo##
niloparentalgene_nilo_m6a_highnilo<-niloparentalgene_nilo_m6a[!is.na(niloparentalgene_nilo_m6a$parentalm6a),]
niloparentalgene_nilo_m6a_highnilo$ratiom6a<-niloparentalgene_nilo_m6a_highnilo$nilom6a/niloparentalgene_nilo_m6a_highnilo$parentalm6a #59
niloparentalgene_nilo_m6a_highnilo<-niloparentalgene_nilo_m6a_highnilo[niloparentalgene_nilo_m6a_highnilo$ratiom6a>1,] #39

write.csv(niloparentalgene_nilo_m6a_onlynilo,"niloparentalgene_nilo_m6a_onlynilo.csv",sep=",",col.names = TRUE,row.names = FALSE)
write.csv(niloparentalgene_nilo_m6a_highnilo,"niloparentalgene_nilo_m6a_highnilo.csv",sep=",",col.names = TRUE,row.names = FALSE)

##RNA-seq differential expression genes parental.degs.txt##
parentalDEG<-read.delim("parental.degs.txt")
parentalDEGsub<-parentalDEG[parentalDEG$adj.P.Val>0.05,]
parentalDEGsub$V10<-row.names(parentalDEGsub)

##nilo gene ensemble id##
nilogeneid<-nilo[,c(7,10)]%>%
  distinct()

nilogeneidDEG<-left_join(nilogeneid,parentalDEGsub,by="V10")%>%
  rename(name=V7)

##exitron in nilo,m6a score higher or only in nilo##
niloexitron_m6a<-rbind(niloparentalgene_nilo_m6a_onlynilo,niloparentalgene_nilo_m6a_highnilo[,-10])
niloexitron_m6anoDEG<-left_join(niloexitron_m6a,nilogeneidDEG,by="name")
niloexitron_m6anoDEG_sub<-na.omit(niloexitron_m6anoDEG)

