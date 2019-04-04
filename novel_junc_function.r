####################2019-04-03#########
##TCGA junction gz file################
##part intron retained junctions stat##
#######################################
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(magrittr)

novel_junc<-function(input){

df = read.table(gzfile(input),header = TRUE)

##canonical splicing site##
splicesite<-c("GT-AG","GC-AG","AT-AC")
df<-df[df$score>=3&toupper(df$splice_site)%in%splicesite,]
##subset by strand##
df_pos<-df[df$strand=="+",]
df_neg<-df[df$strand=="-",]

##positive strand (+)##
##DA,D and A df##
df_pos.DA<-df_pos[df_pos$anchor=="DA",]
df_pos.D<-df_pos[df_pos$anchor=="D",]
df_pos.A<-df_pos[df_pos$anchor=="A",]

df_pos.DAwithD.A<-df_pos.DA[(df_pos.DA$start%in%df_pos.D$start)&(df_pos.DA$end%in%df_pos.A$end),]
df_pos.D.withDA<-df_pos.D[df_pos.D$start%in%df_pos.DAwithD.A$start,]
df_pos.A.withDA<-df_pos.A[df_pos.A$end%in%df_pos.DAwithD.A$end,]

##combine above three together##
df_pos.D.DA.A<-rbind(df_pos.D.withDA,df_pos.DAwithD.A,df_pos.A.withDA)%>%
               group_by(genes)%>%
               arrange(chrom,start,end)


##choose genes with D,A and DA at least once##
df_pos.D.DA.A$genes<-as.character(df_pos.D.DA.A$genes)
df_pos.D.DA.A$start<-as.integer(as.character(df_pos.D.DA.A$start))
df_pos.D.DA.A$end<-as.integer(as.character(df_pos.D.DA.A$end))
genepos.sub<-NULL
for (i in 1:length(unique(df_pos.D.DA.A$genes))){
  genepos<-unique(df_pos.D.DA.A$genes)
  
  a<-nrow(df_pos.D.DA.A[df_pos.D.DA.A$genes==genepos[i]&df_pos.D.DA.A$anchor=="A",])>=1
  b<-nrow(df_pos.D.DA.A[df_pos.D.DA.A$genes==genepos[i]&df_pos.D.DA.A$anchor=="D",])>=1
  c<-nrow(df_pos.D.DA.A[df_pos.D.DA.A$genes==genepos[i]&df_pos.D.DA.A$anchor=="DA",])>=1
  
  if (a&b&c)genepos.sub[i]<-na.omit(genepos[i])
  
}

df_pos.D.DA.A<-df_pos.D.DA.A[df_pos.D.DA.A$genes%in%genepos.sub,]

##D end < A start##
genepos.sub<-NULL
for (i in 1:length(unique(df_pos.D.DA.A$genes))){
genepos<-unique(df_pos.D.DA.A$genes)

DAs<-unlist(df_pos.D.DA.A[df_pos.D.DA.A$genes==genepos[i]&df_pos.D.DA.A$anchor=="DA","start"])
DAe<-unlist(df_pos.D.DA.A[df_pos.D.DA.A$genes==genepos[i]&df_pos.D.DA.A$anchor=="DA","end"])
Ds<-unlist(df_pos.D.DA.A[df_pos.D.DA.A$genes==genepos[i]&df_pos.D.DA.A$anchor=="D","start"])
De<-unlist(df_pos.D.DA.A[df_pos.D.DA.A$genes==genepos[i]&df_pos.D.DA.A$anchor=="D","end"])
As<-unlist(df_pos.D.DA.A[df_pos.D.DA.A$genes==genepos[i]&df_pos.D.DA.A$anchor=="A","start"])
Ae<-unlist(df_pos.D.DA.A[df_pos.D.DA.A$genes==genepos[i]&df_pos.D.DA.A$anchor=="A","end"])

if (De<max(As))genepos.sub[i]<-genepos[i]
}
genepos.sub<-na.omit(genepos.sub)

##negtive strand(-)##
##DA,D and A df##
df_neg.DA<-df_neg[df_neg$anchor=="DA",]
df_neg.D<-df_neg[df_neg$anchor=="D",]
df_neg.A<-df_neg[df_neg$anchor=="A",]

df_neg.DAwithD.A<-df_neg.DA[(df_neg.DA$start%in%df_neg.A$start)&(df_neg.DA$end%in%df_neg.D$end),]
df_neg.D.withDA<-df_neg.D[df_neg.D$end%in%df_neg.DAwithD.A$end,]
df_neg.A.withDA<-df_neg.A[df_neg.A$start%in%df_neg.DAwithD.A$start,]

##combine above three together##
df_neg.D.DA.A<-rbind(df_neg.A.withDA,df_neg.DAwithD.A,df_neg.D.withDA)%>%
  group_by(genes)%>%
  arrange(chrom,start,end)


##choose genes with D,A and DA at least once##
df_neg.D.DA.A$genes<-as.character(df_neg.D.DA.A$genes)
df_neg.D.DA.A$start<-as.integer(as.character(df_neg.D.DA.A$start))
df_neg.D.DA.A$end<-as.integer(as.character(df_neg.D.DA.A$end))
geneneg.sub<-NULL
for (i in 1:length(unique(df_neg.D.DA.A$genes))){
  geneneg<-unique(df_neg.D.DA.A$genes)
  
  a<-nrow(df_neg.D.DA.A[df_neg.D.DA.A$genes==geneneg[i]&df_neg.D.DA.A$anchor=="A",])>=1
  b<-nrow(df_neg.D.DA.A[df_neg.D.DA.A$genes==geneneg[i]&df_neg.D.DA.A$anchor=="D",])>=1
  c<-nrow(df_neg.D.DA.A[df_neg.D.DA.A$genes==geneneg[i]&df_neg.D.DA.A$anchor=="DA",])>=1
  
  if (a&b&c)geneneg.sub[i]<-na.omit(geneneg[i])
  
}

df_neg.D.DA.A<-df_neg.D.DA.A[df_neg.D.DA.A$genes%in%geneneg.sub,]

##D end < A start##
geneneg.sub<-NULL
for (i in 1:length(unique(df_neg.D.DA.A$genes))){
  geneneg<-unique(df_neg.D.DA.A$genes)
  
  DAs<-unlist(df_neg.D.DA.A[df_neg.D.DA.A$genes==geneneg[i]&df_neg.D.DA.A$anchor=="DA","start"])
  DAe<-unlist(df_neg.D.DA.A[df_neg.D.DA.A$genes==geneneg[i]&df_neg.D.DA.A$anchor=="DA","end"])
  Ds<-unlist(df_neg.D.DA.A[df_neg.D.DA.A$genes==geneneg[i]&df_neg.D.DA.A$anchor=="D","start"])
  De<-unlist(df_neg.D.DA.A[df_neg.D.DA.A$genes==geneneg[i]&df_neg.D.DA.A$anchor=="D","end"])
  As<-unlist(df_neg.D.DA.A[df_neg.D.DA.A$genes==geneneg[i]&df_neg.D.DA.A$anchor=="A","start"])
  Ae<-unlist(df_neg.D.DA.A[df_neg.D.DA.A$genes==geneneg[i]&df_neg.D.DA.A$anchor=="A","end"])
  
  if (Ae<max(Ds))geneneg.sub[i]<-geneneg[i]
}
geneneg.sub<-na.omit(geneneg.sub)

##genes that have intron retain splicing##
genecomb<-c(genepos.sub,geneneg.sub)

genejunc<-data_frame(sample=input,number=length(genecomb))

return(genejunc)
}

test<-novel_junc("fab484b9-2ceb-42f8-a8d4-4308c674966e.janno.gz")




