###############2019-04-05###################
##without for loop to select junctions######
############################################
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(magrittr)


df = read.table(gzfile("fab484b9-2ceb-42f8-a8d4-4308c674966e.janno.gz"),header = TRUE)

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

##junctions by DA##  
##D end < A start##
df_pos.D.DA.A.subDA<-df_pos.D.DA.A[df_pos.D.DA.A$anchor=="DA",1:4]%>%
  rename(chromDA=chrom,startDA=start,endDA=end,nameDA=name)

df_pos.D.DA.A.subD<-df_pos.D.DA.A[df_pos.D.DA.A$anchor=="D",1:4]%>%
  rename(chromD=chrom,startD=start,endD=end,nameD=name)

df_pos.D.DA.A.subA<-df_pos.D.DA.A[df_pos.D.DA.A$anchor=="A",1:4]%>%
  rename(chromA=chrom,startA=start,endA=end,nameA=name)

pos_DA.D<-left_join(df_pos.D.DA.A.subDA, df_pos.D.DA.A.subD, by = c("chromDA" = "chromD", "startDA" = "startD"))

pos_DA.D.A<-left_join(pos_DA.D,df_pos.D.DA.A.subA,by=c("chromDA" = "chromA", "endDA" = "endA"))


junctions<-NULL
for (i in 1:nrow(pos_DA.D.A)){
  if (pos_DA.D.A[i,"endD"]<pos_DA.D.A[i,"startA"]) junctions[i]<-i
}
pos_DA.D.A.sub<-na.omit(pos_DA.D.A[junctions,])






