##############2019-04-04#########################
##find intron retain splicing junction###########
##TCGA junction gz file##########################
##/home/tywang/Projects/Exitron/TCGA/ACC/jannos##
#################################################

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(magrittr)

#novel_junc<-function(input){
  
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
df_pos.D.DA.A.subDA<-df_pos.D.DA.A[df_pos.D.DA.A$anchor=="DA",]
df_pos.D.DA.A.subD<-df_pos.D.DA.A[df_pos.D.DA.A$anchor=="D",]
df_pos.D.DA.A.subA<-df_pos.D.DA.A[df_pos.D.DA.A$anchor=="A",]

junction.pos.DA<-NULL
junction.pos.D<-NULL
junction.pos.A<-NULL
  for (i in 1:nrow(df_pos.D.DA.A.subDA)){
    for (j in 1:nrow(df_pos.D.DA.A.subD)){
      for (k in 1:nrow(df_pos.D.DA.A.subA)){
        DAs<-df_pos.D.DA.A.subDA[i,"start"]
        DAe<-df_pos.D.DA.A.subDA[i,"end"]
        Ds<-df_pos.D.DA.A.subD[j,"start"]
        De<-df_pos.D.DA.A.subD[j,"end"]
        As<-df_pos.D.DA.A.subA[k,"start"]
        Ae<-df_pos.D.DA.A.subA[k,"end"]
        if((DAs==Ds)&(DAe==Ae)&(De<As)) {junction.pos.DA[i]<-df_pos.D.DA.A.subDA[i,"name"];
                                        junction.pos.D[j]<-df_pos.D.DA.A.subD[j,"name"];
                                        junction.pos.A[k]<-df_pos.D.DA.A.subA[k,"name"];
                                      # print(c(i,j,k));
                                        
                                        }
        else junction.pos[i]<-NA
        
    }
  }
}
 
##use output of print(i,j,k),select rows of DA,D and A file##
pos.id<-read.delim2("pos.junc.id.txt") #nrow=211 (DA,D,A unique combination meet the above criteria)
pos.DA.junc.id<-unique(pos.id$i)
pos.D.junc.id<-unique(pos.id$j)
pos.A.junc.id<-unique(pos.id$k)

df_pos.D.DA.A.subDA.id<-df_pos.D.DA.A.subDA[pos.DA.junc.id,]#150
df_pos.D.DA.A.subD.id<-df_pos.D.DA.A.subD[pos.D.junc.id,]
df_pos.D.DA.A.subA.id<-df_pos.D.DA.A.subA[pos.A.junc.id,]

######################negtive strand(-)##############################
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
  
##subset DA,Aand D##
  df_neg.D.DA.A.subDA<-df_neg.D.DA.A[df_neg.D.DA.A$anchor=="DA",]
  df_neg.D.DA.A.subD<-df_neg.D.DA.A[df_neg.D.DA.A$anchor=="D",]
  df_neg.D.DA.A.subA<-df_neg.D.DA.A[df_neg.D.DA.A$anchor=="A",]
  
##D end < A start##
junction.neg.DA<-NULL
junction.neg.D<-NULL
junction.neg.A<-NULL
sink("test.pos.junc.id.txt")
  for (i in 1:nrow(df_neg.D.DA.A.subDA)){
    for (j in 1:nrow(df_neg.D.DA.A.subD)){
      for (k in 1:nrow(df_neg.D.DA.A.subA)){
        DAs<-df_neg.D.DA.A.subDA[i,"start"]
        DAe<-df_neg.D.DA.A.subDA[i,"end"]
        Ds<-df_neg.D.DA.A.subD[j,"start"]
        De<-df_neg.D.DA.A.subD[j,"end"]
        As<-df_neg.D.DA.A.subA[k,"start"]
        Ae<-df_neg.D.DA.A.subA[k,"end"]
        if((DAs==As)&(DAe==De)&(Ae<Ds)) {
        junction.neg.DA[i]<-df_neg.D.DA.A.subDA[i,"name"];
        junction.neg.D[j]<-df_neg.D.DA.A.subD[j,"name"];
        junction.neg.A[k]<-df_neg.D.DA.A.subA[k,"name"];
        print(c(i,j,k))
        }
        #else junction.neg[i]<-NA
        
      }
    }
  }
sink()
  
 # genejunc<-data_frame(sample=input,number=length(genecomb))
  
 # return(genejunc)
#}



