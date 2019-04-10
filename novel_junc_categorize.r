###############2019-04-10#########################
##categorize the novel junctions to two###########
##First: no N; second:with N between D and A######
#"fab484b9-2ceb-42f8-a8d4-4308c674966e.janno.gz"##
##################################################

#library(readr)
library(dplyr)
#library(ggplot2)
library(tidyr)
#library(magrittr)


args <- commandArgs(TRUE)
  
df = read.table(gzfile(args[1]),header = TRUE)

##canonical splicing site##
splicesite<-c("GT-AG","GC-AG","AT-AC")
df<-df[df$score>=3&toupper(df$splice_site)%in%splicesite,]
df$chrom<-as.character(df$chrom)
##subset by strand##
df_pos<-df[df$strand=="+",]
df_neg<-df[df$strand=="-",]
##positive strand (+)##
##DA,D, N and A df## NOTE: N between D and A##
  df_pos.DA<-df_pos[df_pos$anchor=="DA",c(1:4,11)]%>%
    rename(chromDA=chrom,startDA=start,endDA=end,nameDA=name,anchorDA=anchor)
  df_pos.D<-df_pos[df_pos$anchor=="D",c(1:4,11)]%>%
    rename(chromD=chrom,startD=start,endD=end,nameD=name,anchorD=anchor)
  df_pos.A<-df_pos[df_pos$anchor=="A",c(1:4,11)]%>%
    rename(chromA=chrom,startA=start,endA=end,nameA=name,anchorA=anchor)
  df_pos.N<-df_pos[df_pos$anchor=="N",c(1:4,11)]%>%
    rename(chromN=chrom,startN=start,endN=end,nameN=name,anchorN=anchor)
  ##match Dstart with DAstart and then match Aend with DAend##
  pos_DA.D<-left_join(df_pos.DA, df_pos.D, by = c("chromDA" = "chromD", "startDA" = "startD"))
  pos_DA.D.A<-left_join(pos_DA.D,df_pos.A,by=c("chromDA" = "chromA", "endDA" = "endA"))
  pos.final<-pos_DA.D.A%>% drop_na()
  
  ##Dend < Astart##
  junctions<-NULL
  for (i in 1:nrow(pos.final)){
    if (pos.final[i,"endD"]<pos.final[i,"startA"]) junctions[i]<-i
  }
  ##positive strand,junction combinations##
  pos.final.sub<-pos.final[junctions,]%>% drop_na()
  
##Dend <Nstart<Nend<Astart#############
##DA,D and A unique combination number##
  junctionDA<-NULL
  junctionN<-NULL
  for (i in 1:nrow(pos.final.sub)){
    for (j in 1:nrow(df_pos.N)){
      ifelse(pos.final.sub$chromDA[i]==df_pos.N$chromN[j]&pos.final.sub$endD[i]<df_pos.N$startN[j]&pos.final.sub$startA[i]>df_pos.N$endN[j],junctionDA[i]<-i,NA)
      ifelse(pos.final.sub$chromDA[i]==df_pos.N$chromN[j]&pos.final.sub$endD[i]<df_pos.N$startN[j]&pos.final.sub$startA[i]>df_pos.N$endN[j],junctionN[j]<-j,NA)
      
    }
  }

pos.final.sub.N<-na.omit(pos.final.sub[junctionDA,]) #21 (use this number)
df_pos.N.DA<-na.omit(df_pos.N[junctionN,]) #13

pos.final.sub.N.all<-left_join(pos.final.sub.N,df_pos.N.DA,by = c("chromDA" = "chromN")) #152 (not this number)

junctionDAN<-NULL
for (i in 1:nrow(pos.final.sub.N.all)){
  junctionDAN[i]<-ifelse(pos.final.sub.N.all$endD[i]<pos.final.sub.N.all$startN[i]&pos.final.sub.N.all$startA[i]>pos.final.sub.N.all$endN[i],i,NA)
}
junctionDAN<-na.omit(junctionDAN)
pos.final.sub.N.all<-pos.final.sub.N.all[junctionDAN,] #63 number with N
  
##################################negative strand (-)####################################
  ##DA,D and A df##
  df_neg.DA<-df_neg[df_neg$anchor=="DA",c(1:4,11)]%>%
    rename(chromDA=chrom,startDA=start,endDA=end,nameDA=name,anchorDA=anchor)
  df_neg.D<-df_neg[df_neg$anchor=="D",c(1:4,11)]%>%
    rename(chromD=chrom,startD=start,endD=end,nameD=name,anchorD=anchor)
  df_neg.A<-df_neg[df_neg$anchor=="A",c(1:4,11)]%>%
    rename(chromA=chrom,startA=start,endA=end,nameA=name,anchorA=anchor)
  df_neg.N<-df_neg[df_neg$anchor=="N",c(1:4,11)]%>%
    rename(chromN=chrom,startN=start,endN=end,nameN=name,anchorN=anchor)
  ##match Dstart with DAstart and then match Aend with DAend##
  neg_DA.A<-left_join(df_neg.DA, df_neg.A, by = c("chromDA" = "chromA", "startDA" = "startA"))
  neg_DA.D.A<-left_join(neg_DA.A,df_neg.D,by=c("chromDA" = "chromD", "endDA" = "endD"))
  neg.final<-neg_DA.D.A%>% drop_na()
  
  ##Dend < Astart##
  junctions<-NULL
  for (i in 1:nrow(neg.final)){
    if (neg.final[i,"endA"]<neg.final[i,"startD"]) junctions[i]<-i
  }
  ##negative strand,junction combinations##
  neg.final.sub<-neg.final[junctions,]%>% drop_na()
  
  ##Aend <Nstart<Nend<Dstart#############
  ##DA,D and A unique combination number##
  junctionDA<-NULL
  junctionN<-NULL
  for (i in 1:nrow(neg.final.sub)){
    for (j in 1:nrow(df_neg.N)){
      ifelse(neg.final.sub$chromDA[i]==df_neg.N$chromN[j]&neg.final.sub$endA[i]<df_neg.N$startN[j]&neg.final.sub$startD[i]>df_neg.N$endN[j],junctionDA[i]<-i,NA)
      ifelse(neg.final.sub$chromDA[i]==df_neg.N$chromN[j]&neg.final.sub$endA[i]<df_neg.N$startN[j]&neg.final.sub$startD[i]>df_neg.N$endN[j],junctionN[j]<-j,NA)
      
    }
  }
  
  neg.final.sub.N<-na.omit(neg.final.sub[junctionDA,]) #3 (use this number)
  df_neg.N.DA<-na.omit(df_neg.N[junctionN,]) #3
  
  neg.final.sub.N.all<-left_join(neg.final.sub.N,df_neg.N.DA,by = c("chromDA" = "chromN")) #152 (not this number)
  
  junctionDAN<-NULL
  for (i in 1:nrow(neg.final.sub.N.all)){
    junctionDAN[i]<-ifelse(neg.final.sub.N.all$endA[i]<neg.final.sub.N.all$startN[i]&neg.final.sub.N.all$startD[i]>neg.final.sub.N.all$endN[i],i,NA)
  }
  junctionDAN<-na.omit(junctionDAN)
  neg.final.sub.N.all<-neg.final.sub.N.all[junctionDAN,] #4 number with N
  
  
##all novel junction##  
junc.number<-nrow(pos.final.sub)+nrow(neg.final.sub) #289
##junction with N between D and A##
junc.N.number<-nrow(pos.final.sub.N)+nrow(neg.final.sub.N) #24
##output##
junc.df<-dplyr::tibble(sample=strsplit(args[1],"/")[[1]][9],novel=junc.number,novel.N=junc.N.number,total=nrow(df))
output<-paste0(strsplit(args[1],"/")[[1]][9],".N.txt")
write.table(junc.df,output,quote=FALSE,col.names=FALSE,row.names=FALSE)




