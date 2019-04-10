###############2019-04-05###################
##without for loop to select junctions######
#"fab484b9-2ceb-42f8-a8d4-4308c674966e.janno.gz"
############################################
#library(readr)
library(dplyr)
#library(ggplot2)
library(tidyr)
#library(magrittr)


novel_juncfinal<-function(input){

df = read.table(gzfile(input),header = TRUE)

##canonical splicing site##
splicesite<-c("GT-AG","GC-AG","AT-AC")
df<-df[df$score>=3&toupper(df$splice_site)%in%splicesite,]
df$chrom<-as.character(df$chrom)
##subset by strand##
df_pos<-df[df$strand=="+",]
df_neg<-df[df$strand=="-",]
##positive strand (+)##
##DA,D and A df##
df_pos.DA<-df_pos[df_pos$anchor=="DA",c(1:4,11)]%>%
  rename(chromDA=chrom,startDA=start,endDA=end,nameDA=name,anchorDA=anchor)
df_pos.D<-df_pos[df_pos$anchor=="D",c(1:4,11)]%>%
  rename(chromD=chrom,startD=start,endD=end,nameD=name,anchorD=anchor)
df_pos.A<-df_pos[df_pos$anchor=="A",c(1:4,11)]%>%
  rename(chromA=chrom,startA=start,endA=end,nameA=name,anchorA=anchor)
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


##negative strand (-)##
##DA,D and A df##
df_neg.DA<-df_neg[df_neg$anchor=="DA",c(1:4,11)]%>%
  rename(chromDA=chrom,startDA=start,endDA=end,nameDA=name,anchorDA=anchor)
df_neg.D<-df_neg[df_neg$anchor=="D",c(1:4,11)]%>%
  rename(chromD=chrom,startD=start,endD=end,nameD=name,anchorD=anchor)
df_neg.A<-df_neg[df_neg$anchor=="A",c(1:4,11)]%>%
  rename(chromA=chrom,startA=start,endA=end,nameA=name,anchorA=anchor)
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

junc.number<-nrow(pos.final.sub)+nrow(neg.final.sub) #289
junc.df<-data_frame(sample=input,novel=junc.number,total=nrow(df))
return(junc.df)

}

#test<-novel_juncfinal("fab484b9-2ceb-42f8-a8d4-4308c674966e.janno.gz")





