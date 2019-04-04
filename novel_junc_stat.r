####################2019-03-27#########
##TCGA junction gz file################
##part intron retained junctions stat##
#######################################
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(magrittr)

df = read.table(gzfile("fab484b9-2ceb-42f8-a8d4-4308c674966e.janno.gz"),header = TRUE)
df_1<-df%>%
  group_by(chrom)%>%
  arrange(chrom,start,end)
##canonical splicing site##
splicesite<-c("GT-AG","GC-AG","AT-AC")

df_chr1<-df[df$chrom=="chr1",]
df_chr1<-df_chr1[df_chr1$score>=3&toupper(df_chr1$splice_site)%in%splicesite,]

##int (start and end) to factor##
cols <- c("start", "end")
df_chr1 %<>%
  mutate_each_(funs(factor(.)),cols)
str(df_chr1)

##how many duplciated start##
sum(table(df_chr1$start)>=2) #2285
length(unique(df_chr1dup.s$start)) #2285

#village[village$Names %in% names(which(table(village$Names) > 1)), ]

##start and end duplicated files##
df_chr1dup.s<-df_chr1[df_chr1$start %in% names(which(table(df_chr1$start)>1)),]
df_chr1dup.e<-df_chr1[df_chr1$end %in% names(which(table(df_chr1$end)>1)),]
##combined##
df_chr1dup<-rbind(df_chr1dup.s,df_chr1dup.e)%>%
  arrange(start,end,name)



##select splice_site in splicesite rows##
##start file##
df_chr1dup.s.sub<-df_chr1dup.s[toupper(df_chr1dup.s$splice_site)%in%splicesite,]
##count how many different anchor for each start position##
df_chr1dup.s.sub.summary<-df_chr1dup.s.sub%>%group_by(start)%>%count(anchor)
##select rows with start position duplicated##
df_chr1dup.s.sub.summary<-df_chr1dup.s.sub.summary[df_chr1dup.s.sub.summary$start %in% names(which(table(df_chr1dup.s.sub.summary$start)>1)),]

##end file##
df_chr1dup.e.sub<-df_chr1dup.e[toupper(df_chr1dup.e$splice_site)%in%splicesite,]
##count how many different anchor for each start position##
df_chr1dup.e.sub.summary<-df_chr1dup.e.sub%>%group_by(end)%>%count(anchor)
##select rows with start position duplicated##
df_chr1dup.e.sub.summary<-df_chr1dup.e.sub.summary[df_chr1dup.e.sub.summary$end %in% names(which(table(df_chr1dup.e.sub.summary$end)>1)),]

##find junctions have overlap with two summary files##
start.unique<-unique(df_chr1dup.s.sub.summary$start)
end.unique<-unique(df_chr1dup.e.sub.summary$end)

df_chr1dup.junc<-df_chr1dup[df_chr1dup$start%in%start.unique|df_chr1dup$end%in%end.unique,]
df_chr1dup.junc<-unique(df_chr1dup.junc)

##find GENE junctions with DA, A and D##

df_chr1dup.junc.sub<-df_chr1dup.junc[df_chr1dup.junc$genes %in% names(which(table(df_chr1dup.junc$genes)>2)),]
##how many different anchor for each gene##
df_chr1dup.junc.sub.summary<-df_chr1dup.junc.sub%>%group_by(genes)%>%count(anchor)

##keep genes at least 3 anchors##
names(which(table(df_chr1dup.junc.sub.summary$genes)>2))
df_chr1dup.junc.sub.summary.sub<-df_chr1dup.junc.sub.summary[df_chr1dup.junc.sub.summary$genes%in%names(which(table(df_chr1dup.junc.sub.summary$genes)>2)),]

##select genes with three kinds of anchors:D,A,DA##
genesub<-NULL
df_chr1dup.junc.sub.summary.sub$genes<-as.character(df_chr1dup.junc.sub.summary.sub$genes)
for (i in 1:length(unique(df_chr1dup.junc.sub.summary.sub$genes))){
  genes<-unique(df_chr1dup.junc.sub.summary.sub$genes)
a<-nrow(df_chr1dup.junc.sub.summary.sub[df_chr1dup.junc.sub.summary.sub$genes==genes[i]&df_chr1dup.junc.sub.summary.sub$anchor=="A",])>=1
b<-nrow(df_chr1dup.junc.sub.summary.sub[df_chr1dup.junc.sub.summary.sub$genes==genes[i]&df_chr1dup.junc.sub.summary.sub$anchor=="D",])>=1
c<-nrow(df_chr1dup.junc.sub.summary.sub[df_chr1dup.junc.sub.summary.sub$genes==genes[i]&df_chr1dup.junc.sub.summary.sub$anchor=="DA",])>=1
#df_chr1dup.junc.sub.summary.sub[which(a&b&c),]
if(a&b&c) genesub[i]<-genes[i]
genesub<-na.omit(genesub)

}

#View(df_chr1dup.junc.sub.summary.sub[df_chr1dup.junc.sub.summary.sub$genes%in%genesub,])

##find junctions in df_chr1dup.junc.sub overlap with genesub##
##for further isolation of genes with genecomb list###########
df_chr1dup.junc.sub.geneselect<-df_chr1dup.junc.sub[df_chr1dup.junc.sub$genes%in%genesub,]%>%group_by(genes)%>%arrange(start,end)
  
##find genes with junctions: if (+): DAstart==Dstart,DAend==Aend, DAstart<Dend<Astart<DAend:
df_chr1dup.junc.sub.geneselect.pos<-df_chr1dup.junc.sub.geneselect[df_chr1dup.junc.sub.geneselect$strand=="+",]
##pattern D,DA,A####Need to think about it ################
genepos.sub<-NULL
df_chr1dup.junc.sub.geneselect.pos$start<-as.integer(as.character(df_chr1dup.junc.sub.geneselect.pos$start))
df_chr1dup.junc.sub.geneselect.pos$end<-as.integer(as.character(df_chr1dup.junc.sub.geneselect.pos$end))
for (i in 1:length(unique(df_chr1dup.junc.sub.geneselect.pos$genes))){
  genepos<-unique(df_chr1dup.junc.sub.geneselect.pos$genes)
  DAs<-df_chr1dup.junc.sub.geneselect.pos[df_chr1dup.junc.sub.geneselect.pos$genes==genepos[i]&df_chr1dup.junc.sub.geneselect.pos$anchor=="DA","start"]
  DAe<-df_chr1dup.junc.sub.geneselect.pos[df_chr1dup.junc.sub.geneselect.pos$genes==genepos[i]&df_chr1dup.junc.sub.geneselect.pos$anchor=="DA","end"]
  Ds<-df_chr1dup.junc.sub.geneselect.pos[df_chr1dup.junc.sub.geneselect.pos$genes==genepos[i]&df_chr1dup.junc.sub.geneselect.pos$anchor=="D","start"]
  De<-df_chr1dup.junc.sub.geneselect.pos[df_chr1dup.junc.sub.geneselect.pos$genes==genepos[i]&df_chr1dup.junc.sub.geneselect.pos$anchor=="D","end"]
  As<-df_chr1dup.junc.sub.geneselect.pos[df_chr1dup.junc.sub.geneselect.pos$genes==genepos[i]&df_chr1dup.junc.sub.geneselect.pos$anchor=="A","start"]
  Ae<-df_chr1dup.junc.sub.geneselect.pos[df_chr1dup.junc.sub.geneselect.pos$genes==genepos[i]&df_chr1dup.junc.sub.geneselect.pos$anchor=="A","end"]
  
  if(DAs==Ds&DAe==Ae&De<As) genepos.sub[i]<-genepos[i]
}
###################################################
df_chr1dup.junc.sub.geneselect.pos.DA<-df_chr1dup.junc.sub.geneselect.pos[df_chr1dup.junc.sub.geneselect.pos$anchor=="DA",]
df_chr1dup.junc.sub.geneselect.pos.D<-df_chr1dup.junc.sub.geneselect.pos[df_chr1dup.junc.sub.geneselect.pos$anchor=="D",]
df_chr1dup.junc.sub.geneselect.pos.A<-df_chr1dup.junc.sub.geneselect.pos[df_chr1dup.junc.sub.geneselect.pos$anchor=="A",]


df_chr1dup.junc.sub.geneselect.pos.DAwithD.A<-df_chr1dup.junc.sub.geneselect.pos.DA[(df_chr1dup.junc.sub.geneselect.pos.DA$start%in%df_chr1dup.junc.sub.geneselect.pos.D$start)&(df_chr1dup.junc.sub.geneselect.pos.DA$end%in%df_chr1dup.junc.sub.geneselect.pos.A$end),]
df_chr1dup.junc.sub.geneselect.pos.D.withDA<-df_chr1dup.junc.sub.geneselect.pos.D[df_chr1dup.junc.sub.geneselect.pos.D$start%in%df_chr1dup.junc.sub.geneselect.pos.DAwithD.A$start,]
df_chr1dup.junc.sub.geneselect.pos.A.withDA<-df_chr1dup.junc.sub.geneselect.pos.A[df_chr1dup.junc.sub.geneselect.pos.A$end%in%df_chr1dup.junc.sub.geneselect.pos.DAwithD.A$end,]

##combine above three together##
df_chr1dup.junc.sub.geneselect.pos.D.DA.A<-rbind(df_chr1dup.junc.sub.geneselect.pos.D.withDA,df_chr1dup.junc.sub.geneselect.pos.DAwithD.A,df_chr1dup.junc.sub.geneselect.pos.A.withDA)

df_chr1dup.junc.sub.geneselect.pos.D.DA.A<-df_chr1dup.junc.sub.geneselect.pos.D.DA.A%>%
  group_by(genes)%>%
  arrange(start,end)%>%
#  filter(min(df_chr1dup.junc.sub.geneselect.pos.D.DA.A[[df_chr1dup.junc.sub.geneselect.pos.D.DA.A$anchor=="A","start"]])>max(df_chr1dup.junc.sub.geneselect.pos.D.DA.A[[df_chr1dup.junc.sub.geneselect.pos.D.DA.A$anchor=="D","end"]]))
  
##D end < A start##
df_chr1dup.junc.sub.geneselect.pos.D.DA.A$start<-as.integer(as.character(df_chr1dup.junc.sub.geneselect.pos.D.DA.A$start))
df_chr1dup.junc.sub.geneselect.pos.D.DA.A$end<-as.integer(as.character(df_chr1dup.junc.sub.geneselect.pos.D.DA.A$end))
df_chr1dup.junc.sub.geneselect.pos.D.DA.A$genes<-as.character(df_chr1dup.junc.sub.geneselect.pos.D.DA.A$genes)
genepos.sub<-NULL
for (i in 1:length(unique(df_chr1dup.junc.sub.geneselect.pos.D.DA.A$genes))){
  genepos<-unique(df_chr1dup.junc.sub.geneselect.pos.D.DA.A$genes)
  DAs<-unlist(df_chr1dup.junc.sub.geneselect.pos.D.DA.A[df_chr1dup.junc.sub.geneselect.pos.D.DA.A$genes==genepos[i]&df_chr1dup.junc.sub.geneselect.pos.D.DA.A$anchor=="DA","start"])
  DAe<-unlist(df_chr1dup.junc.sub.geneselect.pos.D.DA.A[df_chr1dup.junc.sub.geneselect.pos.D.DA.A$genes==genepos[i]&df_chr1dup.junc.sub.geneselect.pos.D.DA.A$anchor=="DA","end"])
  Ds<-unlist(df_chr1dup.junc.sub.geneselect.pos.D.DA.A[df_chr1dup.junc.sub.geneselect.pos.D.DA.A$genes==genepos[i]&df_chr1dup.junc.sub.geneselect.pos.D.DA.A$anchor=="D","start"])
  De<-unlist(df_chr1dup.junc.sub.geneselect.pos.D.DA.A[df_chr1dup.junc.sub.geneselect.pos.D.DA.A$genes==genepos[i]&df_chr1dup.junc.sub.geneselect.pos.D.DA.A$anchor=="D","end"])
  As<-unlist(df_chr1dup.junc.sub.geneselect.pos.D.DA.A[df_chr1dup.junc.sub.geneselect.pos.D.DA.A$genes==genepos[i]&df_chr1dup.junc.sub.geneselect.pos.D.DA.A$anchor=="A","start"])
  Ae<-unlist(df_chr1dup.junc.sub.geneselect.pos.D.DA.A[df_chr1dup.junc.sub.geneselect.pos.D.DA.A$genes==genepos[i]&df_chr1dup.junc.sub.geneselect.pos.D.DA.A$anchor=="A","end"])
  if (De<max(As))genepos.sub[i]<-genepos[i]
  
}

##################strand (-)#####################################################

##find genes with junctions: if (-): DAstart==Astart,DAend==Dend, DAstart<Aend<Dstart<DAend:
df_chr1dup.junc.sub.geneselect.neg<-df_chr1dup.junc.sub.geneselect[df_chr1dup.junc.sub.geneselect$strand=="-",]
##make sure start and end are integer ################

df_chr1dup.junc.sub.geneselect.neg$start<-as.integer(as.character(df_chr1dup.junc.sub.geneselect.neg$start))
df_chr1dup.junc.sub.geneselect.neg$end<-as.integer(as.character(df_chr1dup.junc.sub.geneselect.neg$end))

##subset neg with DA,D,and A
df_chr1dup.junc.sub.geneselect.neg.DA<-df_chr1dup.junc.sub.geneselect.neg[df_chr1dup.junc.sub.geneselect.neg$anchor=="DA",]
df_chr1dup.junc.sub.geneselect.neg.D<-df_chr1dup.junc.sub.geneselect.neg[df_chr1dup.junc.sub.geneselect.neg$anchor=="D",]
df_chr1dup.junc.sub.geneselect.neg.A<-df_chr1dup.junc.sub.geneselect.neg[df_chr1dup.junc.sub.geneselect.neg$anchor=="A",]

##find overlap DA D and A with start and end position##
df_chr1dup.junc.sub.geneselect.neg.DAwithD.A<-df_chr1dup.junc.sub.geneselect.neg.DA[(df_chr1dup.junc.sub.geneselect.neg.DA$start%in%df_chr1dup.junc.sub.geneselect.neg.A$start)&(df_chr1dup.junc.sub.geneselect.neg.DA$end%in%df_chr1dup.junc.sub.geneselect.neg.D$end),]
df_chr1dup.junc.sub.geneselect.neg.D.withDA<-df_chr1dup.junc.sub.geneselect.neg.D[df_chr1dup.junc.sub.geneselect.neg.D$end%in%df_chr1dup.junc.sub.geneselect.neg.DAwithD.A$end,]
df_chr1dup.junc.sub.geneselect.neg.A.withDA<-df_chr1dup.junc.sub.geneselect.neg.A[df_chr1dup.junc.sub.geneselect.neg.A$start%in%df_chr1dup.junc.sub.geneselect.neg.DAwithD.A$start,]

##combine above three together:A,DA,D##
df_chr1dup.junc.sub.geneselect.neg.D.DA.A<-rbind(df_chr1dup.junc.sub.geneselect.neg.A.withDA,df_chr1dup.junc.sub.geneselect.neg.DAwithD.A,df_chr1dup.junc.sub.geneselect.neg.D.withDA)

##order##
df_chr1dup.junc.sub.geneselect.neg.D.DA.A<-df_chr1dup.junc.sub.geneselect.neg.D.DA.A%>%
  group_by(genes)%>%
  arrange(start,end)%>%

##A end < D start##
df_chr1dup.junc.sub.geneselect.neg.D.DA.A$start<-as.integer(as.character(df_chr1dup.junc.sub.geneselect.neg.D.DA.A$start))
df_chr1dup.junc.sub.geneselect.neg.D.DA.A$end<-as.integer(as.character(df_chr1dup.junc.sub.geneselect.neg.D.DA.A$end))
df_chr1dup.junc.sub.geneselect.neg.D.DA.A$genes<-as.character(df_chr1dup.junc.sub.geneselect.neg.D.DA.A$genes)
geneneg.sub<-NULL
for (i in 1:length(unique(df_chr1dup.junc.sub.geneselect.neg.D.DA.A$genes))){
  geneneg<-unique(df_chr1dup.junc.sub.geneselect.neg.D.DA.A$genes)
  DAs<-unlist(df_chr1dup.junc.sub.geneselect.neg.D.DA.A[df_chr1dup.junc.sub.geneselect.neg.D.DA.A$genes==geneneg[i]&df_chr1dup.junc.sub.geneselect.neg.D.DA.A$anchor=="DA","start"])
  DAe<-unlist(df_chr1dup.junc.sub.geneselect.neg.D.DA.A[df_chr1dup.junc.sub.geneselect.neg.D.DA.A$genes==geneneg[i]&df_chr1dup.junc.sub.geneselect.neg.D.DA.A$anchor=="DA","end"])
  Ds<-unlist(df_chr1dup.junc.sub.geneselect.neg.D.DA.A[df_chr1dup.junc.sub.geneselect.neg.D.DA.A$genes==geneneg[i]&df_chr1dup.junc.sub.geneselect.neg.D.DA.A$anchor=="D","start"])
  De<-unlist(df_chr1dup.junc.sub.geneselect.neg.D.DA.A[df_chr1dup.junc.sub.geneselect.neg.D.DA.A$genes==geneneg[i]&df_chr1dup.junc.sub.geneselect.neg.D.DA.A$anchor=="D","end"])
  As<-unlist(df_chr1dup.junc.sub.geneselect.neg.D.DA.A[df_chr1dup.junc.sub.geneselect.neg.D.DA.A$genes==geneneg[i]&df_chr1dup.junc.sub.geneselect.neg.D.DA.A$anchor=="A","start"])
  Ae<-unlist(df_chr1dup.junc.sub.geneselect.neg.D.DA.A[df_chr1dup.junc.sub.geneselect.neg.D.DA.A$genes==geneneg[i]&df_chr1dup.junc.sub.geneselect.neg.D.DA.A$anchor=="A","end"])
  if (Ae<max(Ds))geneneg.sub[i]<-geneneg[i]
  
}



genecomb<-na.omit(c(genepos.sub,geneneg.sub)) #13





