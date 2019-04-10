#############2019-04-09####################
##combine junction number and sample info##
##e.g. ACC.txt and ACC.info.txt############
##for each cohort##########################
##and then all cohort######################
###########################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

##TCGA cohort##
cohort.all<-c("ACC","BLCA","BRCA","CESC", "CHOL", "COAD", "DLBC",  "ESCA",  "GBM",  "HNSC",  "KICH",  "KIRC", "KIRP",  "LAML",  "LGG",  "LIHC",  "LUAD",  "LUSC",  "MESO",  "OV",  "PAAD",  "PCPG",  "PRAD",  "READ",  "SARC","SKCM",  "STAD",  "TGCT",  "THCA",  "THYM",  "UCEC",  "UCS",  "UVM")

##function to combine junction and info##
dat4test6<-function(input1,input2){
  dat<-read.delim(input1,header = TRUE,sep = " ",stringsAsFactors = FALSE)%>%rename(FILE_ID=sample)
  for (i in 1:nrow(dat)){
    dat$FILE_ID[i]<-strsplit(dat$FILE_ID[i],"[.]")[[1]][1]
  }
  info<-read.delim(input2,header=TRUE,sep="\t",stringsAsFactors = FALSE)
  comb<-full_join(dat,info,by="FILE_ID")
  return(comb)
  
}

dat4test<-function(input1,input2){
  dat<-read.delim(input1,header = FALSE,sep = " ",stringsAsFactors = FALSE)%>%rename(FILE_ID=V1,novel=V2,total=V3)
  for (i in 1:nrow(dat)){
    dat$FILE_ID[i]<-strsplit(dat$FILE_ID[i],"[.]")[[1]][1]
  }
  info<-read.delim(input2,header=TRUE,sep="\t",stringsAsFactors = FALSE)
  comb<-full_join(dat,info,by="FILE_ID")
  return(comb)
  
}

##how to apply function to multiple lists##
#cohort<-c("ACC","DLBC","TGCT","THYM", "UCS",  "UVM")
#cohort.comb<-array(,dim=c(,6,6))
#for (i in 1:length(cohort)){
# cohort.comb[i]<-array(,dim=c(nrow(paste0(cohort[i],".txt")),6,6))
# cohort.comb[i]<-dat4test(paste0(cohort[i],".txt"),paste0(cohort[i],".info.txt"))
#}


##six cohort produced in R##
ACC<-dat4test6("ACC.txt","ACC.info.txt")
DLBC<-dat4test6("DLBC.txt","DLBC.info.txt")
TGCT<-dat4test6("TGCT.txt","TGCT.info.txt")
THYM<-dat4test6("THYM.txt","THYM.info.txt")
UCS<-dat4test6("UCS.txt","UCS.info.txt")
UVM<-dat4test6("UVM.txt","UVM.info.txt")

##27cohort produced in linux##
BLCA<-dat4test("BLCA.janno.txt","BLCA.info.txt")
BRCA<-dat4test("BRCA.janno.txt","BRCA.info.txt")
CESC<-dat4test("CESC.janno.txt","CESC.info.txt")
CHOL<-dat4test("CHOL.janno.txt","CHOL.info.txt")
COAD<-dat4test("COAD.janno.txt","COAD.info.txt")
ESCA<-dat4test("ESCA.janno.txt","ESCA.info.txt")
GBM<-dat4test("GBM.janno.txt","GBM.info.txt")
HNSC<-dat4test("HNSC.janno.txt","HNSC.info.txt")
KICH<-dat4test("KICH.janno.txt","KICH.info.txt")
KIRC<-dat4test("KIRC.janno.txt","KIRC.info.txt")
KIRP<-dat4test("KIRP.janno.txt","KIRP.info.txt")
LAML<-dat4test("LAML.janno.txt","LAML.info.txt")
LGG<-dat4test("LGG.janno.txt","LGG.info.txt")
LIHC<-dat4test("LIHC.janno.txt","LIHC.info.txt")
LUAD<-dat4test("LUAD.janno.txt","LUAD.info.txt")
LUSC<-dat4test("LUSC.janno.txt","LUSC.info.txt")
MESO<-dat4test("MESO.janno.txt","MESO.info.txt")
OV<-dat4test("OV.janno.txt","OV.info.txt")
PAAD<-dat4test("PAAD.janno.txt","PAAD.info.txt")
PCPG<-dat4test("PCPG.janno.txt","PCPG.info.txt")
PRAD<-dat4test("PRAD.janno.txt","PRAD.info.txt")
READ<-dat4test("READ.janno.txt","READ.info.txt")
SARC<-dat4test("SARC.janno.txt","SARC.info.txt")
SKCM<-dat4test("SKCM.janno.txt","SKCM.info.txt")
STAD<-dat4test("STAD.janno.txt","STAD.info.txt")
THCA<-dat4test("THCA.janno.txt","THCA.info.txt")
UCEC<-dat4test("UCEC.janno.txt","UCEC.info.txt")


cohortdf<-rbind(ACC,BLCA,BRCA,CESC, CHOL, COAD, DLBC,  ESCA,  GBM,  HNSC,  KICH,  KIRC, KIRP,  LAML,  LGG,  LIHC,  LUAD,  LUSC,  MESO,  OV,  PAAD,  PCPG,  PRAD,  READ,  SARC,SKCM,  STAD,  TGCT,  THCA,  THYM,  UCEC,  UCS,  UVM)

for (i in 1:nrow(cohortdf)){
  cohortdf$PROJECT_f[i]<-strsplit(cohortdf$PROJECT[i],"-")[[1]][2]
}

cohortdf$PROJECT_f<-as.factor(cohortdf$PROJECT_f)
cohortdf$TUMOR_TYPE<-as.factor(cohortdf$TUMOR_TYPE)

table(cohortdf$TUMOR_TYPE)
#Additional - New Primary                           Additional Metastatic 
#11                                               1 
#Metastatic Primary Blood Derived Cancer - Peripheral Blood 
#392                                             151 
#Primary Tumor                                 Recurrent Tumor 
#9760                                              47 
#Solid Tissue Normal 
#730 

##separate sample into tumor and normal##
cohortdf<-na.omit(cohortdf)
cohortdf$status<-NULL
for (i in 1:nrow(cohortdf)){
  cohortdf$status[i]<-ifelse(cohortdf$TUMOR_TYPE[i]=="Solid Tissue Normal",0,1)
}
cohortdf$status<-as.factor(cohortdf$status)

##poisson regression for novel##
fit1<-glm(novel~total+PROJECT_f+status,family="poisson", data=cohortdf)
sink("poisson-fit1.txt")
summary(fit1)
sink()

##plot##
cohort<-c('ACC','BLCA','BRCA', 'CESC', 'CHOL','COAD','DLBC','ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC',  'LUAD', 'LUSC',
          'MESO', 'OV', 'PAAD','PCPG', 'PRAD','READ','SARC','SKCM','STAD', 'TGCT', 'THCA', 'THYM','UCEC', 'UCS', 'UVM')
color<-c('#7b6a4a', '#8049d8',  '#5dd74f', '#cf45cc','#dde73f','#646acd', '#a0d64b','#d3458d', '#67d88e','#dd4529', '#7fdfcd', '#d34058', '#568734',
         '#cd83d9','#bebe49', '#964891', '#dfb23f', '#729bd9','#db842d','#51acc0','#b25c35','#a7cbda', '#8f7b30', '#516590','#c5d994', '#9d5356',
         '#569973', '#dd8999', '#4b7267','#daa476','#d0aed7', '#ccbeaa','#886b82')
cohortcolor<-data_frame(cohort=cohort,color=color)
cohortdf<-cohortdf%>%rename(cohort=PROJECT_f)
cohortdf$cohort<-as.character(cohortdf$cohort)

junctionplot<-left_join(cohortdf,cohortcolor,by="cohort") #11096
##tumor##
junctionplot$status<-as.numeric(as.character(junctionplot$status))
junctionplot<-junctionplot[junctionplot$status==1,] #10215
##normal##
junctionplot<-junctionplot[junctionplot$status==0,]#730

##create rank of number and generate median of it within each cohort##
junctionplot<-junctionplot %>%
  group_by(cohort) %>%
  mutate(my_ranks = order(order(novel, decreasing=FALSE)))%>%
  mutate(my_median = median(novel,na.rm=TRUE))%>%
  mutate(my_mean_ranks=mean(my_ranks))
##order by median##
junctionplot<-as.data.frame(junctionplot)
junctionplot<-junctionplot[order(junctionplot$my_median),]
##create factor variable of cohort, rename color##
junctionplot$cohort_f<-factor(junctionplot$cohort,levels=unique(junctionplot$cohort))
junctionplot<-rename(junctionplot,my_color=color)

##plot##
p<-ggplot(junctionplot,aes(x=my_ranks, y=novel,color=cohort_f))+
  geom_point(size=0.5)+
  scale_colour_manual(values=unique(junctionplot$my_color))+
  labs(y="Number of junction")


p1<-p+facet_grid(.~cohort_f,switch="x",scales = "free_x")+
  #geom_hline(data=data2,aes(yintercept=my_median),linetype="dotted", color = "red")+
  geom_segment(aes(x=median(my_ranks)-3,xend=median(my_ranks)+3,y=my_median,yend=my_median),colour="red",size=1)+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(),legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.spacing = unit(0, "mm"),panel.border = element_blank(),                      
        strip.background = element_blank(),axis.line.x = element_line(),axis.line.y = element_line(),strip.text.x = element_text(angle = 90),strip.placement = "outside")

p1







