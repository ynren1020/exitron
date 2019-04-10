##############2019-04-05##################
##apply novel_juncfinal function##########
##to a list of janno.gz files#############
##########################################

source("./novel_junc_finalsimple.r")

files<-list.files(pattern="*.janno.gz")
myfiles = do.call(rbind, lapply(files, novel_juncfinal))

write.table(myfiles,"THYM.txt",quote=FALSE,col.names = TRUE,row.names = FALSE)