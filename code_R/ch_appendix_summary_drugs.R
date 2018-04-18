# load data_master and working directory first...
# selecting drugs: Summary statistics:
colnames(traindata$DrugResponse) <- paste0("drug ", 1:31)
stat.desc(traindata$DrugResponse[1:35,], basic=TRUE)[c(1,3:6,8:9,12:13),]
xtable(stat.desc(traindata$DrugResponse[1:35,], basic=TRUE)[c(1,3:6,8:9,12:13),1:16])
xtable(stat.desc(traindata$DrugResponse[1:35,], basic=TRUE)[c(1,3:6,8:9,12:13),1:16], digits=0)
#
xtable(stat.desc(traindata$DrugResponse[1:35,], basic=TRUE)[c(1,3:6,8:9,12:13),17:31])
xtable(stat.desc(traindata$DrugResponse[1:35,], basic=TRUE)[c(1,3:6,8:9,12:13),17:31], digits=0)
