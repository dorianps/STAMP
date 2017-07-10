args=commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text=args[i]))
}



load( DATAfile )

lapply(c('ggplot2', 'permute', 'randomForest', 'reshape2', 'magrittr','multtest','plyr'),require, character.only=T)
source('/data/jag/dpustina/Code/APHASIA/AC_pRF.R')
seed = sample(100:10000,1)
setwd(savedir)
suppressMessages( p.test<-pRF(response=y,mtry=mtry, predictors=x,n.perms=nperms,ntree=ntree,seed=seed, type="regression",alpha=alpha) )
thisresult =  rep(NA, ncol(x))
thisresult[order(p.test$Res.table$p.value)] = 1:ncol(x)

sig = p.test$Res.table$p.value < 0.05
thissignificant =  rep(0, ncol(x))
thissignificant[sig] = 1
runid = paste(sample(0:9,8, replace=T), collapse='')

write.table(thisresult, col.names=F, row.names=F, file = file.path(savedir,paste0('selectVarsRF.RESULT.result.',jobid, '.',runid, '.csv'))  )
write.table(thissignificant, col.names=F, row.names=F, file = file.path(savedir,paste0('selectVarsRF.RESULT.significant.',jobid, '.',runid, '.csv'))  )

