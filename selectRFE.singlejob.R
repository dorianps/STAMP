args=commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text=args[i]))
}



load( DATAfile )

lapply(c('caret', 'randomForest'),require, character.only=T)

if (RFE$ctrl$allowParallel==T) {
  library(doMC)
  registerDoMC( cores=RFE$cores )
}

RFE$subsets = RFE$subsets[which(RFE$subsets > 1)] # remove subsets==1

rfProfile = rfe(x=RFE$x,
                y=RFE$y,
                rfeControl=RFE$ctrl,
                sizes = RFE$subsets,
                metric=RFE$metric)

RFE$rfProfile = rfProfile
save(RFE, file = DATAfile)

