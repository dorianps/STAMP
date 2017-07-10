#' this should be the file used to create full split cross validations

stop('in case is sourced by mistake')
library(caret)

behav = SD$PNTcorrect
# superfolds = createFolds(behav,k=10)

subsets = c(2:30, seq(32,100,by=3), seq(100,268,by=5))
metric = 'RMSE'

rfefiles = rfefiles.WEIGHTED.REST

rfeCV = list()
dummyrfefiles = rfefiles
dummyrfefiles$Deg = dummyrfefiles$Deg[1]
dummyrfefiles$Bwn = dummyrfefiles$Bwn[1]
dummyrfefiles$LoT = dummyrfefiles$LoT[1]
dummyrfefiles$Eff = dummyrfefiles$Eff[1]
for (ffold in 1:10) {
  inTrain = (1:53)[ -superfolds[[ffold]] ]
  eval(parse(text=paste0('rfeCV$Fold',ffold, '=dummyrfefiles') ))
  
  for (mes in 1:length(rfefiles)) {
    for (dens in 1:10) {
      load(rfefiles[[mes]]$PNTcorrect[[dens]]$DATAfile)
      
      output  = selectRFE(x=RFE$x[inTrain,], y=RFE$y[inTrain], ctrl = RFE$ctrl, subsets = subsets, metric = metric, cores=1, noprint=T, memory = '5G')
      rfeCV[[ffold]][[mes]]$PNTcorrect[[dens]] = output
      cat(paste(dens,''))      
    }
  }
}

#rfeCV.DTI = rfeCV
#rfeCV.REST = rfeCV
save(rfeCV.DTI, rfeCV.REST, superfolds, file='superfolds.RData')



# damage rfeCV for PNTcorrect
######################################
rfdata = groupdmg
subsets = c(1:30,seq(33,72, by=3))
b = 1 # PNTcorrect or use for loop

rfeCV.DMG = list()
behav=eval(parse(text=paste0('SD$',behavs[b])))
for (ffold in 1:10) {
  
  inTrain = (1:53)[ -superfolds[[ffold]] ]
  
  output = selectRFE(x=rfdata[inTrain,],
                     y=behav[inTrain],
                     ctrl = ctrl, 
                     subsets = subsets, 
                     metric = metric, 
                     cores=1, 
                     noprint=T, memory = '2G')

  variable = paste0('rfeCV.DMG$Fold',ffold, '$',behavs[b],'=dummyrfefiles')
  eval(parse(text=paste0( variable, ' = output' )))
}

save(rfeCV.DTI, rfeCV.REST, rfeCV.DMG, superfolds, file='superfolds.RData')
##################################




#################################
subsets = c(seq(5,30, by=5),
            seq(40,100, by=10), 
            seq(130,210,by=30),
            seq(300, 1000, by=100),
            seq(3000,10000,by=2000),
            seq(10000,35000, by=5000))
ctrl2 = ctrl
ctrl2$allowParallel = T

b=1 # for PNTcorrect
behav=eval(parse(text=paste0('SD$',behavs[b])))


rfdata = groupmatDTI[ , apply(groupmatDTI,2,sd)>0 ]
rfdata = groupmat

rfeCV.groupsel = list()
for (ffold in 1:10) {
  inTrain = (1:53)[ -superfolds[[ffold]] ]
  
  output = selectRFE(x=rfdata[inTrain,], 
                     y=behav[inTrain], 
                     ctrl = ctrl2, 
                     subsets = subsets, 
                     metric = metric, 
                     cores=2,
                     noprint=T,
                     memory = '5G')
  
  
  variable = paste0('rfeCV.groupsel$Fold',ffold, '$',behavs[b])
  eval(parse(text=paste0( variable, ' = output' )))
}

rfeCV.groupsel.DTI = rfeCV.groupsel
rfeCV.groupsel.REST = rfeCV.groupsel

save(rfeCV.groupsel.DTI, rfeCV.groupsel.REST, rfeCV.DTI, rfeCV.REST, rfeCV.DMG, superfolds, file='superfolds.RData')
##########################################







# full CV prediction
######################################

getCVpreds = function(RFEmaster, loadfile, inTest, recalcTolerance = F) {
  load(loadfile)
  
  if (recalcTolerance) {
    optsize = RFE$ctrl$functions$selectSize(RFE$rfProfile$results, 'RMSE', tol = 1, maximize = F )
    varsel = RFE$ctrl$functions$selectVar(RFE$rfProfile$variables, optsize)
  } else {
    varsel = RFE$rfProfile$optVariables
  }
  
  rfm = randomForest(x = RFEmaster$x[ -inTest , varsel  ],
                     y = RFEmaster$y[ -inTest ],
                     ntree=500)
  output = predict(rfm, newdata = RFEmaster$x[ inTest , varsel ])
  return(output)
}

getResultsAllThreshCV <- function(behavname, dtirest, superfolds, mes, densities, SD, ...) {

  if (dtirest=='DTI') {
    rfefiles=rfefiles.WEIGHTED.DTI
    rfefiles.fold = rfeCV.DTI
  }
  if (dtirest=='REST') {
    rfefiles=rfefiles.WEIGHTED.REST
    rfefiles.fold = rfeCV.REST
  }
  
  # first work on raw values and get a final prediction from all densities
  denspreds = matrix(NA, nrow=nrow(SD), ncol=length(densities))
  colnames(denspreds) = paste0('Dens_', densities)
  for (d in 1:length(densities)) {
    dens=densities[d]
    load( eval(parse(text=paste0('rfefiles$',mes,'$',behavname,'$Dens',dens,'$DATAfile') )) )
    RFEmaster = RFE
    
    for (ffold in 1:length(superfolds)) {
      denspreds[ superfolds[[ffold]] , d ] =
      getCVpreds(RFEmaster,
                 loadfile=eval(parse(text=paste0('rfefiles.fold$Fold', ffold, '$',mes,'$',behavname,'$Dens',dens,'$DATAfile'))), 
                 inTest = superfolds[[ffold]], 
                 recalcTolerance = recalcTolerance )
    }
  }
  
  resp = eval(parse(text=paste0('SD$',behavname)))
  dfram = data.frame(resp=resp, denspreds)
  finpred = rep(NA,nrow(dfram))
  for (fold in superfolds) {
    rfm = randomForest(resp ~ ., data=dfram[-fold,])
    finpred[fold] = predict(rfm, newdata = dfram[fold,])
  }
  return(finpred)
}


# checking all are done
measures=c('Deg', 'Bwn', 'LoT', 'Eff')
rfeCV = rfeCV.DTI
for (ii in 1:10) {
  for (dd in 1:10) {
    for (mes in measures) {
      load( eval(parse(text=paste0('rfeCV$Fold', ii ,'$', mes ,'$PNTcorrect$Dens', densities[dd], '$DATAfile'))) )
      if (!'rfProfile' %in% names(RFE)) {
        print(paste('missing fold', ii, 'measure', mes, 'density', densities[dd] ))
       # system(RFE$shcmd)
      }
    }
  }
}
  

i=1 # PNTcorrect
behav=eval(parse(text=paste0('SD$',behavs[i])))

preds = data.frame(
  DTIdeg=rep(NA,length(behav)),
  DTIbwn=rep(NA,length(behav)),
  DTIlot=rep(NA,length(behav)),
  DTIeff=rep(NA,length(behav)),
  DTImat=rep(NA,length(behav)),
  RESTdeg=rep(NA,length(behav)),
  RESTbwn=rep(NA,length(behav)),
  RESTlot=rep(NA,length(behav)),
  RESTeff=rep(NA,length(behav)),
  RESTmat=rep(NA,length(behav)),
  DMG=rep(NA,length(behav)),
  LesionSize=les.size
)



recalcTolerance = T

# damage
load( eval(parse(text=paste0('rfefiles.DMG$',behavs[i],'$DATAfile') ))  )
RFEmaster = RFE
for (ffold in 1:10) 
  preds$DMG[superfolds[[ffold]]] = 
  getCVpreds(RFEmaster, loadfile=eval(parse(text=paste0('rfeCV.DMG$Fold', ffold, '$PNTcorrect$DATAfile'))), inTest = superfolds[[ffold]], recalcTolerance = recalcTolerance )


load( eval(parse(text=paste0('groupselDTI$',behavs[i]) )) )
RFEmaster = RFE
for (ffold in 1:10) 
  preds$DTImat[superfolds[[ffold]]] = 
  getCVpreds(RFEmaster, loadfile=eval(parse(text=paste0('rfeCV.groupsel.DTI$Fold', ffold, '$PNTcorrect$DATAfile'))), inTest = superfolds[[ffold]], recalcTolerance = recalcTolerance )


load( eval(parse(text=paste0('groupselREST$',behavs[i]) )) )
RFEmaster = RFE
for (ffold in 1:10) 
  preds$RESTmat[superfolds[[ffold]]] = 
  getCVpreds(RFEmaster, loadfile=eval(parse(text=paste0('rfeCV.groupsel.REST$Fold', ffold, '$PNTcorrect$DATAfile'))), inTest = superfolds[[ffold]], recalcTolerance = recalcTolerance )


dtirest='DTI'; mes='Deg'
eval(parse(text= paste0('preds$', dtirest, tolower(mes), ' = getResultsAllThreshCV(behavs[i], dtirest, superfolds, mes, densities, SD)') ))

dtirest='DTI'; mes='Bwn'
eval(parse(text= paste0('preds$', dtirest, tolower(mes), ' = getResultsAllThreshCV(behavs[i], dtirest, superfolds, mes, densities, SD)') ))

dtirest='DTI'; mes='LoT'
eval(parse(text= paste0('preds$', dtirest, tolower(mes), ' = getResultsAllThreshCV(behavs[i], dtirest, superfolds, mes, densities, SD)') ))

dtirest='DTI'; mes='Eff'
eval(parse(text= paste0('preds$', dtirest, tolower(mes), ' = getResultsAllThreshCV(behavs[i], dtirest, superfolds, mes, densities, SD)') ))

dtirest='REST'; mes='Deg'
eval(parse(text= paste0('preds$', dtirest, tolower(mes), ' = getResultsAllThreshCV(behavs[i], dtirest, superfolds, mes, densities, SD)') ))

dtirest='REST'; mes='Bwn'
eval(parse(text= paste0('preds$', dtirest, tolower(mes), ' = getResultsAllThreshCV(behavs[i], dtirest, superfolds, mes, densities, SD)') ))

dtirest='REST'; mes='LoT'
eval(parse(text= paste0('preds$', dtirest, tolower(mes), ' = getResultsAllThreshCV(behavs[i], dtirest, superfolds, mes, densities, SD)') ))

dtirest='REST'; mes='Eff'
eval(parse(text= paste0('preds$', dtirest, tolower(mes), ' = getResultsAllThreshCV(behavs[i], dtirest, superfolds, mes, densities, SD)') ))



finpred = rep(NA, nrow(preds))
for (fold in superfolds) {
  rfm = randomForest(x=preds[-fold,], y=behav[-fold], ntree=500)
  finpred[fold] = predict(rfm,preds[fold,])
}

cor(behav,preds)
cor(behav,finpred)


# with regular tolerance=4
# > cor(behav,finpred)
# [1] 0.5158103
# > cor(behav,preds)
# DTIdeg    DTIbwn    DTIlot    DTIeff   DTImat   RESTdeg   RESTbwn   RESTlot    RESTeff    RESTmat       DMG LesionSize
# [1,] 0.2011923 0.2477569 0.2035954 0.4030301 0.370943 0.4923537 0.2374056 0.1041428 0.06697471 -0.1013272 0.2581951 -0.4691368

# with reduced tolerance=2
# > cor(behav,preds)
# DTIdeg    DTIbwn    DTIlot    DTIeff    DTImat   RESTdeg   RESTbwn   RESTlot   RESTeff    RESTmat       DMG LesionSize
# [1,] 0.362159 0.2594989 0.3798503 0.6294335 0.4191676 0.5306434 0.3697277 0.1314792 0.2005833 -0.1465812 0.3823072 -0.4691368
# > cor(behav,finpred)
# [1] 0.6586547