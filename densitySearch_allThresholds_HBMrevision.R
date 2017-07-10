#' I think this needs some variables from the original
#'  "restAnalysis_v2" file.
#'  
#'  here is the first load code from it
source('/data/jag/dpustina/Code/APHASIA/loadExcel.R', echo=F)
restppl =  SD$rest1==1 & !is.na(SD$rest1)
SD = SD[ restppl, ]
restnames = SD$NameOnDisk



# just control and other RFE parameters
#####
require('caret')
require(Metrics)
selectVarRerank <- function (y, size) 
{
  y = y[ y$Variables==size,]
  finalImp <- ddply(y[, c("Overall", "var")], 
                    .(var), 
                    function(x) mean(x$Overall, na.rm = TRUE)
  )
  names(finalImp)[2] <- "Overall"
  finalImp <- finalImp[order(finalImp$Overall, decreasing = TRUE), ]
  as.character(finalImp$var[1:size])
}

rfFuncs2 <- rfFuncs
rfFuncs2$fit <- function (x, y, first, last, ...) {
  loadNamespace("randomForest")
  randomForest::randomForest(x, y, importance=T,...)
}
rfFuncs2$selectVar=selectVarRerank
rfFuncs2$selectSize=pickSizeTolerance
formals(rfFuncs2$selectSize)$tol = 4

repeats=10#500

ctrl = rfeControl(functions = rfFuncs2,
                  method = 'repeatedcv',
                  repeats=repeats,
                  saveDetails=T,
                  returnResamp='all',
                  verbose=F,
                  rerank=T, # recomputes importance after each elimination step
                  allowParallel=F)
metric='RMSE'

source('/data/jag/dpustina/Code/APHASIA/makeGraph.Stroke.R')

# variables
subsets = c(1:30, seq(32,100,by=3), seq(100,268,by=5))
densities = seq(0.1,1,by=0.1)
checks = c('Deg', 'Bwn', 'LoT', 'Eff')
behavs = c('PNTcorrect', 'WABAQ', 'WABrep', 'WABcomp')

rfefiles = list()
for (aa in checks) {
  eval(parse(text=paste0('rfefiles$',aa,'=list()')))
  for (bb in behavs) {
    eval(parse(text=paste0('rfefiles$',aa,'$', bb ,'=list()')))
    for (ss in densities) {
      eval(parse(text=paste0('rfefiles$',aa,'$', bb ,'$Dens', ss , '=list()')))
    }
  }
}
  

#####



# big for density here SEND ALL JOBS
#####
for (dd in 1:length(densities)) {
  density=densities[dd]
  # build graph theyr matrix

  graphmatDTI = list()
  for (i in 1:nrow(SD)) {
    temp = makeGraph(abs(submats[[i]]), graphdensity = density, getEfficiency = T)
    graphmatDTI[[i]] = c(
      temp$degree, 
      temp$betweeness, 
      temp$localtransitivity,
      temp$effinv,
      temp$globalTransitivity
    )
    cat(paste(i,''))
  }
  graphmatpredDTI = t(sapply(graphmatDTI, rbind ))
  graphmatpredDTI[is.na(graphmatpredDTI)]=0 # some some regions have local transitivity= NA
  colnames(graphmatpredDTI) = c(
    paste0('Deg_',1:nnodes),
    paste0('Bwn_',1:nnodes),
    paste0('LoT_',1:nnodes),
    paste0('Eff_',1:nnodes),
    paste0('GloT',1)
  )
  # MATRIX BUILT
  
  
  # loop through graph measures
  
  # for ch here
  for (ch in 1:length(checks)) {
    rfdata=graphmatpredDTI
    
    
    checkfor = checks[ch]
    sel= grep(checkfor, colnames(rfdata) )
    rfdata = as.matrix(rfdata[,sel])
    cat(paste("\nSD rfdata", sd(rfdata)))
    if (ncol(rfdata) != 268) stop('Non 268 columns')
    
    # for behav here
    for (b in 1:length(behavs) ) {
      behav=eval(parse(text=paste0('SD$',behavs[b])))
      inTrain = (1:53)
    
      rafaja = selectRFE(x=rfdata[inTrain,],
                               y=behav[inTrain],
                               ctrl = ctrl, 
                               subsets = subsets, 
                               metric = metric, 
                               cores=1, 
                               noprint=T, memory = '1.5G')
      variable = paste0('rfefiles$',checkfor,'$', behavs[b] ,'$Dens', density)
      eval(parse(text=paste0( variable, ' = rafaja' )))

      
    }
  }
}
rfefiles.WEIGHTED.DTI = rfefiles
rfefiles.WEIGHTED.REST = rfefiles
rfefiles.BINARY.REST = rfefiles
save(rfefiles.WEIGHTED.DTI,rfefiles.WEIGHTED.REST,rfefiles.BINARY.REST,
     file='rfefiles.WEIGHTED.Rdata')
#####



# checking if they finished
#####
rfefiles=rfefiles.WEIGHTED.DTI
for (aa in 1:length(rfefiles))
  for (bb in 1:length(rfefiles[[aa]]))
    for (cc in 1:length(rfefiles[[aa]][[bb]]))
      system(rfefiles[[aa]][[bb]][[cc]]$shcmd)


load(rfefiles.WEIGHTED.DTI[[4]][[4]][[9]]$DATAfile)
print(paste(RFE$rfProfile$optsize, paste(RFE$rfProfile$optVariables,collapse=' ') ))
plot(RFE$rfProfile, type=c('g','o'), main=paste('Optimal -', rfProfile$optsize, 'variables') )


rferesults = function(RFE) {
  require(caret)
  require(randomForest)
  varsel=RFE$rfProfile$optVariables
  predgraph.repeated = matrix(NA,nrow=nrow(RFE$x), ncol=10)
  for (i in 1:10) {
    predgraph = rep(NA, nrow(RFE$x))
    for (folds in createFolds(RFE$y,k = 10,list=T)) {
      inTest = folds
      rfm = randomForest(x = as.matrix(RFE$x[ -inTest, varsel]),
                         y = RFE$y[-inTest], 
                         ntree = 500,
      )
      predgraph[inTest] = predict(rfm, newdata = as.matrix(RFE$x[ inTest, varsel])  )      
    }
    predgraph.repeated[ , i] = predgraph
  }
  
  predgraphFULL = rep(NA, nrow(RFE$x))
  varsel=1:ncol(RFE$x)
  for (folds in createFolds(RFE$y,k = 10,list=T)) {
    inTest = folds
    rfm = randomForest(x = as.matrix(RFE$x[ -inTest, varsel]),
                       y = RFE$y[-inTest], 
                       ntree = 500,
    )
    predgraphFULL[inTest] = predict(rfm, newdata = as.matrix(RFE$x[ inTest, varsel])  )
    
  }
  
  require(Metrics)
  rmseSelected =   mean( apply(predgraph.repeated,2, function(x) rmse(RFE$y,x)) )
  rmseFULL = rmse(RFE$y, predgraphFULL)
  
  corSelected = mean(cor(RFE$y, predgraph.repeated))
  corFULL = cor(RFE$y, predgraphFULL)
  corFULLSelected = mean( cor(predgraph.repeated,predgraphFULL) )
  
  require(psych)
  p.FULLSelected=paired.r(corSelected,corFULL, corFULLSelected,n = nrow(RFE$x))$p
  
  
  return(list(
    predgraph=predgraph.repeated,
    predgraphFULL=predgraphFULL,
    rmseSelected=rmseSelected,
    rmseFULL=rmseFULL,
    corSelected=corSelected,
    corFULL=corFULL,
    corFULLSelected=corFULLSelected,
    p.FULLSelected=p.FULLSelected
    ))
}

cor(RFE$y, predgraph)
rmse(RFE$y, predgraph)
plot(RFE$y, predgraph, ylim=range(RFE$y), pch=19); abline(a=0,b=1)
#####




# compute results and resave
#####
rfefiles=rfefiles.WEIGHTED.REST
for (behav in behavs) {
  for (mes in checks) {
    densdata = data.frame(Sel=rep(NA,length(densities)), Unsel=rep(NA,length(densities)) ); rownames(densdata)=densities
    for (dens in densities) {
      loadfile = eval(parse(text=paste0('rfefiles$',mes,'$',behav,'$Dens',dens,'$DATAfile') ))
      eval(parse(text=paste0('load(\'', loadfile, '\')' ) ))
      if (!is.null(RFE$results)) next
      RFE$results=rferesults(RFE)
      RFE$behav = behav
      RFE$measure=mes
      RFE$density=dens
      save(RFE, file=loadfile)

      cat(paste(behav, mes, dens, '\n'))
    }
  }
}
#####



# density list
#####
rfefiles=rfefiles.WEIGHTED.DTI
denslist = list()
for (behav in behavs) {
  eval(parse(text=paste0('denslist$', behav, '=list()')))
  for (mes in checks) {
    eval(parse(text=paste0('denslist$', behav, '$',mes,'=list()')))
    densdata = data.frame(Sel=rep(NA,length(densities)), 
                          Unsel=rep(NA,length(densities)),
                          pval=rep(NA,length(densities)) ) 
    rownames(densdata)=densities
    for (dens in densities) {
      loadfile = eval(parse(text=paste0('rfefiles$',mes,'$',behav,'$Dens',dens,'$DATAfile') ))
      eval(parse(text=paste0('load(\'', loadfile, '\')' ) ))
      dchar=as.character(dens)
      densdata[dchar,'Sel'] = RFE$results$corSelected
      densdata[dchar,'Unsel'] = RFE$results$corFULL
      densdata[dchar,'pval'] = RFE$results$p.FULLSelected
      
      cat(paste(behav, mes, dens, '\n'))
    }
    
    eval(parse(text=paste0('denslist$', behav, '$',mes,'$densdata=densdata' )))
    eval(parse(text=paste0('denslist$', behav, '$',mes,'$behav=behav' )))
    eval(parse(text=paste0('denslist$', behav, '$',mes,'$mes=mes' )))
    
  }
}
#####




# plot the shit
#####
par(mfrow=c(4,4))
for (i in 1:4) {
  mmax = 0
  mmin = 1
  for (j in 1:4) {
    mmax = max(mmax, max(denslist[[i]][[j]]$densdata[,1:2]))
    mmin = min(mmin, min(denslist[[i]][[j]]$densdata[1:2]))
  }
  for (j in 1:4) {
    
    densdata = denslist[[i]][[j]]$densdata
    behav = denslist[[i]][[j]]$behav
    mes = denslist[[i]][[j]]$mes
    
    par(mar=c(3,3,2,1))
    plot(rownames(densdata),densdata$Sel, type='l', col='red',ylim=c(0.1,0.9), #c(mmin-0.1,mmax+0.1), 
         lwd=3, xlab='', ylab='', cex.axis=1.4, yaxt='n',
         main=paste(behav, '-', mes) )
    axis(2, labels=c('', 0.4, '', 0.8), at=c(0.2,0.4,0.6,0.8), cex.axis=1.4, las=2)
    points(x=rownames(densdata), y=densdata$Sel, pch=19, col='black')
    lines(rownames(densdata),densdata$Unsel, col='black', lty='dashed')
    points(x=rownames(densdata), y=densdata$Unsel, pch=1)
    text(x=rownames(densdata)[densdata$pval<0.05], y=0.70, "*", pos=3, cex=1.8,col='brown')
    abline(a=0.27,b=0,lty='dotted')
    
  }
}
#####



# how about combining all predictions from all densities... seriously?
# it doesn't do better than the maximum correlation
#####
rfefiles=rfefiles.WEIGHTED.DTI
for (behav in behavs) {
  eval(parse(text=paste0('denslist$', behav, '=list()')))
  for (mes in checks) {
  dfram=matrix(NA,ncol=length(densities),nrow=53); colnames(dfram) = densities
    for (dens in densities) {
      loadfile = eval(parse(text=paste0('rfefiles$',mes,'$',behav,'$Dens',dens,'$DATAfile') ))
      eval(parse(text=paste0('load(\'', loadfile, '\')' ) ))
      
      dchar=as.character(dens)
      dfram[ , dchar] = rowMeans(RFE$results$predgraph) # [,1]
    }
    dfram = data.frame(dfram, resp=RFE$y)
    
    finpred = rep(NA,nrow(dfram))
    for (fold in createFolds(y=RFE$y,k=10)) {
      rfm = randomForest(resp ~ ., data=dfram[-fold,],mtry=10)
      finpred[fold] = predict(rfm, newdata = dfram[fold,])
    }
    
    print(paste(behav, mes, '- Finpred:',   round(cor(dfram[,'resp'],finpred),2), '- Max:', round(max(cor(dfram$resp, dfram[,1:10])),2) ))

    
  }
}
# DTI weighted
# [1] "PNTcorrect Deg - Finpred: 0.57 - Max: 0.64"
# [1] "PNTcorrect Bwn - Finpred: 0.55 - Max: 0.71"
# [1] "PNTcorrect LoT - Finpred: 0.73 - Max: 0.79"
# [1] "PNTcorrect Eff - Finpred: 0.65 - Max: 0.66"
# [1] "WABAQ Deg - Finpred: 0.56 - Max: 0.59"
# [1] "WABAQ Bwn - Finpred: 0.67 - Max: 0.72"
# [1] "WABAQ LoT - Finpred: 0.71 - Max: 0.72"
# [1] "WABAQ Eff - Finpred: 0.76 - Max: 0.76"
# [1] "WABrep Deg - Finpred: 0.61 - Max: 0.69"
# [1] "WABrep Bwn - Finpred: 0.67 - Max: 0.67"
# [1] "WABrep LoT - Finpred: 0.62 - Max: 0.69"
# [1] "WABrep Eff - Finpred: 0.78 - Max: 0.78"
# [1] "WABcomp Deg - Finpred: 0.37 - Max: 0.56"
# [1] "WABcomp Bwn - Finpred: 0.73 - Max: 0.69"
# [1] "WABcomp LoT - Finpred: 0.56 - Max: 0.71"
# [1] "WABcomp Eff - Finpred: 0.63 - Max: 0.71"
# REST weighted
# [1] "PNTcorrect Deg - Finpred: 0.32 - Max: 0.48"
# [1] "PNTcorrect Bwn - Finpred: 0.71 - Max: 0.63"
# [1] "PNTcorrect LoT - Finpred: 0.24 - Max: 0.56"
# [1] "PNTcorrect Eff - Finpred: 0.54 - Max: 0.49"
# [1] "WABAQ Deg - Finpred: 0.57 - Max: 0.6"
# [1] "WABAQ Bwn - Finpred: 0.68 - Max: 0.6"
# [1] "WABAQ LoT - Finpred: 0.59 - Max: 0.62"
# [1] "WABAQ Eff - Finpred: 0.46 - Max: 0.44"
# [1] "WABrep Deg - Finpred: 0.69 - Max: 0.57"
# [1] "WABrep Bwn - Finpred: 0.79 - Max: 0.71"
# [1] "WABrep LoT - Finpred: 0.64 - Max: 0.62"
# [1] "WABrep Eff - Finpred: 0.47 - Max: 0.47"
# [1] "WABcomp Deg - Finpred: 0.49 - Max: 0.49"
# [1] "WABcomp Bwn - Finpred: 0.65 - Max: 0.64"
# [1] "WABcomp LoT - Finpred: 0.53 - Max: 0.56"
# [1] "WABcomp Eff - Finpred: 0.46 - Max: 0.5"
#####





# how about combining, not predictions, but raw graph values, as selected by RFE
#####
rfefiles=rfefiles.WEIGHTED.REST
for (behav in behavs) {
  eval(parse(text=paste0('denslist$', behav, '=list()')))
  for (mes in checks) {
    dfram=matrix(NA,nrow=53)
    maxcor = 0
    for (dens in densities) {
      loadfile = eval(parse(text=paste0('rfefiles$',mes,'$',behav,'$Dens',dens,'$DATAfile') ))
      eval(parse(text=paste0('load(\'', loadfile, '\')' ) ))
      
      addpreds = RFE$x[ , RFE$rfProfile$optVariables]
      colnames(addpreds) = paste0('D', dens, '_', colnames(addpreds))
      dfram = cbind(dfram, addpreds)
      
      maxcor = ifelse(RFE$results$corSelected > maxcor, RFE$results$corSelected, maxcor)
    }
    dfram = dfram[,-1]
    dfram = data.frame(dfram, resp=RFE$y)
    
    finpred = rep(NA,nrow(dfram))
    for (fold in createFolds(y=RFE$y,k=10)) {
      rfm = randomForest(resp ~ ., data=dfram[-fold,],mtry=10)
      finpred[fold] = predict(rfm, newdata = dfram[fold,])
    }
    
    print(paste(behav, mes, '- Finpred:',   round(cor(RFE$y,finpred),2),
                '- Maxcor:', round(maxcor,2),
                '-', (ncol(dfram)-1) , 'predictors'  ))
  }
}
# # DTIweighted
# [1] "PNTcorrect Deg - Finpred: 0.53 - Maxcor: 0.64 - 115 predictors"
# [1] "PNTcorrect Bwn - Finpred: 0.67 - Maxcor: 0.7 - 227 predictors"
# [1] "PNTcorrect LoT - Finpred: 0.73 - Maxcor: 0.78 - 45 predictors"
# [1] "PNTcorrect Eff - Finpred: 0.57 - Maxcor: 0.65 - 147 predictors"
# [1] "WABAQ Deg - Finpred: 0.59 - Maxcor: 0.59 - 84 predictors"
# [1] "WABAQ Bwn - Finpred: 0.71 - Maxcor: 0.72 - 211 predictors"
# [1] "WABAQ LoT - Finpred: 0.69 - Maxcor: 0.73 - 174 predictors"
# [1] "WABAQ Eff - Finpred: 0.68 - Maxcor: 0.75 - 66 predictors"
# [1] "WABrep Deg - Finpred: 0.64 - Maxcor: 0.69 - 42 predictors"
# [1] "WABrep Bwn - Finpred: 0.67 - Maxcor: 0.67 - 58 predictors"
# [1] "WABrep LoT - Finpred: 0.57 - Maxcor: 0.73 - 163 predictors"
# [1] "WABrep Eff - Finpred: 0.59 - Maxcor: 0.77 - 91 predictors"
# [1] "WABcomp Deg - Finpred: 0.48 - Maxcor: 0.6 - 155 predictors"
# [1] "WABcomp Bwn - Finpred: 0.69 - Maxcor: 0.7 - 103 predictors"
# [1] "WABcomp LoT - Finpred: 0.61 - Maxcor: 0.69 - 203 predictors"
# [1] "WABcomp Eff - Finpred: 0.71 - Maxcor: 0.67 - 32 predictors" #####<
# # REST weighted
# [1] "PNTcorrect Deg - Finpred: 0.19 - Maxcor: 0.46 - 602 predictors"
# [1] "PNTcorrect Bwn - Finpred: 0.74 - Maxcor: 0.57 - 217 predictors" #####<
# [1] "PNTcorrect LoT - Finpred: 0.52 - Maxcor: 0.54 - 199 predictors"
# [1] "PNTcorrect Eff - Finpred: 0.43 - Maxcor: 0.49 - 203 predictors"
# [1] "WABAQ Deg - Finpred: 0.4 - Maxcor: 0.59 - 340 predictors"
# [1] "WABAQ Bwn - Finpred: 0.8 - Maxcor: 0.63 - 228 predictors" #####<
# [1] "WABAQ LoT - Finpred: 0.57 - Maxcor: 0.6 - 173 predictors"
# [1] "WABAQ Eff - Finpred: 0.31 - Maxcor: 0.39 - 305 predictors"
# [1] "WABrep Deg - Finpred: 0.33 - Maxcor: 0.59 - 474 predictors"
# [1] "WABrep Bwn - Finpred: 0.75 - Maxcor: 0.7 - 121 predictors" #####<
# [1] "WABrep LoT - Finpred: 0.53 - Maxcor: 0.59 - 145 predictors"
# [1] "WABrep Eff - Finpred: 0.39 - Maxcor: 0.43 - 159 predictors"
# [1] "WABcomp Deg - Finpred: 0.46 - Maxcor: 0.54 - 136 predictors"
# [1] "WABcomp Bwn - Finpred: 0.76 - Maxcor: 0.56 - 379 predictors" #####<
# [1] "WABcomp LoT - Finpred: 0.44 - Maxcor: 0.55 - 319 predictors"
# [1] "WABcomp Eff - Finpred: 0.19 - Maxcor: 0.48 - 129 predictors"
#####





# combine raw graphs, then input prediction to density preds
# not much benefit
#####
rfefiles=rfefiles.WEIGHTED.REST
for (behav in behavs) {
  eval(parse(text=paste0('denslist$', behav, '=list()')))
  FOLDS = createFolds(y=RFE$y,k=10)
  
  for (mes in checks) {
    # first work on raw values and get a final prediction from all densities
    dfram=matrix(NA,nrow=53)
    maxcor = 0
    for (dens in densities) {
      loadfile = eval(parse(text=paste0('rfefiles$',mes,'$',behav,'$Dens',dens,'$DATAfile') ))
      eval(parse(text=paste0('load(\'', loadfile, '\')' ) ))
      
      addpreds = RFE$x[ , RFE$rfProfile$optVariables]
      colnames(addpreds) = paste0('D', dens, '_', colnames(addpreds))
      dfram = cbind(dfram, addpreds)
      
      maxcor = ifelse(RFE$results$corSelected > maxcor, RFE$results$corSelected, maxcor)
    }
    dfram = dfram[,-1]
    dfram = data.frame(dfram, resp=RFE$y)
    
    finpred = rep(NA,nrow(dfram))
    for (fold in FOLDS) {
      rfm = randomForest(resp ~ ., data=dfram[-fold,],mtry=10)
      finpred[fold] = predict(rfm, newdata = dfram[fold,])
    }
    
        # second, add prediction to density predictions
    dfram=matrix(NA,ncol=length(densities),nrow=53); colnames(dfram) = densities
    for (dens in densities) {
      loadfile = eval(parse(text=paste0('rfefiles$',mes,'$',behav,'$Dens',dens,'$DATAfile') ))
      eval(parse(text=paste0('load(\'', loadfile, '\')' ) ))
      
      dchar=as.character(dens)
      dfram[ , dchar] = rowMeans(RFE$results$predgraph) # [,1]
    }
    dfram = data.frame(dfram, resp=RFE$y)
    dfram = cbind(dfram, rawpred=finpred)
    # this is the final prediction
    finpred = rep(NA,nrow(dfram))
    for (fold in FOLDS) {
      rfm = randomForest(resp ~ ., data=dfram[-fold,],mtry=10)
      finpred[fold] = predict(rfm, newdata = dfram[fold,])
    }
    
    # second, input the prediction to raw 
    
    print(paste(behav, mes, '- Finpred:',   round(cor(RFE$y,finpred),2),
                '- Maxcor:', round(maxcor,2),
                '-', (ncol(dfram)-1) , 'predictors'  ))
  }
}
# DTI weighted
# [1] "PNTcorrect Deg - Finpred: 0.53 - Maxcor: 0.64 - 11 predictors"
# [1] "PNTcorrect Bwn - Finpred: 0.63 - Maxcor: 0.7 - 11 predictors"
# [1] "PNTcorrect LoT - Finpred: 0.77 - Maxcor: 0.78 - 11 predictors"
# [1] "PNTcorrect Eff - Finpred: 0.65 - Maxcor: 0.65 - 11 predictors"
# [1] "WABAQ Deg - Finpred: 0.58 - Maxcor: 0.59 - 11 predictors"
# [1] "WABAQ Bwn - Finpred: 0.67 - Maxcor: 0.72 - 11 predictors"
# [1] "WABAQ LoT - Finpred: 0.7 - Maxcor: 0.73 - 11 predictors"
# [1] "WABAQ Eff - Finpred: 0.69 - Maxcor: 0.75 - 11 predictors"
# [1] "WABrep Deg - Finpred: 0.6 - Maxcor: 0.69 - 11 predictors"
# [1] "WABrep Bwn - Finpred: 0.71 - Maxcor: 0.67 - 11 predictors"
# [1] "WABrep LoT - Finpred: 0.58 - Maxcor: 0.73 - 11 predictors"
# [1] "WABrep Eff - Finpred: 0.72 - Maxcor: 0.77 - 11 predictors"
# [1] "WABcomp Deg - Finpred: 0.53 - Maxcor: 0.6 - 11 predictors"
# [1] "WABcomp Bwn - Finpred: 0.7 - Maxcor: 0.7 - 11 predictors"
# [1] "WABcomp LoT - Finpred: 0.6 - Maxcor: 0.69 - 11 predictors"
# [1] "WABcomp Eff - Finpred: 0.7 - Maxcor: 0.67 - 11 predictors"
# REST weighted
# [1] "PNTcorrect Deg - Finpred: 0.6 - Maxcor: 0.46 - 11 predictors"
# [1] "PNTcorrect Bwn - Finpred: 0.65 - Maxcor: 0.57 - 11 predictors"
# [1] "PNTcorrect LoT - Finpred: 0.34 - Maxcor: 0.54 - 11 predictors"
# [1] "PNTcorrect Eff - Finpred: 0.58 - Maxcor: 0.49 - 11 predictors"
# [1] "WABAQ Deg - Finpred: 0.68 - Maxcor: 0.59 - 11 predictors"
# [1] "WABAQ Bwn - Finpred: 0.75 - Maxcor: 0.63 - 11 predictors"
# [1] "WABAQ LoT - Finpred: 0.54 - Maxcor: 0.6 - 11 predictors"
# [1] "WABAQ Eff - Finpred: 0.29 - Maxcor: 0.39 - 11 predictors"
# [1] "WABrep Deg - Finpred: 0.71 - Maxcor: 0.59 - 11 predictors"
# [1] "WABrep Bwn - Finpred: 0.84 - Maxcor: 0.7 - 11 predictors"
# [1] "WABrep LoT - Finpred: 0.65 - Maxcor: 0.59 - 11 predictors"
# [1] "WABrep Eff - Finpred: 0.42 - Maxcor: 0.43 - 11 predictors"
# [1] "WABcomp Deg - Finpred: 0.63 - Maxcor: 0.54 - 11 predictors"
# [1] "WABcomp Bwn - Finpred: 0.64 - Maxcor: 0.56 - 11 predictors"
# [1] "WABcomp LoT - Finpred: 0.54 - Maxcor: 0.55 - 11 predictors"
# [1] "WABcomp Eff - Finpred: 0.43 - Maxcor: 0.48 - 11 predictors"
#####





# RFE of 72 damaged
#####
rfdata = groupdmg
subsets = c(1:30,seq(33,72, by=3))
rfefiles=list()

for (b in 1:length(behavs) ) {
  behav=eval(parse(text=paste0('SD$',behavs[b])))
  inTrain = (1:53)
  
  rafaja = selectRFE(x=rfdata[inTrain,],
                     y=behav[inTrain],
                     ctrl = ctrl, 
                     subsets = subsets, 
                     metric = metric, 
                     cores=4, 
                     noprint=T, memory = '1.5G')
  variable = paste0('rfefiles$', behavs[b] )
  eval(parse(text=paste0( variable, ' = rafaja' )))
}
rfefiles.DMG=rfefiles


for (behav in behavs) {
      loadfile = eval(parse(text=paste0('rfefiles.DMG$',behav,'$DATAfile') ))
      eval(parse(text=paste0('load(\'', loadfile, '\')' ) ))
      RFE$results=rferesults(RFE)
      RFE$behav = behav
      save(RFE, file=loadfile)
      print(paste( behav, RFE$rfProfile$optsize,'Selected:', round(RFE$results$corSelected,2),
                   'FULL:', round(RFE$results$corFULL,2),'p', round(RFE$results$p.FULLSelected,3) ))
}
rfefiles.DMG$PNTcorrect = list(DATAfile="/data/jag/dpustina/TEMPRFE/selectRFE.CSZJBUCVNG.Rdata")
rfefiles.DMG$WABAQ = list(DATAfile="/data/jag/dpustina/TEMPRFE/selectRFE.OWYEZGFNAU.Rdata")
rfefiles.DMG$WABrep = list(DATAfile="/data/jag/dpustina/TEMPRFE/selectRFE.QJMHCUVUQA.Rdata")
rfefiles.DMG$WABcomp = list(DATAfile="/data/jag/dpustina/TEMPRFE/selectRFE.EKQGZCQXMQ.Rdata")


#####




# function for all-threshold combination
#####
# this will use raw graph measures to create predictions,
# then combine predictions together, all within the same folds
getResultsAllThresh <- function(behavname, dtirest, FOLDS, mes, densities, SD) {
  
  if (dtirest=='DTI') rfefiles=rfefiles.WEIGHTED.DTI
  if (dtirest=='REST') rfefiles=rfefiles.WEIGHTED.REST

  denslist = list()
  eval(parse(text=paste0('denslist$', behavname, '=list()')))

  # first work on raw values and get a final prediction from all densities
  denspreds = matrix(NA, nrow=nrow(SD), ncol=length(densities))
  colnames(denspreds) = paste0('Dens_', densities)
  for (d in 1:length(densities)) {
    dens=densities[d]
    loadfile = eval(parse(text=paste0('rfefiles$',mes,'$',behavname,'$Dens',dens,'$DATAfile') ))
    eval(parse(text=paste0('load(\'', loadfile, '\')' ) ))
    
    inputvars = RFE$x[ , RFE$rfProfile$optVariables]
    colnames(inputvars) = paste0('D', dens, '_', colnames(inputvars))
    resp = eval(parse(text=paste0('SD$',behavname)))
    thispred = rep(NA,length(resp))
    dfram = data.frame(resp=resp, inputvars)
    for (fold in FOLDS) {
      rfm = randomForest(resp ~ ., data=dfram[-fold,])
      thispred[fold] = predict(rfm, newdata = dfram[fold,])
    }
    denspreds[,d] = thispred
  }
  
  dfram = data.frame(resp=resp, denspreds)
  finpred = rep(NA,nrow(dfram))
  for (fold in FOLDS) {
    rfm = randomForest(resp ~ ., data=dfram[-fold,])
    finpred[fold] = predict(rfm, newdata = dfram[fold,])
  }
  return(finpred)
}
#####




###########
# HBM REVISION START HERE
###########



# Predictions for all behavs with all thresholds
#################################
source('/data/jag/dpustina/Code/APHASIA/translable.R')
load('/data/jag/dpustina/restAnalysisSessionJun2.Rdata') # added this to make it work for HBM revision
saveit=F

for (i in 1:length(behavs)) {
  cat(paste('\n', format(Sys.time(), "%H:%M") , "Running", behavs[i],"... \n"))
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
                     LesionSize=les.size #,
#                      MPOWAB=SD$MPOWAB,
#                      Age=SD$AgeatMRI
                     )

  finaldata=data.frame(cor=double(),name=factor())
  finaldata.rmse=data.frame(rmse=double(),name=factor())
  finaldata.preds2 = list()
  foldlist = list()
  rm(rfefinal)
  for (reploop in 1:20) {
    folds = createFolds(y=behav,k=10)
    foldlist[[reploop]] = folds
  
    if (saveit) writeVLSM(SD, behavs[i], paste0(behavs[i],'_VLSM.nii.gz') )
  
    # damage
    loadfile = eval(parse(text=paste0('rfefiles.DMG$',behavs[i],'$DATAfile') ))
    eval(parse(text=paste0('load(\'', loadfile, '\')' ) ))
    preds$DMG = getResults(RFE, folds)
    if (saveit) writeRFElabs(RFE, paste0(behavs[i],'_DMG.txt') )
  
    # DTI connectome
    loadfile = eval(parse(text=paste0('groupselDTI$',behavs[i]) ))
    eval(parse(text=paste0('load(\'', loadfile, '\')' ) ))
    preds$DTImat = getResults(RFE, folds)
    if (saveit) writeRFEcons(RFE, paste0(behavs[i],'_DTIconnect.edge') )
    
    # REST connectome
    loadfile = eval(parse(text=paste0('groupselREST$',behavs[i]) ))
    eval(parse(text=paste0('load(\'', loadfile, '\')' ) ))
    preds$RESTmat = getResults(RFE, folds)
    if (saveit) writeRFEcons(RFE, paste0(behavs[i],'_RESTconnect.edge') )
    
    #DTI deg
    dtirest='DTI'; mes='Deg'
    preds$DTIdeg = getResultsAllThresh(behavname=behavs[i], dtirest=dtirest, FOLDS=folds, mes=mes,
                                       densities=densities, SD=SD)
    
    #DTI bwn
    dtirest='DTI'; mes='Bwn'
    preds$DTIbwn = getResultsAllThresh(behavname=behavs[i], dtirest=dtirest, FOLDS=folds, mes=mes,
                                       densities=densities, SD=SD)
    
    #DTI lot
    dtirest='DTI'; mes='LoT'
    preds$DTIlot = getResultsAllThresh(behavname=behavs[i], dtirest=dtirest, FOLDS=folds, mes=mes,
                                       densities=densities, SD=SD)
    
    #DTI eff
    dtirest='DTI'; mes='Eff'
    preds$DTIeff = getResultsAllThresh(behavname=behavs[i], dtirest=dtirest, FOLDS=folds, mes=mes,
                                       densities=densities, SD=SD)
    
    #REST deg
    dtirest='REST'; mes='Deg'
    preds$RESTdeg = getResultsAllThresh(behavname=behavs[i], dtirest=dtirest, FOLDS=folds, mes=mes,
                                        densities=densities, SD=SD)
    
    #REST bwn
    dtirest='REST'; mes='Bwn'
    preds$RESTbwn = getResultsAllThresh(behavname=behavs[i], dtirest=dtirest, FOLDS=folds, mes=mes,
                                        densities=densities, SD=SD)
  
    #REST lot
    dtirest='REST'; mes='LoT'
    preds$RESTlot = getResultsAllThresh(behavname=behavs[i], dtirest=dtirest, FOLDS=folds, mes=mes,
                                        densities=densities, SD=SD)
    
    #REST eff
    dtirest='REST'; mes='Eff'
    preds$RESTeff = getResultsAllThresh(behavname=behavs[i], dtirest=dtirest, FOLDS=folds, mes=mes,
                                        densities=densities, SD=SD)
    
    finpred = rep(NA, nrow(preds))
    for (fold in folds) {
      rfm = randomForest(x=preds[-fold,], y=behav[-fold], ntree=500)
      finpred[fold] = predict(rfm,preds[fold,])
    }
  
    # variable seleciton from the 11 multimodal predictors
    if (!exists('rfefinal')) {
      rfefinal = rfe(x = preds, y=behav,metric = metric,
                   rfeControl = ctrl,
                   subsets=2:(ncol(preds)-1) )
    }
  
    # prediction from the best selected multimodals
    finpredRFE = rep(NA, nrow(preds))
    for (fold in folds) {
      rfm = randomForest(x=preds[-fold,rfefinal$optVariables], y=behav[-fold], ntree=500)
      finpredRFE[fold] = predict(rfm,preds[fold,rfefinal$optVariables])
    }
    
    
    # new entry, required by reviewer, prediction only with
    # DTImat and RESTmat
    finpredMAT = rep(NA, nrow(preds))
    thesecols = c('DTImat', 'RESTmat')
    for (fold in folds) {
      rfm = randomForest(x=preds[-fold,thesecols], y=behav[-fold], ntree=500)
      finpredMAT[fold] = predict(rfm,preds[fold,thesecols])
    }
    
    
    preds2 = data.frame(preds, FinalAll=finpred, FinalRFE=finpredRFE, FinalMAT=finpredMAT)
    thisresults = abs(cor(preds2,behav))
    thisresults = data.frame(cor=thisresults, name=rownames(thisresults))
    finaldata=rbind(finaldata,thisresults)
    thisresults.rmse = cbind(apply(preds2, 2, function(x) rmse(x,behav)))
    thisresults.rmse = data.frame(rmse=thisresults.rmse, name=rownames(thisresults.rmse))
    finaldata.rmse=rbind(finaldata.rmse,thisresults.rmse)
    finaldata.preds2[[reploop]] = preds2
    
    cat(paste(reploop, ''))
  } # repetion loop


  finaldata$name = gsub('DTI','DTI_', finaldata$name)
  finaldata$name = gsub('REST','REST_', finaldata$name)
  finaldata$name = gsub('DMG','Parc_Damage', finaldata$name)
  finaldata$name = gsub('LesionSize','Les_Size', finaldata$name)
  finaldata$name = gsub('Final','Final_', finaldata$name)

  levorder=c(
             'DTI_deg', 'DTI_bwn', 'DTI_lot', 'DTI_eff','DTI_mat', 
             'REST_deg', 'REST_bwn', 'REST_lot', 'REST_eff','REST_mat',
             'Les_Size', 'Parc_Damage',
             'Final_All', 'Final_RFE', 'Final_MAT')
  finaldata$name = factor(finaldata$name, levels=levorder)

  fills = rep('black', ncol(preds2))
  fills[grep("Final_All",levorder)] = 'red'
  fills[grep("Final_RFE",levorder)] = 'red4'
  fills[grep("Final_MAT",levorder)] = 'white'
  fills[grep("Parc_",levorder)] = 'orange4'
  fills[grep("Les_",levorder)] = 'orange4'
  fills[grep("DTI_",levorder)] = 'tomato'
  fills[grep("REST_",levorder)] = 'royalblue'
  axcol = fills
  axcol[is.na(axcol)] = 'black'
  #axcol[-grep("Final_",levorder)] = 'black'
  selvars.x=rfefinal$optVariables
  #selvars.y=rep(max(finaldata$cor),length(selvars.x))
  selvars.x = gsub('DTI','DTI_', selvars.x)
  selvars.x = gsub('REST','REST_', selvars.x)
  selvars.x = gsub('DMG','Parc_Damage', selvars.x)
  selvars.x = gsub('LesionSize','Les_Size', selvars.x)

  ggplot(finaldata, aes(factor(name), cor, fill=name)) + 
    geom_boxplot() + # ylim(0.2,0.95) +
    scale_y_continuous(breaks=seq(0.2,0.9,0.1), limits=c(0.2,0.95), minor_breaks = NULL) +
    ylab('Correlation') +
    ggtitle(behavs[i]) +
    scale_fill_manual('', values=fills) +
    guides(fill=F) +
    annotate('text',x=selvars.x, y=rep(0.95, length(selvars.x)),
             label='*', color='gray36', size=12) +
    annotate('segment', x=12.5, xend=12.5, y=0.2,yend=0.95, color='red',alpha=0.4,lwd=1.5,lty=1) + 
    annotate('segment', x=14.5, xend=14.5, y=0.2,yend=0.95, color='red4',alpha=0.4,lwd=1.5,lty=1) + # annotate("rect", xmin=12.5, xmax=14.5, ymin=0.2, ymax=0.95, alpha=0.2) +
    theme(
      axis.text.x=element_text(angle=45,hjust=1,size=16,color=axcol),
      axis.text.y=element_text(size=14),
      axis.title.x=element_blank(),
      axis.title.y=element_text(size=18),
      plot.title=element_text(size=20))
  ggsave(width=8.5, height=6.5, 
         filename=file.path('/data/jag/dpustina/APHASIA/STROKE/analyses/papier',paste0(behavs[i],'_SMAPplot_allThresh_HBMrevision.png') ))
  
# this contains the correlations from all 10 runs
  save(finaldata, rfefinal, finaldata.rmse, finaldata.preds2, foldlist,
       file=paste0('/data/jag/dpustina/APHASIA/STROKE/analyses/papier/',behavs[i],'_finaldata_HBMrevision.Rdata'))

  colnames(preds2) = gsub('DTI','DTI_', colnames(preds2))
  colnames(preds2) = gsub('REST','REST_', colnames(preds2))
  colnames(preds2) = gsub('DMG','Parc_Damage', colnames(preds2))
  colnames(preds2) = gsub('LesionSize','Les_Size', colnames(preds2))
  colnames(preds2) = gsub('Final','Final_', colnames(preds2))
  for (me in 1:ncol(preds2)) {
    png(filename = file.path('/data/jag/dpustina/APHASIA/STROKE/analyses/papier',
                             paste0(behavs[i],'_', colnames(preds2)[me],'_scatter_AllThresh_HBMrevision.png')),
        width=600, height=550,pointsize=20)
    plot(behav,preds2[,me], ylim=range(behav),
         pch=19, xlab='',ylab='',
         main=paste(behavs[i],'-',colnames(preds2)[me]) )
    mtext(paste0(
      'r=', round(cor(preds2[,me],behav),2),' ',
      'RMSE=', round(rmse(preds2[me],behav),2)), side=3)
    abline(a=0,b=1)
    
    dev.off()
    Sys.sleep(2)
  }
}
#############################

# 
# #####################################
# # save SMAP plots one more time with fixed FinalMAT
# for (i in 1:length(behavs)) {
#   
#   cat(paste('\n', format(Sys.time(), "%H:%M") , "Running", behavs[i],"... \n"))
#   behav=eval(parse(text=paste0('SD$',behavs[i])))
#   
#   load(paste0('/data/jag/dpustina/APHASIA/STROKE/analyses/papier/',behavs[i],'_finaldata_HBMrevision.Rdata'))
#   
#   finaldata=data.frame(cor=double(),name=factor())
#   for (reploop in 1:length(finaldata.preds2)) {
#     preds2 = finaldata.preds2[[reploop]]
#     thisresults = abs(cor(preds2,behav))
#     thisresults = data.frame(cor=thisresults, name=rownames(thisresults))
#     finaldata=rbind(finaldata,thisresults)
#   }
#    
#   
#   finaldata$name = gsub('DTI','DTI_', finaldata$name)
#   finaldata$name = gsub('REST','REST_', finaldata$name)
#   finaldata$name = gsub('DMG','Parc_Damage', finaldata$name)
#   finaldata$name = gsub('LesionSize','Les_Size', finaldata$name)
#   finaldata$name = gsub('Final','Final_', finaldata$name)
#   
#   levorder=c(
#     'DTI_deg', 'DTI_bwn', 'DTI_lot', 'DTI_eff','DTI_mat', 
#     'REST_deg', 'REST_bwn', 'REST_lot', 'REST_eff','REST_mat',
#     'Les_Size', 'Parc_Damage',
#     'Final_All', 'Final_RFE', 'Final_MAT')
#   finaldata$name = factor(finaldata$name, levels=levorder)
#   
#   
#   fills = rep('black', ncol(preds2))
#   fills[grep("Final_All",levorder)] = 'red'
#   fills[grep("Final_RFE",levorder)] = 'red4'
#   fills[grep("Final_MAT",levorder)] = NA
#   fills[grep("Parc_",levorder)] = 'orange4'
#   fills[grep("Les_",levorder)] = 'orange4'
#   fills[grep("DTI_",levorder)] = 'tomato'
#   fills[grep("REST_",levorder)] = 'royalblue'
#   axcol = fills
#   axcol[is.na(axcol)] = 'black'
#   #axcol[-grep("Final_",levorder)] = 'black'
#   selvars.x=rfefinal$optVariables
#   #selvars.y=rep(max(finaldata$cor),length(selvars.x))
#   selvars.x = gsub('DTI','DTI_', selvars.x)
#   selvars.x = gsub('REST','REST_', selvars.x)
#   selvars.x = gsub('DMG','Parc_Damage', selvars.x)
#   selvars.x = gsub('LesionSize','Les_Size', selvars.x)
#   
#   ggplot(finaldata, aes(factor(name), cor, fill=name)) + 
#     geom_boxplot() + # ylim(0.2,0.95) +
#     scale_y_continuous(breaks=seq(0.2,0.9,0.1), limits=c(0.2,0.95), minor_breaks = NULL) +
#     ylab('Correlation') +
#     ggtitle(behavs[i]) +
#     scale_fill_manual('', values=fills) +
#     guides(fill=F) +
#     annotate('text',x=selvars.x, y=rep(0.95, length(selvars.x)),
#              label='*', color='gray36', size=12) +
#     annotate('segment', x=12.5, xend=12.5, y=0.2,yend=0.95, color='red',alpha=0.4,lwd=1.5,lty=1) + 
#     annotate('segment', x=14.5, xend=14.5, y=0.2,yend=0.95, color='red4',alpha=0.4,lwd=1.5,lty=1) + # annotate("rect", xmin=12.5, xmax=14.5, ymin=0.2, ymax=0.95, alpha=0.2) +
#     theme(
#       axis.text.x=element_text(angle=45,hjust=1,size=16,color=axcol),
#       axis.text.y=element_text(size=14),
#       axis.title.x=element_blank(),
#       axis.title.y=element_text(size=18),
#       plot.title=element_text(size=20))
#   
#   ggsave(width=8.5, height=6.5, 
#          filename=file.path('/data/jag/dpustina/APHASIA/STROKE/analyses/papier',paste0(behavs[i],'_SMAPplot_allThresh_HBMrevision.png') ))
#   
#   save(finaldata, rfefinal, finaldata.rmse, finaldata.preds2, foldlist,
#        file=paste0('/data/jag/dpustina/APHASIA/STROKE/analyses/papier/',behavs[i],'_finaldata_HBMrevision.Rdata'))
#   
# }
# ######################################






# get some results from finaldata for the paper
##################
for (i in 1:length(behavs)) {
  load(paste0('/data/jag/dpustina/APHASIA/STROKE/analyses/papier/',behavs[i],'_finaldata_HBMrevision.Rdata'))
  behav=eval(parse(text=paste0('SD$',behavs[i])))
  
  checkvars = unique(finaldata$name)
  for (var in checkvars) {
    cat(
      paste0(behavs[i], ' (range ', paste0(range(behav),collapse='-') , '): ', var,
            ' Corr=', round(mean(finaldata$cor[finaldata$name==var]),2), ' (', round(sd(finaldata$cor[finaldata$name==var]),2), ')',
            ' RMSE=', round(mean(finaldata.rmse$rmse[finaldata$name==var]),2), ' (', round(sd(finaldata.rmse$rmse[finaldata$name==var]),2), ')',
            '\n' )
      )
  }
}

# PNTcorrect (range 1-98): DTI_deg Corr=0.56 (0.04) RMSE=24.12 (0.86)
# PNTcorrect (range 1-98): DTI_bwn Corr=0.6 (0.03) RMSE=23.16 (0.65)
# PNTcorrect (range 1-98): DTI_lot Corr=0.78 (0.02) RMSE=17.94 (0.63)
# PNTcorrect (range 1-98): DTI_eff Corr=0.6 (0.05) RMSE=23.22 (1.18)
# PNTcorrect (range 1-98): DTI_mat Corr=0.71 (0.02) RMSE=20.52 (0.56)
# PNTcorrect (range 1-98): REST_deg Corr=0.46 (0.07) RMSE=25.64 (1.09)
# PNTcorrect (range 1-98): REST_bwn Corr=0.72 (0.04) RMSE=20.75 (0.86)
# PNTcorrect (range 1-98): REST_lot Corr=0.44 (0.08) RMSE=26.07 (1.3)
# PNTcorrect (range 1-98): REST_eff Corr=0.57 (0.08) RMSE=23.64 (1.44)
# PNTcorrect (range 1-98): REST_mat Corr=0.72 (0.04) RMSE=23.95 (0.46)
# PNTcorrect (range 1-98): Parc_Damage Corr=0.58 (0.02) RMSE=23.56 (0.47)
# PNTcorrect (range 1-98): Les_Size Corr=0.47 (0) RMSE=81901.4 (0)
# PNTcorrect (range 1-98): Final_All Corr=0.82 (0.03) RMSE=16.57 (1.19)
# PNTcorrect (range 1-98): Final_RFE Corr=0.85 (0.03) RMSE=15.61 (1.2)
# WABAQ (range 30.4-97.6): DTI_deg Corr=0.59 (0.03) RMSE=15.08 (0.46)
# WABAQ (range 30.4-97.6): DTI_bwn Corr=0.69 (0.04) RMSE=13.46 (0.66)
# WABAQ (range 30.4-97.6): DTI_lot Corr=0.71 (0.02) RMSE=13.23 (0.47)
# WABAQ (range 30.4-97.6): DTI_eff Corr=0.71 (0.02) RMSE=13.23 (0.43)
# WABAQ (range 30.4-97.6): DTI_mat Corr=0.7 (0.02) RMSE=13.5 (0.27)
# WABAQ (range 30.4-97.6): REST_deg Corr=0.64 (0.04) RMSE=14.55 (0.54)
# WABAQ (range 30.4-97.6): REST_bwn Corr=0.81 (0.03) RMSE=11.44 (0.64)
# WABAQ (range 30.4-97.6): REST_lot Corr=0.65 (0.04) RMSE=14.31 (0.51)
# WABAQ (range 30.4-97.6): REST_eff Corr=0.35 (0.07) RMSE=17.6 (0.61)
# WABAQ (range 30.4-97.6): REST_mat Corr=0.74 (0.02) RMSE=13.75 (0.27)
# WABAQ (range 30.4-97.6): Parc_Damage Corr=0.69 (0.02) RMSE=13.63 (0.33)
# WABAQ (range 30.4-97.6): Les_Size Corr=0.61 (0) RMSE=81891.46 (0)
# WABAQ (range 30.4-97.6): Final_All Corr=0.88 (0.02) RMSE=9.45 (0.46)
# WABAQ (range 30.4-97.6): Final_RFE Corr=0.89 (0.02) RMSE=9.02 (0.47)
# WABrep (range 1.1-10): DTI_deg Corr=0.6 (0.08) RMSE=1.84 (0.15)
# WABrep (range 1.1-10): DTI_bwn Corr=0.71 (0.02) RMSE=1.62 (0.05)
# WABrep (range 1.1-10): DTI_lot Corr=0.62 (0.03) RMSE=1.8 (0.06)
# WABrep (range 1.1-10): DTI_eff Corr=0.71 (0.03) RMSE=1.63 (0.06)
# WABrep (range 1.1-10): DTI_mat Corr=0.77 (0.01) RMSE=1.5 (0.03)
# WABrep (range 1.1-10): REST_deg Corr=0.72 (0.03) RMSE=1.65 (0.07)
# WABrep (range 1.1-10): REST_bwn Corr=0.86 (0.02) RMSE=1.3 (0.05)
# WABrep (range 1.1-10): REST_lot Corr=0.66 (0.03) RMSE=1.74 (0.06)
# WABrep (range 1.1-10): REST_eff Corr=0.43 (0.1) RMSE=2.09 (0.14)
# WABrep (range 1.1-10): REST_mat Corr=0.72 (0.02) RMSE=1.64 (0.03)
# WABrep (range 1.1-10): Parc_Damage Corr=0.64 (0.02) RMSE=1.78 (0.04)
# WABrep (range 1.1-10): Les_Size Corr=0.51 (0) RMSE=81936.37 (0)
# WABrep (range 1.1-10): Final_All Corr=0.87 (0.02) RMSE=1.16 (0.06)
# WABrep (range 1.1-10): Final_RFE Corr=0.87 (0.02) RMSE=1.13 (0.06)
# WABcomp (range 3.75-10): DTI_deg Corr=0.54 (0.05) RMSE=1.26 (0.06)
# WABcomp (range 3.75-10): DTI_bwn Corr=0.73 (0.05) RMSE=1.02 (0.07)
# WABcomp (range 3.75-10): DTI_lot Corr=0.59 (0.06) RMSE=1.21 (0.07)
# WABcomp (range 3.75-10): DTI_eff Corr=0.69 (0.05) RMSE=1.08 (0.07)
# WABcomp (range 3.75-10): DTI_mat Corr=0.62 (0.02) RMSE=1.17 (0.02)
# WABcomp (range 3.75-10): REST_deg Corr=0.6 (0.04) RMSE=1.2 (0.05)
# WABcomp (range 3.75-10): REST_bwn Corr=0.62 (0.07) RMSE=1.18 (0.07)
# WABcomp (range 3.75-10): REST_lot Corr=0.52 (0.05) RMSE=1.27 (0.05)
# WABcomp (range 3.75-10): REST_eff Corr=0.41 (0.06) RMSE=1.36 (0.05)
# WABcomp (range 3.75-10): REST_mat Corr=0.5 (0.09) RMSE=1.38 (0.03)
# WABcomp (range 3.75-10): Parc_Damage Corr=0.55 (0.04) RMSE=1.25 (0.04)
# WABcomp (range 3.75-10): Les_Size Corr=0.54 (0) RMSE=81935.11 (0)
# WABcomp (range 3.75-10): Final_All Corr=0.79 (0.03) RMSE=0.93 (0.06)
# WABcomp (range 3.75-10): Final_RFE Corr=0.78 (0.03) RMSE=0.95 (0.05)
##################################################




# produce results for BrainNet plots and stats from all thresholds
############################################
for (behavname in behavs) {
  for (dtirest in c('DTI', 'REST')) {
    for (mes in c('Deg', 'Bwn', 'LoT', 'Eff')) {
      allThreshNiftii(behavname, dtirest,  mes, densities, SD)
    }
  }
}

allThreshNiftii <- function(behavname, dtirest,  mes, densities, SD) {
  if (dtirest=='DTI') rfefiles=rfefiles.WEIGHTED.DTI
  if (dtirest=='REST') rfefiles=rfefiles.WEIGHTED.REST
  
  labvals = rep(0,268)
  nvarsdens = rep(0,length(densities))
  for (d in 1:length(densities)) {
    dens=densities[d]
    loadfile = eval(parse(text=paste0('rfefiles$',mes,'$',behavname,'$Dens',dens,'$DATAfile') ))
    eval(parse(text=paste0('load(\'', loadfile, '\')' ) ))
    vars = RFE$rfProfile$optVariables
    vars = as.numeric(gsub("[^\\d]+", "", vars, perl=TRUE))
    labvals[vars] = labvals[vars]+1
    nvarsdens[d] = length(vars)
  }
  cat(paste0('\n', behavname, ' ', dtirest,tolower(mes), ' has ', round(mean(nvarsdens),1), ' (+/-', round(sd(nvarsdens),0), ') average variables, and ', sum(labvals>0), ' regions overall.'))
  
  fname = paste0('/data/jag/dpustina/APHASIA/STROKE/analyses/papier/', behavname,'_',dtirest,tolower(mes),'_MNI_allThresh.nii.gz' )
  translabelMNI(labvals, file = fname) # translabelmerge for ITKsnap view
}

# PNTcorrect DTIdeg has 11.5 (+/-7) average variables, and 62 regions overall.
# PNTcorrect DTIbwn has 22.7 (+/-11) average variables, and 118 regions overall.
# PNTcorrect DTIlot has 4.5 (+/-1) average variables, and 23 regions overall.
# PNTcorrect DTIeff has 14.7 (+/-10) average variables, and 90 regions overall.
# PNTcorrect RESTdeg has 60.2 (+/-38) average variables, and 209 regions overall.
# PNTcorrect RESTbwn has 21.7 (+/-21) average variables, and 149 regions overall.
# PNTcorrect RESTlot has 19.9 (+/-18) average variables, and 123 regions overall.
# PNTcorrect RESTeff has 20.3 (+/-21) average variables, and 138 regions overall.
# WABAQ DTIdeg has 8.4 (+/-11) average variables, and 64 regions overall.
# WABAQ DTIbwn has 21.1 (+/-11) average variables, and 93 regions overall.
# WABAQ DTIlot has 17.4 (+/-8) average variables, and 80 regions overall.
# WABAQ DTIeff has 6.6 (+/-4) average variables, and 44 regions overall.
# WABAQ RESTdeg has 34 (+/-30) average variables, and 153 regions overall.
# WABAQ RESTbwn has 22.8 (+/-20) average variables, and 149 regions overall.
# WABAQ RESTlot has 17.3 (+/-16) average variables, and 116 regions overall.
# WABAQ RESTeff has 30.5 (+/-24) average variables, and 186 regions overall.
# WABrep DTIdeg has 4.2 (+/-1) average variables, and 12 regions overall.
# WABrep DTIbwn has 5.8 (+/-2) average variables, and 36 regions overall.
# WABrep DTIlot has 16.3 (+/-10) average variables, and 60 regions overall.
# WABrep DTIeff has 9.1 (+/-11) average variables, and 67 regions overall.
# WABrep RESTdeg has 47.4 (+/-47) average variables, and 190 regions overall.
# WABrep RESTbwn has 12.1 (+/-22) average variables, and 105 regions overall.
# WABrep RESTlot has 14.5 (+/-15) average variables, and 87 regions overall.
# WABrep RESTeff has 15.9 (+/-14) average variables, and 118 regions overall.
# WABcomp DTIdeg has 15.5 (+/-14) average variables, and 92 regions overall.
# WABcomp DTIbwn has 10.3 (+/-6) average variables, and 49 regions overall.
# WABcomp DTIlot has 20.3 (+/-15) average variables, and 96 regions overall.
# WABcomp DTIeff has 3.2 (+/-3) average variables, and 21 regions overall.
# WABcomp RESTdeg has 13.6 (+/-18) average variables, and 88 regions overall.
# WABcomp RESTbwn has 37.9 (+/-27) average variables, and 208 regions overall.
# WABcomp RESTlot has 31.9 (+/-27) average variables, and 170 regions overall.
# WABcomp RESTeff has 12.9 (+/-19) average variables, and 110 regions overall.
#######################################################





getResults <- function(RFE, folds) {
  require(randomForest)
  output = rep(NA, nrow(RFE$x))
  for (fold in folds) {
    rfm = randomForest(y=RFE$y[-fold], x=RFE$x[ -fold , RFE$rfProfile$optVariables], ntree=500)
    output[fold] = predict(rfm, newdata=RFE$x[ fold, RFE$rfProfile$optVariables])
  }
  return(output)
}

temp=rfe(x = preds, y=behav,sizes = 2:11,metric = 'RMSE',rfeControl = ctrl,maximize = F)
plot(temp, type=c('g','o'), main=paste(temp$optsize, 'optimal variables') )




# fisher Z of Final_All vs. max correlation
#############################################
library(psych)

i=1
behav=eval(parse(text=paste0('SD$',behavs[i])))
load(paste0('/data/jag/dpustina/APHASIA/STROKE/analyses/papier/',behavs[i],'_finaldata_HBMrevision.Rdata'))

corrs = data.frame(cormax = rep(NA,20), corfin = rep(NA,20))
better = 0

bestvar = rep(NA,20)
for (jj in 1:20) bestvar[jj] = which.max( cor(behav, finaldata.preds2[[jj]])[1:11])
ux <- unique(bestvar)
maxcor = ux[which.max(tabulate(match(bestvar, ux)))]
fincor = 13

for (jj in 1:20) {

  
  xy = cor(behav,finaldata.preds2[[jj]][,maxcor])
  xz = cor(behav,finaldata.preds2[[jj]][,fincor])
  yz = cor(finaldata.preds2[[jj]][,fincor],finaldata.preds2[[jj]][,maxcor])
  nsub = nrow(finaldata.preds2[[jj]])
  pdiff = paired.r(xy = xy, xz = xz, yz = yz, n=nsub)$p
  if (pdiff < 0.05) better=better+1
  corrs$cormax[jj] = xy
  corrs$corfin[jj] = xz
  
  noooo = ''
  if (xy > xz) noooo ='<--------------'
  
  print(paste0(behavs[i], ' Final_All (', round(xz,2), ') vs ', maxcor,
              colnames(finaldata.preds2[[jj]])[maxcor], ' (', round(xy,2), '): ', round(pdiff,3), noooo
          )
        )
  
}

corrs2=fisherz(corrs)
runsdiff = t.test(x = corrs2$corfin, y = corrs2$cormax, paired = T, alternative="greater")$p.value
print(paste(behavs[i], 'Overall, Final > Max in', round(better/20*100,0), '% of times, and fisherz yields p =', runsdiff))

