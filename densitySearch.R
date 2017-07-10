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
plot(RFE$rfProfile, type=c('g','o'), main=paste('Optimal -', RFE$rfProfile$optsize, 'variables') )


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




# preds at select densities
#####
seldens = list()
dumdum=matrix(1,nrow=1,ncol=4); colnames(dumdum) = checks
seldens$PNTcorrect = c(0.4,0.4,0.3,0.5)*dumdum
seldens$WABAQ = c(0.1,0.1,0.3,0.1)*dumdum
seldens$WABrep = c(0.7,0.4,0.1,0.7)*dumdum
seldens$WABcomp = c(0.1,0.6,0.3,0.1)*dumdum
seldensDTI = seldens
seldens = list()
dumdum=matrix(1,nrow=1,ncol=4); colnames(dumdum) = checks
seldens$PNTcorrect = c(0.6,0.4,0.8,1)*dumdum
seldens$WABAQ = c(0.9,1,1,1)*dumdum
seldens$WABrep = c(0.9,0.3,0.5,1)*dumdum
seldens$WABcomp = c(0.1,0.2,0.4,1)*dumdum
seldensREST=seldens

# groupmat selections
groupsel = list()
groupsel$PNTcorrect=rfelist.ALL.reranksimple.PNT[[1]]$DATAfile
groupsel$WABAQ=rfelist.ALL.reranksimple.WABAQ[[1]]$DATAfile
groupsel$WABrep=rfelist.ALL.reranksimple.WABrep[[1]]$DATAfile
groupsel$WABcomp=rfelist.ALL.reranksimple.WABcomp[[1]]$DATAfile
groupselDTI=groupsel
groupsel = list()
groupsel$PNTcorrect=rfelist.ALLbold.reranksimple.PNT[[1]]$DATAfile
groupsel$WABAQ=rfelist.ALLbold.reranksimple.WABAQ[[1]]$DATAfile
groupsel$WABrep=rfelist.ALLbold.reranksimple.WABrep[[1]]$DATAfile
groupsel$WABcomp=rfelist.ALLbold.reranksimple.WABcomp[[1]]$DATAfile
groupselREST=groupsel

source('/data/jag/dpustina/Code/APHASIA/translable.R')
saveit=F

for (i in 1:length(behavs)) {
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
  rm(rfefinal)
  for (reploop in 1:10) {
    folds = createFolds(y=behav,k=10)
  
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
    dens = eval(parse(text=paste0('seldensDTI$', behavs[i], '[ ,\'Deg\']') ))
    loadfile = eval(parse(text=paste0('rfefiles.WEIGHTED.DTI$Deg$', behavs[i], '$Dens', dens, '$DATAfile') ))
    eval(parse(text=paste0('load(\'', loadfile, '\')' ) ))
    preds$DTIdeg = getResults(RFE, folds)
    if (saveit) writeRFElabs(RFE, what=paste0(behavs[i],'_DTIdeg.txt') )
    
    #DTI bwn
    dens = eval(parse(text=paste0('seldensDTI$', behavs[i], '[ ,\'Bwn\']') ))
    loadfile = eval(parse(text=paste0('rfefiles.WEIGHTED.DTI$Bwn$', behavs[i], '$Dens', dens, '$DATAfile') ))
    eval(parse(text=paste0('load(\'', loadfile, '\')' ) ))
    preds$DTIbwn = getResults(RFE, folds)
    if (saveit) writeRFElabs(RFE, paste0(behavs[i],'_DTIbwn.txt') )
    
    #DTI lot
    dens = eval(parse(text=paste0('seldensDTI$', behavs[i], '[ ,\'LoT\']') ))
    loadfile = eval(parse(text=paste0('rfefiles.WEIGHTED.DTI$LoT$', behavs[i], '$Dens', dens, '$DATAfile') ))
    eval(parse(text=paste0('load(\'', loadfile, '\')' ) ))
    preds$DTIlot = getResults(RFE, folds)
    if (saveit) writeRFElabs(RFE, paste0(behavs[i],'_DTIlot.txt') )
    
    #DTI eff
    dens = eval(parse(text=paste0('seldensDTI$', behavs[i], '[ ,\'Eff\']') ))
    loadfile = eval(parse(text=paste0('rfefiles.WEIGHTED.DTI$Eff$', behavs[i], '$Dens', dens, '$DATAfile') ))
    eval(parse(text=paste0('load(\'', loadfile, '\')' ) ))
    preds$DTIeff = getResults(RFE, folds)
    if (saveit) writeRFElabs(RFE, paste0(behavs[i],'_DTIeff.txt') )
    
    #REST deg
    dens = eval(parse(text=paste0('seldensREST$', behavs[i], '[ ,\'Deg\']') ))
    loadfile = eval(parse(text=paste0('rfefiles.WEIGHTED.REST$Deg$', behavs[i], '$Dens', dens, '$DATAfile') ))
    eval(parse(text=paste0('load(\'', loadfile, '\')' ) ))
    preds$RESTdeg = getResults(RFE, folds)
    if (saveit) writeRFElabs(RFE, paste0(behavs[i],'_RESTdeg.txt') )
    
    #REST bwn
    dens = eval(parse(text=paste0('seldensREST$', behavs[i], '[ ,\'Bwn\']') ))
    loadfile = eval(parse(text=paste0('rfefiles.WEIGHTED.REST$Bwn$', behavs[i], '$Dens', dens, '$DATAfile') ))
    eval(parse(text=paste0('load(\'', loadfile, '\')' ) ))
    preds$RESTbwn = getResults(RFE, folds)
    if (saveit) writeRFElabs(RFE, paste0(behavs[i],'_RESTbwn.txt') )
  
    #REST lot
    dens = eval(parse(text=paste0('seldensREST$', behavs[i], '[ ,\'LoT\']') ))
    loadfile = eval(parse(text=paste0('rfefiles.WEIGHTED.REST$LoT$', behavs[i], '$Dens', dens, '$DATAfile') ))
    eval(parse(text=paste0('load(\'', loadfile, '\')' ) ))
    preds$RESTlot = getResults(RFE, folds)
    if (saveit) writeRFElabs(RFE, paste0(behavs[i],'_RESTlot.txt') )
    
    #REST eff
    dens = eval(parse(text=paste0('seldensREST$', behavs[i], '[ ,\'Eff\']') ))
    loadfile = eval(parse(text=paste0('rfefiles.WEIGHTED.REST$Eff$', behavs[i], '$Dens', dens, '$DATAfile') ))
    eval(parse(text=paste0('load(\'', loadfile, '\')' ) ))
    preds$RESTeff = getResults(RFE, folds)
    if (saveit) writeRFElabs(RFE, paste0(behavs[i],'_RESTeff.txt') )
    
    finpred = rep(NA, nrow(preds))
    for (fold in folds) {
      rfm = randomForest(x=preds[-fold,], y=behav[-fold], ntree=500)
      finpred[fold] = predict(rfm,preds[fold,])
    }
  
    if (!exists('rfefinal')) {
      rfefinal = rfe(x = preds, y=behav,metric = metric,
                   rfeControl = ctrl,
                   subsets=2:(ncol(preds)-1) )
    }

    finpredRFE = rep(NA, nrow(preds))
    for (fold in folds) {
      rfm = randomForest(x=preds[-fold,rfefinal$optVariables], y=behav[-fold], ntree=500)
      finpredRFE[fold] = predict(rfm,preds[fold,rfefinal$optVariables])
    }
    
    
    preds2 = data.frame(preds, FinalAll=finpred, FinalRFE=finpredRFE)
    thisresults = abs(cor(preds2,behav))
    thisresults = data.frame(cor=thisresults, name=rownames(thisresults))
    finaldata=rbind(finaldata,thisresults)
  } # repetion loop


  levorder=c('LesionSize', 'DMG', 
             'DTIdeg', 'DTIbwn', 'DTIlot', 'DTIeff',
             'RESTdeg', 'RESTbwn', 'RESTlot', 'RESTeff',
             'DTImat', 'RESTmat',
             'FinalAll', 'FinalRFE')
  finaldata$name = factor(finaldata$name, levels=levorder)

  fills = rep('lightpink', ncol(preds2))
  fills[13] = 'red'
  fills[14] = 'red4'
  axcol = fills
  axcol[1:12] = 'black'
  selvars.x=rfefinal$optVariables
  selvars.y=rep(max(finaldata$cor),length(selvars.x))

  ggplot(finaldata, aes(factor(name), cor, fill=name)) + 
    geom_boxplot() +
    ylab('Correlation') +
    ggtitle(behavs[i]) +
    scale_fill_manual('', values=fills) +
    guides(fill=F) +
    annotate('text',x=selvars.x, y=selvars.y,
             label='*', color=fills[14], size=15) +
    theme(
      axis.text.x=element_text(angle=45,hjust=1,size=16,color=axcol),
      axis.text.y=element_text(size=14),
      axis.title.x=element_blank(),
      axis.title.y=element_text(size=18),
      plot.title=element_text(size=20))
  ggsave(width=8, height=6.5, 
         filename=file.path('/data/jag/dpustina/APHASIA/STROKE/analyses/papier',paste0(behavs[i],'_SMAPplot.png') ))
  

  for (me in 1:ncol(preds2)) {
    png(filename = file.path('/data/jag/dpustina/APHASIA/STROKE/analyses/papier',
                             paste0(behavs[i],'_', colnames(preds2)[me],'_scatter.png')),
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



