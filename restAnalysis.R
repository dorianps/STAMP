source('/data/jag/dpustina/Code/APHASIA/loadExcel.R', echo=F)

restppl =  SD$rest1==1 & !is.na(SD$rest1)
SD = SD[ restppl, ]
restnames = SD$NameOnDisk

restfold = '/data/jag/dpustina/APHASIA/STROKE/processing/REST/'

parcel='FIN'
nnodes = ifelse(parcel=='AAL',116,268)
nodefile=ifelse(parcel=='AAL','/AALnetwork116.txt','/FINnetwork268.txt')
# load one for template
temp = matrix(NA, nrow=nnodes, ncol=nnodes)
fcn = length(temp[ upper.tri(temp) ])
for (i in 1:nnodes) for (j in 1:nnodes)
  temp[i,j] = paste0('V',i,'_',j)
groupmat = matrix(NA, nrow=nrow(SD), ncol=fcn)
groupdmg = matrix(NA, nrow=nrow(SD), ncol=nnodes)
colnames(groupmat) = temp[upper.tri(temp)]
submats = list()

# matrix of upper triangles
for (i in 1:nrow(SD)) {
  thisdate = list.files(paste0(restfold,restnames[i],'/'))
  finf = Sys.glob(paste0(restfold,restnames[i],'/',thisdate,nodefile))
  svec = as.matrix(read.table(finf))
  submats[[i]] = svec
  groupmat[ i, ] = svec[ upper.tri(svec) ]
  
  finf = Sys.glob(paste0(restfold,restnames[i],'/',thisdate,'/strokeRest/ANTSLOAD.Rdata'))
  load(finf)
  groupdmg[i,] = ANTSLOAD$rsnetwork$nodeDamage[[ifelse(parcel=='AAL',1,2)]]
  cat(paste(i,''))
}

groupdmg = groupdmg[ , colSums(groupdmg) > 1]



######################################################
# DTI things
#
# matrix of upper triangles
nnodes = 268
dirDTI = '/data/jag/dpustina/APHASIA/STROKE/processing/IITconnectome'
temp = matrix(NA, nrow=nnodes, ncol=nnodes)
fcn = length(temp[ upper.tri(temp) ])
for (i in 1:nnodes) for (j in 1:nnodes)
  temp[i,j] = paste0('V',i,'_',j)
groupmatDTI = matrix(NA, nrow=nrow(SD), ncol=fcn)
colnames(groupmatDTI) = temp[upper.tri(temp)]
submatsDTI = list()
for (i in 1:nrow(SD)) {
  finf = Sys.glob(file.path(dirDTI,SD$NameOnDisk[i], '*NormalizedFinnConnectome.csv'))
  svec = as.matrix(read.table(finf))
  submatsDTI[[i]] = svec
  groupmatDTI[ i, ] = svec[ upper.tri(svec) ]
}

# graph measures
graphmatDTI = list()
for (i in 1:nrow(SD)) {
  temp = makeGraph(abs(submatsDTI[[i]]), graphdensity = 0.25, getEfficiency = T)
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

# this might be good to show why we selected the above measures:
# get quickly pheatmaps for everbody and create a movie, 
# will show stable high correlations for some, and we selected those least correlated
# mm = cor(cbind(temp$centrality, temp$closeness, temp$pagerank, temp$degree, temp$betweeness, temp$localtransitivity, temp$strength, temp$degcent, temp$hubScore, temp$effinv))
# pheatmap(mm, cluster_rows = F, cluster_cols = F, show_rownames = T)

# RFs with no selection, or post selection
cat(paste(format(Sys.time(), "%H:%M") , 'Start RF'))
ntrees=400
sel= grep("Eff", colnames(graphmatpredDTI) )
rfdata = as.matrix(graphmatpredDTI[,sel])[,effSelection$selTopSig]
predgraphDTI = rep(NA, nrow(SD))
require(randomForest)
for (i in 1:nrow(SD)) {
  cat(paste(i," "))
  rfm = randomForest(x = as.matrix(rfdata[ -i, ]),
                     y = SD$WABAQ[-i], 
                     ntree = ntrees)
  predgraphDTI[i] = predict(rfm, newdata = rfdata[ i, ])
}

cor(SD$WABAQ, predgraphDTI)
plot(SD$WABAQ, predgraphDTI)
# r=0.48 for all 805 measures
# r=0.45   for degree
# r=0.47 for betweenness
# r=0.46 for LoTrans
# r=0.25 for Glo Trans



### pRF selection
behav=SD$WABAQ
nruns = 200
nperms=500 # some run with 1000, but maybe just lost time
ntree=500 #500 for graphs, 1000 for raw corrs

source('/data/jag/dpustina/Code/APHASIA/selectVarsRF.R',echo=F)

# graphs
checkfor = 'Eff'
sel= grep(checkfor, colnames(graphmatpredDTI) )
predictors = data.frame(graphmatpredDTI[,sel])
# raw
# predictors = data.frame(groupmatDTI)

selectVarsRF(x = predictors, y=behav,
              ntree=ntree, nperms=nperms,
             nruns=nruns, memory='1G')

degreeSelection = selectVarsRF.result('EBRDMVGYPL')
bwnSelection = selectVarsRF.result('MEVUAOQSFW')
LoTSelection = selectVarsRF.result('HOTZDPYCVK')
effSelection = selectVarsRF.result('LIBMDJQTGU')
save(degreeSelection, bwnSelection, LoTSelection, effSelection, file = '/data/jag/dpustina/APHASIA/STROKE/processing/IITconnectome/VariableSelections.Rdata')

rawSelection = selectVarsRF.result('JPMBTYORZV')


#
#################################################





########################## graph theory things
# densities = seq(0.05,0.90, by=0.05)
# denscor = rep(NA, length(densities)); names(denscor) = paste0('V',densities)
# for (i in 1:length(densities)) { ####################
#                                  dens = densities[i]
                          
graphmat = list()
for (i in 1:nrow(SD)) {
  temp = makeGraph(abs(submats[[i]]), graphdensity = 0.25)
  graphmat[[i]] = c(
                        temp$degree, 
                        temp$betweeness, 
                        temp$localtransitivity,
                        temp$globalTransitivity
                    )
}
graphmatpred = t(sapply(graphmat, rbind ))
# colnames(graphmatpred) = paste0('V',1:268)

cat(paste(format(Sys.time(), "%H:%M") , 'Start RF'))
ntrees=400
rfdata = as.matrix(graphmatpred[,sel])
predgraph = rep(NA, nrow(SD))
require(randomForest)
for (i in 1:nrow(SD)) {
  cat(paste(i," "))
  rfm = randomForest(x = as.matrix(rfdata[ -i, ]),
                     y = SD$WABAQ[-i], 
                     ntree = ntrees)
  predgraph[i] = predict(rfm, newdata = rfdata[ i, ])
}

cor(SD$WABAQ, predgraph)
# } ################################
plot(SD$WABAQ, predgraph)


#####
## variable selection with pAC
graphmat = list()
for (i in 1:nrow(SD)) {
  temp = makeGraph(abs(submats[[i]]), graphdensity = 0.25)
  graphmat[[i]] = c(
#     temp$degree
#     temp$betweeness) 
     temp$localtransitivity)
#     temp$globalTransitivity)
}
graphmatpred = t(sapply(graphmat, rbind ))
colnames(graphmatpred) = paste0('V',1:268)

behav=SD$WABAQ
source('/data/jag/dpustina/Code/APHASIA/AC_pRF.R')
lapply(c('ggplot2', 'permute', 'randomForest', 'reshape2', 'magrittr','multtest','plyr'),require, character.only=T)
runs = 100
ntree=1000 #400 for graphs, 1000 for raw corrs
predictors=data.frame(graphmatpred)
predictors=data.frame(groupmat)
results = matrix(NA, nrow = runs, ncol = ncol(predictors))
colnames(results) = paste0('V',1:ncol(predictors))
significants = rep(0, ncol(predictors))
maxsig = rep(0,runs)
mtry = floor(ncol(predictors)/3)

for (kot in 1:runs) {
  seed = sample(100:10000,1)
  tic=proc.time()
  suppressMessages(
    p.test<-pRF(response=behav,mtry=mtry,
                predictors=predictors,n.perms=300,ntree=ntree,seed=seed,
                type="regression",alpha=0.05)
  )
  proc.time() - tic
  results[kot,order(p.test$Res.table$p.value)] = 1:ncol(predictors)
  sig = p.test$Res.table$p.value < 0.05
  maxsig[kot] = sum(sig)
  significants[sig] = significants[sig]+1
  cat(paste(kot,''))
}

thisorder = order(significants, decreasing= T)
print(paste('Average significant variables:', mean(maxsig)))
MYpredicts = data.frame(AverageOrder = colMeans(results)[thisorder] , SignificanceRatio=(significants/runs)[thisorder])
MYpredicts

degreepredicts = MYpredicts
degreemaxsig = maxsig
degreesignificants = significants
degreeresults=results
degreesel = c(49,239,128,  173, 148, 260, 264, 196, 160, 19, 12)
save(degreepredicts,degreemaxsig,degreesel, degreesignificants, degreeresults,file = 'degreesel.Rdata')

bwnesspredicts = MYpredicts
bwnessmaxsig = maxsig
bwnesssignificants = significants
bwnessresults=results
bwnesssel = c(41,262,  253, 233, 252, 217, 42, 238, 168, 101, 12)
save(bwnesspredicts,bwnessmaxsig,bwnesssel, bwnesssignificants, bwnessresults,file = 'bwnesssel.Rdata')

lotranspredicts = MYpredicts
lotransmaxsig = maxsig
lotranssignificants = significants
lotransresults=results
lotranssel = c(222,173, 42, 220,  255, 24, 84, 94, 230, 38, 192)
save(lotranspredicts,lotransmaxsig,lotranssel, lotranssignificants, lotransresults,file = 'lotranssel.Rdata')

rownames(MYpredicts) = colnames(groupdmg)[thisorder]
MYpredicts
dmgpredicts = MYpredicts
dmgmaxsig = maxsig
dmgsignificants = significants
dmgresults=results
dmgsel = c(171, 181, 158, 165,   184)
save(dmgpredicts,dmgmaxsig,dmgsel, dmgsignificants, dmgresults,file = 'dmgsel.Rdata')

##### 
repeats=30#500
source('/data/jag/dpustina/Code/APHASIA/dorFuncs.R',echo=F)
ctrl = rfeControl(functions = dorFuncs, method = 'repeatedcv', repeats=repeats, verbose=F)
#subsets=c(1:5,seq(10,268, by=10))
subsets = seq(1, 30001, by=5000)
tic=proc.time()
rfProfile = rfe(x=groupmat,
                y=SD$WABAQ,
                rfeControl=ctrl,
                sizes = subsets,
                metric='Cor')

proc.time()-tic
plot(rfProfile, type=c('g','o'))
rfProfile



betweenesssel = as.numeric(gsub('V', '', predictors(rfProfile)))[1:40]
# degreesel=as.numeric(gsub('V', '', predictors(rfProfile)))[1:11] # 49 239 128 173 260 148 264 160 196 149 240





# predict with all selected vars
graphmat = list()
for (i in 1:nrow(SD)) {
  temp = makeGraph(abs(submats[[i]]), graphdensity = 0.25)
  graphmat[[i]] = c(
    temp$degree[degreesel], 
    temp$betweeness[bwnesssel], 
    temp$localtransitivity[lotranssel],
    temp$globalTransitivity
  )
}
graphmatpred = t(sapply(graphmat, rbind ))
colnames(graphmatpred) = c( paste0('Deg',degreesel),
                            paste0('Bwn',bwnesssel),
                            paste0('LoT',lotranssel),
                            'GloT')

cat(paste(format(Sys.time(), "%H:%M") , 'Start RF'))
ntrees=1000

rfdata = as.matrix(graphmatpred)
# rfdata = as.matrix(  rfdata[, grep('GloT', colnames(rfdata))]  ) # activate this for selection
predgraph = rep(NA, nrow(SD))
require(randomForest)
for (i in 1:nrow(SD)) {
  cat(paste(i," "))
  rfm = randomForest(x = as.matrix(rfdata[ -i, ]),
                     y = SD$WABAQ[-i], 
                     ntree = ntrees)
  predgraph[i] = predict(rfm, newdata = rfdata[ i, ])
}

cor(SD$WABAQ, predgraph,method='spearman')
# [1] 0.6433938

predALL = predgraph
predGlot = predgraph

pred1 = predgraph # degree
pred2 = predgraph # bwness
pred3 = predgraph # lotrans
lmdata = data.frame(pred1, # degree
                    pred2, # bwness
                    pred3, # lotrans
                    pred4 = graphmatpred[,34], # glotrans
                    behav = SD$WABAQ
)

lmpred = rep(NA, nrow(SD))
for (i in 1:nrow(SD)) {
  lmod = randomForest(behav ~ ., data=lmdata[-i,])
  lmpred[i] = predict(lmod, newdata=lmdata[i,])
}
cor(lmdata$behav, lmpred,method='spearman')
sqrt(mean((lmdata$behav - lmpred)^2)) # RMSE
mean(abs((lmdata$behav - lmpred)))

plot(SD$WABAQ,lmpred)
###############################

behav = 'WABAQ' ; # WABAQ WABcomp WABrep PNTcorrect
# predicted vector
predvec = predvecdmg = rep(NA,nrow(SD))
for (i in 1:nrow(SD)) {
  tempbehav = eval(parse(text=paste0('SD$',behav,'[-i]')))

  cors = cor(groupmat[ i, ] , t( groupmat[ -i, ] ), method='spearman' )
  predvec[i] = weighted.mean(tempbehav, cors)
  
  cors = cor(groupdmg[ i, ] , t( groupdmg[ -i, ] ), method='spearman' )
  predvecdmg[i] = weighted.mean(tempbehav, cors)
}

sel = !is.na(predvec); print(sum(sel))
bdata = eval(parse(text=paste('SD$',behav,'[sel]')))
corvec = cor(bdata, predvec[sel]); corvec
sel = !is.na(predvecdmg); print(sum(sel))
bdata = eval(parse(text=paste('SD$',behav,'[sel]')))
cordmg = cor(bdata, predvecdmg[sel]); cordmg
par(mfrow=c(1,2))
plot(SD$WABAQ[sel], predvec[sel]); title(paste('Functional Vector\n',behav,'Ploras Style\nr=',round(corvec,2)))
plot(SD$WABAQ[sel], predvecdmg[sel]); title(paste('Regional Lesion\n',behav,'Ploras Style\nr=',round(cordmg,2),' regions=',ncol(groupdmg)))
par(mfrow=c(1,1))


# RF using full upper triangle
predvecRF =  predvecdmgRF = rep(NA,nrow(SD))
ntrees = 1000
# rfdata = cbind(fisherz(groupmat), groupdmg)
rfdata = cbind(fisherz(groupmat))
cat(paste(format(Sys.time(), "%H:%M") , 'Start RF'))
require(randomForest); require(psych)
for (i in 1:nrow(SD)) {
  cat(paste(i," "))
  rfm = randomForest(x = rfdata[ -i, ],
                     y = SD$WABAQ[-i], 
                     ntree = ntrees)
  predvecRF[i] = predict(rfm, newdata = rfdata[ i, ])
}
cat(paste(format(Sys.time(), "%H:%M") , 'Stop RF'))
cor(SD$WABAQ, predvecRF)

predvecRF_500treeAALrestonly = predvecRF
predvecRF_500treeAALdmganrest = predvecRF




# RF on damage only
predvecdmgRF = rep(NA,nrow(SD))
ntrees = 100
cat(paste(format(Sys.time(), "%H:%M") , 'Start RF'))
rfdata = cbind(groupdmg, predvec)
rfdata = groupdmg
# rfdata = cbind(groupdmg, SD$MPOPNT, SD$AGEWAB)
require(randomForest)
for (i in 1:nrow(SD)) {
  cat(paste(i," "))
  rfm = randomForest(x = rfdata[ -i, ],
                     y = SD$WABAQ[-i], 
                     ntree = ntrees)
  predvecdmgRF[i] = predict(rfm, newdata = rfdata[ i, ])
}

cor(SD$WABAQ, predvecdmgRF)

behav='WABAQ'
bdata = eval(parse(text=paste('SD$',behav,'[sel]')))
corvec = cor(bdata, predvecRF_500treeAALrestonly); corvec
cordmg = cor(bdata, predvecRF_500treeAALdmganrest); cordmg
cordmgonly = cor(SD$WABAQ, predvecdmgRF); cordmgonly
par(mfrow=c(1,3))
plot(SD$WABAQ, predvecRF_500treeAALrestonly); title(paste('Functional Vector\n',behav,'RF Style\nr=',round(corvec,2),' variables=',ncol(groupmat))); abline(a=0,b=1)
plot(SD$WABAQ,predvecRF_500treeAALdmganrest); title(paste('Functional + damage Lesion\n',behav,'RF Style\nr=',round(cordmg,2),' variables=',ncol(groupmat))); abline(a=0,b=1)
plot(SD$WABAQ,predvecdmgRF); title(paste('Regional Damage only\n',behav,'RF Style\nr=',round(cordmgonly,2),' variables=',ncol(groupdmg))); abline(a=0,b=1)
par(mfrow=c(1,1))




dmg_rest = predvecdmgRF
dmg_only = predvecdmgRF
behav = SD$WABAQ
method='spearman'
paired.r(xy = cor(behav, dmg_only, method=method),
         xz = cor(behav, dmg_rest, method=method),
         yz = cor(dmg_only, dmg_rest, method=method),
         n=53, twotailed=F)


for (i in 1:10) {
  print(i)
  i = i + 2
}


##### combine all three
# predvec + predvecdmgRF + predgraph
rfdata = data.frame(
                    functavg=predvec,
                    dmgRF=predvecdmgRF,
                    #lmpred,
                    lmdata[,1:3],
                    #graphRF=predgraph,
                    #age=SD$AGEWAB,
                    #mpo=SD$MPOWAB,
                    resp = SD$WABAQ)
multipred = multipredRF = rep(NA, nrow(SD))
for (i in 1:nrow(SD)) {
  lmod = lm(resp ~ ., data=rfdata[-i,])
  multipred[i] = predict(lmod, newdata=rfdata[i,])
  rfmod = randomForest(resp ~ ., data=rfdata[-i,], ntree=1000)
  multipredRF[i] = predict(rfmod, newdata=rfdata[i,])
}
method='spearman'
cor(SD$WABAQ, multipred,method=method)
cor(SD$WABAQ, multipredRF,method=method)
plot(SD$WABAQ, multipred,ylim=c(0,100),xlim=c(0,100)); abline(a=0,b=1)
mean(abs(SD$WABAQ - multipred))
mean(abs(SD$WABAQ - multipredRF))
sqrt(mean((SD$WABAQ-multipred)^2))
sqrt(mean((SD$WABAQ-multipredRF)^2))



behav = SD$WABAQ
method='spearman'
z=multipred; y=multipredall # predvecdmgRF # pred3 pred2 #pred1 # predvec
paired.r(xy = cor(behav, y, method=method),
         xz = cor(behav, z, method=method),
         yz = cor(y, z, method=method),
         n=53, twotailed=F)


# abstract for SNL 2016
lims=c(min(behav)-1, max(behav)+1)
method='spearman'

thisvar = predvec
thistit = "Functional connectivity similarity"
png(filename=paste0(thistit,'.png'),width=300,height=300)
par(mgp=c(2,1,0))
plot(thisvar ~ behav, pch=19, xlab='True', ylab='Predicted'); title(paste( thistit, '\nrho =', round(cor(behav,thisvar, method=method),2) ))
dev.off()


thisvar = predvecdmgRF
thistit = "Structural damage prediction"
png(filename=paste0(thistit,'.png'),width=300,height=300)
par(mgp=c(2,1,0))
plot(thisvar ~ behav, pch=19, xlab='True', ylab='Predicted'); title(paste( thistit, '\nrho =', round(cor(behav,thisvar, method=method),2) ))
dev.off()

thisvar = pred1
thistit = "Graph Theory - Degree"
png(filename=paste0(thistit,'.png'),width=300,height=300)
par(mgp=c(2,1,0))
plot(thisvar ~ behav, pch=19, xlab='True', ylab='Predicted'); title(paste( thistit, '\nrho =', round(cor(behav,thisvar, method=method),2) ))
dev.off()


thisvar = pred2
thistit = "Graph Theory - Betweeness"
png(filename=paste0(thistit,'.png'),width=300,height=300)
par(mgp=c(2,1,0))
plot(thisvar ~ behav, pch=19, xlab='True', ylab='Predicted'); title(paste( thistit, '\nrho =', round(cor(behav,thisvar, method=method),2) ))
dev.off()


thisvar = pred3
thistit = "Graph Theory - Local Transitivity"
png(filename=paste0(thistit,'.png'),width=300,height=300)
par(mgp=c(2,1,0))
plot(thisvar ~ behav, pch=19, xlab='True', ylab='Predicted'); title(paste( thistit, '\nrho =', round(cor(behav,thisvar, method=method),2) ))
dev.off()


thisvar = multipred
thistit = "Combined Prediction"
png(filename=paste0(thistit,'.png'),width=500,height=400)
par(mgp=c(3,1,0),mar=c(5,5,5,1))
plot(thisvar ~ behav, pch=19, xlab='True', ylab='Predicted',ylim=c(0,103),xlim=c(0,100), cex=1.9, cex.axis=1.9, cex.lab=2.5); abline(a=0,b=1); 
title(paste( thistit, '\nrho =', round(cor(behav,thisvar, method=method),2) ), cex.main=2)
dev.off()




# common model
graphmatpredall = cbind(graphmatpred[,1:33],
                        groupdmg, predvec)
rfdata=data.frame(graphmatpredall, resp=SD$WABAQ)
multipredALL = rep(NA, nrow(SD))
for (i in 1:nrow(SD)) {
  rfmod = randomForest(resp ~ ., data=rfdata[-i,], ntree=1000)
  multipredALL[i] = predict(rfmod, newdata=rfdata[i,])
}
cor(SD$WABAQ, multipredALL, method='spearman')