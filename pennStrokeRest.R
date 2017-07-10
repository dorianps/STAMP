# run our stroke patients

# template stuff
template = antsImageRead('/data/jag/dpustina/templates/pennTemplate/template.nii.gz')
tempbrain = antsImageRead('/data/jag/dpustina/templates/pennTemplate/templateBrain.nii.gz')
tempmask = antsImageRead('/data/jag/dpustina/templates/pennTemplate/templateBrainMask.nii.gz')
templaal = antsImageRead('/data/jag/dpustina/templates/pennTemplate/labels/AAL/AAL_Labels.nii.gz')
templfin = antsImageRead('/data/jag/dpustina/templates/pennTemplate/shen_1mm_268_parcellation_finn2015.nii.gz')


args=commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text=args[i]))
}
# mr = 'MR0874'
# date = 1

basedir = '/data/jag/dpustina/APHASIA/OluData/Niftii/'
#ls /data/jag/dpustina/APHASIA/OluData/Niftii/$sub | sort -r | head -n 1
cmd = paste0('ls ', basedir, mr, ' | sort ', ifelse(date==2,'-r ','') , '| head -n 1')
thisdate = system(cmd, intern=T)
mprdir = Sys.glob(paste0(basedir,mr,'/',thisdate,'/rawNii/'))
bolddir = Sys.glob(paste0(basedir,mr,'/',thisdate,'/RESTBOLD/'))

# load bold
fn = Sys.glob(paste0(bolddir, '*bold*.nii.gz'))
if (length(fn)==0) stop(paste(mr,'No bold images found!'))
img  = antsImageRead( fn[1])
extraRuns = list()
if (length(fn) > 1) {
  for (e in 2:length(fn)) extraRuns[[e-1]] = antsImageRead(fn[e])
} else {
  extraRuns=NA
}

# load segmentation
sfn = Sys.glob(paste0(mprdir, 'ToPennTemplate025BrainSegmentation.nii.gz')) 
if (length(sfn)==0) stop(paste(mr,'No segmentation found!'))
seg = antsImageRead(sfn)

# load t1
t1fn = Sys.glob(paste0(mprdir, 'ToPennTemplate025ExtractedBrain0N4.nii.gz')) 
if (length(t1fn)==0) stop(paste(mr,'No N4 found!'))
t1 = antsImageRead( t1fn   )

# load lesion
lfn = Sys.glob(paste0(mprdir, 'myMaskNative.nii.gz')) 
if (length(lfn)==0) stop(paste(mr,'No lesion found!'))
lesion = antsImageRead( lfn )
if ( any(! antsGetOrigin(lesion) == antsGetOrigin(seg)) ) antsCopyImageInfo(seg,lesion)  # problems with recent OLU masks

# find transforms
subpennmat = Sys.glob(paste0(
  '/data/jag/dpustina/APHASIA/STROKE/processing/LesionMaskedRegistration/',
  mr, '/pennTemplate_to_', mr, '_0GenericAffine.mat') )
subpennwarp = Sys.glob(paste0(
  '/data/jag/dpustina/APHASIA/STROKE/processing/LesionMaskedRegistration/',
  mr, '/pennTemplate_to_', mr, '_1Warp.nii.gz') )
tempreg = c(subpennwarp,subpennmat)
if (length(tempreg)!=2) stop(paste(mr,'Cannot find transforms!'))




####################
# run strokeRest
######################

source('/data/jag/dpustina/Code/APHASIA/strokeRest.R', echo=F)

steadyT=12
fdthresh = 0.3
repeatMotionEst = 2
nCompCor = 4
polydegree = 4
freqLimits = c(0.009,0.08)
smoothingSigmas = c(5.2, 5.2, 5.2, 0.000000)

rsnetwork = strokeRest(img = img,
                        fdthresh = fdthresh, 
                        repeatMotionEst = repeatMotionEst,
                        polydegree = polydegree,
                        nCompCor = nCompCor,
                        smoothingSigmas = smoothingSigmas,
                        freqLimits = freqLimits,
                        globalMeanNuis = F,
                        
                        structuralImage = t1, 
                        structuralSeg = seg, 
                        structuralLesion = lesion,
                        
                        tempreg  = tempreg,
                        template = template,
                        tempbrain = tempbrain,
                        tempbmask = tempmask,
                        tempparc = list(templaal,templfin),
                        
                        extraRuns = extraRuns,
                        
                        verbose = T)


rsnetwork$fn = fn
rsnetwork$sfn = sfn
rsnetwork$t1fn = t1fn
rsnetwork$lfn = lfn
rsnetwork$mr = mr
rsnetwork$date = date
rsnetwork$thisdate = thisdate

# save output
outdir = paste0('/data/jag/dpustina/APHASIA/STROKE/processing/RESTnoGLOB/',mr,'/',thisdate)
dir.create(outdir, showWarnings = F, recursive = T)
setwd(outdir)
write.table(rsnetwork$nodeCor[[1]], file = 'AALnetwork116.txt', col.names = F, row.names = F)
write.table(rsnetwork$nodeCor[[2]], file = 'FINnetwork268.txt', col.names = F, row.names = F)
save.ANTsR(objects='rsnetwork', filename = './strokeRest', overwrite = T, clonediskfiles = T)



struaal = antsApplyTransforms(fixed=t1, moving = templaal, transformlist = tempreg, interpolator = 'MultiLabel')
strufin = antsApplyTransforms(fixed=t1, moving = templfin, transformlist = tempreg, interpolator = 'MultiLabel')



onefile=T
pdf(file = paste0('Summary_',mr,'_date_',date,'.pdf'),
    width=15, height=7)


#########################
# motion params
#########################
par(mfrow=c(1,1))
moco=rsnetwork$moco
reg_params <- as.matrix(moco$moco_params)
nTimes = nrow(moco$moco_params)
tr = antsGetSpacing(rsnetwork$fusedImg)[4]


require(ggplot2)
orderedBreaks = c("Framewise", "X", "Y", "Z", "Pitch", "Roll", "Yaw",'DVARS' )
moco.dat <- data.frame(Time=rep(1:nTimes, 7)*tr)
moco.dat$Values = c(moco$fd$MeanDisplacement, as.vector(reg_params)  )
moco.dat$Category = c( rep("Displacement", 4*nTimes),rep("Angle", 3*nTimes) )
moco.dat$Type = rep(c("Framewise","X", "Y", "Z","Pitch", "Roll", "Yaw"), each=nTimes)

mocodvar = data.frame(Time=(1:nTimes)*tr,
                      Values=rsnetwork$dvars,
                      Category=(rep('DVARS',nTimes)),
                      Type=(rep('DVARS',nTimes)))
moco.dat=rbind(moco.dat,mocodvar)

#sel=moco.dat$Category=='Displacement'
regPlot <- ggplot(moco.dat, aes(x=Time, y=Values, group=Type, colour=Type) )
regPlot <- regPlot + geom_line(size=0.5)
regPlot <- regPlot + theme(text=element_text(size=15), legend.position="bottom")
regPlot <- regPlot + ggtitle(paste(mr, thisdate, 'date', date,
                                   paste('\n', ifelse(is.na(rsnetwork$badtimes), 'NA', length(unique(rsnetwork$badtimes))),'bad volumes,', nrow(rsnetwork$boldMat), 'left'),
                                   "\nMotion correction parameters"))
regPlot <- regPlot + facet_grid(Category ~ ., scales="free" )
regPlot <- regPlot + scale_color_discrete(breaks=orderedBreaks)
regPlot <- regPlot + geom_hline( yintercept=fdthresh, linetype="dashed")

if (!is.na(rsnetwork$badtimes)) {
  badstarts = rsnetwork$badtimes[seq(1,length(rsnetwork$badtimes),2)]
  badstops = rsnetwork$badtimes[seq(2,length(rsnetwork$badtimes),2)]
  bad.data.rect = data.frame(Start=badstarts*tr, Stop=badstops*tr)
  regPlot <- regPlot + annotate("rect", xmin=badstarts*tr, xmax=badstops*tr, ymin=-Inf, ymax=Inf, fill="black", alpha=0.2)
}


print(regPlot)
####


#########################
# connectivity marices
#########################
par(mar=c(0,0,0,0))

col1 <- colorRampPalette(rev(c("#7F0000","red","#FF7F00","yellow","white", 
                               "cyan", "#007FFF", "blue","#00007F")))
col <- colorRampPalette(c('white', "black"))

require(corrplot)
# corrplot(rsnetwork$nodeCor[[1]], method='square',
#          tl.pos='n',cl.pos='b', col=col1(200),
#          main='AAL Connect Matrix')
dmgvec = rsnetwork$nodeDamage[[1]]
dmgmat = rsnetwork$nodeCor[[1]]*0
dmgmat = apply(dmgmat,2, function(x) x+dmgvec)
dmgmat = apply(dmgmat,1, function(x) x+dmgvec)
dmgmat[dmgmat>1]=1
par(mar=c(0,0,0,0))
# corrplot(dmgmat, method='square',tl.pos='n',cl.pos='b', 
#          col=col(200),is.corr=F,
#          main='AAL node damages')

plot1 = pheatmap::pheatmap( rsnetwork$nodeCor[[1]], cluster_rows = F,
                            cluster_cols=F, border_color=NA,
                            main='AAL node matrix',silent=T)
plot2 = pheatmap::pheatmap(dmgmat, cluster_rows = F,
                           cluster_cols=F, border_color=NA,
                           main='AAL node damage',silent=T,
                           color=col(200))
require(gridExtra)
grid.arrange(plot1$gtable,plot2$gtable,ncol=2)



# corrplot(rsnetwork$nodeCor[[2]], method='square',
#          tl.pos='n',cl.pos='b', col=col1(200),
#          main='Fin Connect Matrix')
dmgvec = rsnetwork$nodeDamage[[2]]
dmgmat = rsnetwork$nodeCor[[2]]*0
dmgmat = apply(dmgmat,2, function(x) x+dmgvec)
dmgmat = apply(dmgmat,1, function(x) x+dmgvec)
dmgmat[dmgmat>1]=1
par(mar=c(0,0,0,0))
# corrplot(dmgmat, method='square',
#          tl.pos='n',cl.pos='b', col=col(200),is.corr=F,
#          main='Fin node damages')


plot1 = pheatmap::pheatmap( rsnetwork$nodeCor[[2]], cluster_rows = F,
                    cluster_cols=F,border_color=NA,
                    main='Fin node matrix',silent=T)
plot2 = pheatmap::pheatmap(dmgmat, cluster_rows = F,
                            cluster_cols=F,border_color=NA,
                            main='Fin node damage',silent=T,
                           color=col(200))
grid.arrange(plot1$gtable,plot2$gtable,ncol=2)






# print images
par(mar=c(1,0,2,0), mfrow=c(1,2))

# t1 + lesion
plot(t1,window.img=c(0.1,max(t1)))
title(paste(mr,'date',date,'- T1'))
plot(t1,lesion,window.img=c(0.1,max(t1)),window.overlay=0:1,alpha=0.7)
title(paste(mr,'date',date,'- Lesion'))

# t1 + parcels
plot(t1, struaal, alpha=0.7, window.overlay=c(1,116))
title(paste(mr,'date',date,'- AAL on T1'))
plot(t1, strufin, alpha=0.7, window.overlay=c(1,268))
title(paste(mr,'date',date,'- Fin on T1'))

# bold +/- T1
plot(rsnetwork$meanbold, rsnetwork$boldmap$warpedmovout, window.img=c(0.1,max(rsnetwork$meanbold)), alpha=0.7, window.overlay=c(0.1,max(rsnetwork$boldmap$warpedmovout)))
title(paste(mr,'date',date,'- T1 on Bold'))
plot(t1, rsnetwork$boldmap$warpedfixout, window.img=c(0.1,max(t1)), alpha=0.7, window.overlay=c(1,max(rsnetwork$boldmap$warpedfixout)))
title(paste(mr,'date',date,'- Bold on T1'))

# bold + lesion
plot(rsnetwork$meanbold,
     window.img=c(0.1,max(rsnetwork$meanbold)))
title(paste(mr,'date',date,'- Mean Bold'))
plot(rsnetwork$meanbold,rsnetwork$boldlesion, 
     window.img=c(0.1,max(rsnetwork$meanbold)),
     window.overlay=0:1,alpha=0.7)
title(paste(mr,'date',date,'- Bold Lesion'))


# bold + segmentation
plot(rsnetwork$meanbold,
     window.img=c(0.1,max(rsnetwork$meanbold)))
title(paste(mr,'date',date,'- Mean Bold'))
plot(rsnetwork$meanbold,rsnetwork$seg2bold, 
     window.img=c(0.1,max(rsnetwork$meanbold)),
     window.overlay=c(0,max(rsnetwork$seg2bold)),alpha=1)
title(paste(mr,'date',date,'- Seg2Bold \n(used for nuisance and global mean)'))


# bold + parcels
plot(rsnetwork$meanbold, rsnetwork$boldlabels[[1]], 
     window.img=c(0.1,max(rsnetwork$meanbold)),alpha=0.7, 
     window.overlay=c(1,116))
title(paste(mr,'date',date,'- AAL on Bold'))
plot(rsnetwork$meanbold, rsnetwork$boldlabels[[2]], 
     window.img=c(0.1,max(rsnetwork$meanbold)), alpha=0.7, 
     window.overlay=c(1,268))
title(paste(mr,'date',date,'- Fin on Bold'))






dev.off()