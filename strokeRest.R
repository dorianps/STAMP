
strokeRest <- function(
  img,  
  extraRuns = NA,
  
  
  polydegree = 4,
  steadyT=12.0,
  repeatMotionEst = 1,
  mocoTxType = "BOLDRigid",
  fdthresh=Inf,
  nCompCor = 4,
  freqLimits = c( 0.01, 0.1 ),
  globalMeanNuis = T,
  
  
  smoothingSigmas = NA,
  
  structuralImage = NA,
  structuralSeg = NA,
  structuralLesion = NA,
  structuralMask=NA,
  
  tempreg  = NA,
  template = NA,
  tempbrain = NA,
  tempbmask = NA,
  tempparc = NA,
  
  verbose = FALSE)
{
  
  if (class(img) != "antsImage" | length(dim(img)) != 4)
    stop("bold must be antsImage with four dimensions")
  if ( !is.na(extraRuns) & class( extraRuns )[[1]] != "list"  )
    stop("extraRuns must be a list of antsImages.")
  if (is.na(structuralImage) | class(structuralImage) != "antsImage")
    stop("structural image required ")
  
  
  output = list()
  output$start = date()
  tic = proc.time()
  
  
  # eliminate non steady state timepoints
  if (steadyT > 0) {
    emp = steadyRemove(img,extraRuns,steadyT)
    img = emp$img
    extraRuns = emp$extraRuns
    if (verbose) cat(paste(format(Sys.time(), "%H:%M"), "Non steady volumes removed...\n"))
  }
  
  
  
  # create runNuis and timevals
  tr = antsGetSpacing(img)[4]
  runNuis = rep(1, dim(img)[4] )
  allTimes = dim(img)[4]
  if ( ! all( is.na( extraRuns ) ) ) {
    for ( i in 1:length( extraRuns ) ) {
      allTimes = allTimes + dim( extraRuns[[i]] )[4]
      runNuis = c( runNuis, rep(i+1, dim(extraRuns[[i]])[4] ) )
    }
  }
  runNuis = factor( runNuis )
  timevals = 0
  for ( runlev in levels( runNuis ) ) {
    if ( length( timevals ) == 1 ) {
      timevals = as.numeric( 1:sum( runNuis == runlev ) ) * tr
    } else {
      timevals = c( timevals, as.numeric( 1:sum( runNuis == runlev ) ) * tr )
    }
  }
  
  
  output$timevals = timevals
  
  
  
  
  
  # motion correct and concatenate
  if (verbose) cat(paste(format(Sys.time(), "%H:%M"), "Running motion correction...\n"))
  emp = mocoFuse(img,extraRuns,mocoTxType, repeatMotionEst)
  moco = emp$moco
  meanbold = emp$meanbold
  
  output$meanbold = meanbold
  output$moco = moco
  if (repeatMotionEst > 0) output$FD = moco$fd$MeanDisplacement
  
  
  
  
  # need to segment wm,gm,csf
  if (is.na(structuralSeg)) {
    if (verbose) cat(paste(format(Sys.time(), "%H:%M"), "Segmenting structural image...\n"))
    n4 = abpN4(t1)
    abp = abpBrainExtraction(img = n4, tem = template, temmask = tempmask)
    seg = kmeansSegmentation(abp$brain, 3)
    structuralSeg = seg$segmentation
    rm(abp,seg)
    output$n4 = n4
    output$structuralSeg = seg
  }    
  
  structmask = thresholdImage( structuralSeg, 1, Inf )
  t1brain = structuralImage * structmask
  
  output$t1brain = t1brain
  output$structmask = structmask
  
  
  # coregister bold to structural, bring back segmentation
  if (verbose) cat(paste(format(Sys.time(), "%H:%M"), "Coregistering bold to structural...\n"))
  boldmap = antsRegistration( meanbold, t1brain,
                              typeofTransform='SyNBoldAff', verbose=FALSE )
  tempstructseg = structuralSeg*1
  if (!is.na(structuralLesion)) tempstructseg = structuralSeg * abs(structuralLesion-1) #remove lesion from WM/CSF segmentation
  seg2bold = antsApplyTransforms( meanbold, tempstructseg, boldmap$fwdtransforms,
                                  interpolator = "MultiLabel" )
  if (is.na(structuralMask)) structuralMask = structmask
  mask2bold = antsApplyTransforms(meanbold, structuralMask, boldmap$fwdtransforms,
                                  interpolator = "NearestNeighbor" )
  
  output$seg2bold = seg2bold # registration maps saved later
  output$mask2bold = mask2bold 
  output$boldmap = boldmap
  
  
  
  
  
  # fuse runs with poly
  if (verbose) cat(paste(format(Sys.time(), "%H:%M"), "Fusing with poly...\n"))
  emp = polyRegres(moco, mask2bold, runNuis, timevals, polydegree)
  fusedImg = emp$fusedImg
  polyNuis = emp$polyNuis
  
  output$polyNuis = polyNuis

  
    
    
  # compute badtimes
  if (repeatMotionEst > 0) {
    if (verbose) cat(paste(format(Sys.time(), "%H:%M"), "Computing bad times...\n"))
    emp = badTimes(fusedImg=fusedImg, mask=mask2bold, moco=moco, fdthresh=fdthresh)
    boldMat = emp$boldMat # with badtimes interpolated
    badtimes = emp$badtimes
    goodtimes = emp$goodtimes
    dvars = emp$dvars
    # done in function if (!any(is.na(badtimes))) badtimes = sort(c(badtimes, badtimes+1))
    if (length(badtimes) > (nrow(boldMat)/10)) warning("Careful! More than 10% bad timepoints.")

  } else {
    goodtimes = 1:dim(fusedImg)[4]
    badtimes = NA
    boldMat = timeseries2matrix(fusedImg, mask2bold)
    dvars <- computeDVARS( timeseries2matrix( fusedImg, mask2bold ) )
  }
    
  output$badtimes = badtimes
  output$goodtimes = goodtimes
  output$dvars = dvars
  output$runID = runNuis
  
    
  # white matter is labeled as 3
  if (verbose) cat(paste(format(Sys.time(), "%H:%M"), "Computing tissue nusiances...\n"))
  wmMask = thresholdImage(seg2bold, 3,3)
  wmMask = iMath( wmMask, "ME", 1)
  wmVox = which(subset(wmMask, mask2bold > 0 )==1)
  wmMean = rowMeans(boldMat[,wmVox])
  
  csfMask = thresholdImage(seg2bold, 1,1)
  csfVox = which(subset(csfMask, mask2bold > 0)==1)
  csfMean= rowMeans(boldMat[,csfVox])
  

  
  tissueNuis = cbind( wmMean, csfMean )
    
  if (globalMeanNuis) {
    nonlesvox = imageListToMatrix(
      imageList = list(thresholdImage(seg2bold,1,Inf)), 
        mask=mask2bold)==1
    globalMean = rowMeans(boldMat[,nonlesvox])
    tissueNuis = cbind( globalMean, tissueNuis )
  }
    
    
  # interpolate badtimes for nuisances
  nTimes=nrow(boldMat)
  if ( any(!is.na(badtimes)) ) {
    for ( v in c(1:dim(tissueNuis)[2]) ) {
      tissueInterp = spline( c(1:nTimes)[goodtimes], tissueNuis[goodtimes,v],
                             method='natural', xout=badtimes )$y
      tissueNuis[badtimes,v]=tissueInterp
    }
  }
  tissueDeriv = rbind( rep(0,dim(tissueNuis)[2]), diff(tissueNuis,1) )
  
  nuisance = cbind( tissueNuis, 
                    tissueDeriv, 
                    dvars=dvars)
  
    
    
  # add moco and compcors to nuissance
  if (repeatMotionEst > 0) {
    reg_params <- as.matrix( moco$moco_params )
    mocoNuis              = reg_params
    mocoNuis2             = reg_params * reg_params
    colnames( mocoNuis2 ) = paste( "MocoSqr", 1:ncol( mocoNuis2 ) , sep='' )
    mocoNuis  = cbind( mocoNuis, mocoNuis2 )
    mocoDeriv = rbind( rep( 0,  dim(mocoNuis)[2] ), diff( mocoNuis, 1 ) )
    colnames( mocoDeriv ) = paste( "MocoDeriv", 1:ncol( mocoDeriv ) , sep='' )
    nuisance = cbind( mocoNuis,
                      mocoDeriv, 
                      nuisance)
  }
  

  if ( ! all( is.na( extraRuns ) ) )  nuisance = cbind(nuisance, runs=runNuis)
  
  if ( nCompCor > 0 ) {
    if (verbose) cat(paste(format(Sys.time(), "%H:%M"), "Computing compcor nuisances...\n"))
    # use seg2bold and moco_img to get a better compcor set "anatomical compcor"
    # http://www.ncbi.nlm.nih.gov/pubmed/25987368
    tempMask = thresholdImage( seg2bold, 1, 1 )
    tempMask = (tempMask + thresholdImage( seg2bold, 3, 3 ) ) %>% iMath( "ME", 1 )
    tempMat = timeseries2matrix( fusedImg, tempMask )
    compcorNuis = compcor( tempMat, nCompCor )
    colnames( compcorNuis ) = paste("compcor",1:ncol(compcorNuis), sep='' )
    nuisance = cbind( nuisance, compcorNuis )
    rm( tempMat  )
    rm( tempMask )
  }
    
  output$nuisance = nuisance  
    
  # regress out nuissances from already fixed boldMat
  boldMat = t(boldMat)
  boldMat[goodtimes,] <- residuals( lm( boldMat[goodtimes,] ~ nuisance[goodtimes,] ) )
  boldMat = t(boldMat)
  
  
  # fill bad times with interpolated points
  if ( any(!is.na(badtimes)) ) {
    nVox = length(which(as.array(mask2bold)==1))
    for ( v in c(1:nVox) ) {
      boldMat[badtimes,v]=spline( c(1:nTimes)[goodtimes], boldMat[goodtimes,v],
                                  method='natural', xout=badtimes )$y
    }
  }
  fusedImg = matrix2timeseries( fusedImg, mask2bold, boldMat )
    

  # smooth images
  if (verbose) cat(paste(format(Sys.time(), "%H:%M"), "Smoothing images...\n"))
  if ( any( is.na( smoothingSigmas ) ) ) {
    sptl    = sqrt( sum( antsGetSpacing(img)[1:3]^2  )) * 1.5
    smoothingSigmas = c( sptl, sptl, sptl, 1.0 )
  }
  fusedImg = smoothImage( fusedImg, smoothingSigmas, FWHM=TRUE )
  
  
  
  # frequency filtering
  if (verbose) cat(paste(format(Sys.time(), "%H:%M"), "Frequency filtering...\n"))
  boldMat = timeseries2matrix( fusedImg, mask2bold )
  if (  ( length( freqLimits ) == 2  ) & ( freqLimits[1] < freqLimits[2] ) ) {
    locruns = levels( runNuis )
    boldMatTemp <- frequencyFilterfMRI( boldMat[ runNuis == locruns[1], ],
                                        tr=tr, freqLo=freqLimits[1], freqHi=freqLimits[2], opt="trig" )
    if ( length( locruns ) > 1 ) {
      for ( myrun in locruns[2:length(locruns)] ) {
        boldMatTemp2 <- frequencyFilterfMRI( boldMat[ runNuis == myrun, ],
                                             tr=tr, freqLo=freqLimits[1], freqHi=freqLimits[2], opt="trig" )
        boldMatTemp = rbind( boldMatTemp , boldMatTemp2 )
      }
    }
    boldMat = boldMatTemp
  }
  
  
  
  # one last global demeaning if necessary
  # only good volumes left after this point
  if (globalMeanNuis) boldMat = residuals( lm( boldMat[goodtimes, ] ~ globalMean[goodtimes] ) )
  
  
  # put bolMat in fusedImg
  fusedImg = matrix2timeseries( fusedImg, mask2bold, boldMat )
  
  
  output$fusedImg = fusedImg
  
  
  # STOP HERE IF NO PARCELLATION FED
  if (any(is.na(tempparc))) return(output)
  
  
  
  ###############################
  # Bring nodes in BOLD space and
  # produce correlation matrix
  ###############################
  
  if (verbose) cat(paste(format(Sys.time(), "%H:%M"), "Node correlations...\n"))
  
  if (is.na(tempbrain) &
        !is.na(template) &
        !is.na(tempbmask)) tempbrain = template*tempbmask
  if (is.na(tempbrain)) stop('Cannot run node correlations without specifying tempbrain')
    
    
  # register template to subject with lesion mask
  if (any(is.na(tempreg))) {
    if (verbose) cat(paste(format(Sys.time(), "%H:%M"), "Template registration...\n"))
    maskreg = thresholdImage(structuralSeg, 1, Inf)
    if (!is.na(structuralLesion)) {
      maskreg[structuralLesion==1] = 0
    }
    tempregister = antsRegistration(fixed=t1brain,moving=tempbrain,typeofTransform = 'SyN',mask=maskreg)
    tempreg  = tempregister$fwdtransforms
    
    output$temp2structmap = tempreg
    output$bold2tempmap = c(tempregister$invtransforms, boldmap$invtransforms)
    output$temp2boldmap = c(boldmap$fwdtransforms , tempreg)
    output$bold2structmap = boldmap$invtransforms
  } else {
    output$temp2structmap = tempreg
    output$temp2boldmap = c(boldmap$fwdtransforms , tempreg)
    output$bold2structmap = boldmap$invtransforms
  }
  
  

  # this is already trimmed and filtered to have only
  # goodtimes but all bold voxels
  output$boldMat = boldMat
  
  
  
  # bring labels in bold space and compute correlations
  if (class(tempparc) == 'antsImage') tempparc = list(tempparc)
  temp2boldreg = c(  boldmap$fwdtransforms , tempreg)
  output$boldlabels = list()
  output$roisignals = list()
  output$nodeCor = list()
  output$nodeDamage = list()
  
  for (nparc in 1:length(tempparc)) {  # loop through parcellation schemes
    # bring this parcellation in bold space
    boldlabel = antsApplyTransforms(fixed = meanbold,
                                    moving = tempparc[[nparc]],
                                    transformlist = temp2boldreg,
                                    interpolator='MultiLabel')
    boldlabel = boldlabel * mask2bold
    output$boldlabels[[nparc]] = boldlabel
    
    # compute roi signals
    roimat = labels2matrix(boldlabel,mask2bold)
    if (any(roimat<0)) stop("cannot accept parcellation with negative value")
    if (any(unique(c(as.numeric(boldlabel))) == 0))  roimat = roimat[-1,] # remove label 0

    roisignals = matrix(NA, nrow=nrow(boldMat), ncol=nrow(roimat))
    for (i in 1:ncol(roisignals)) {
      roisignals[,i] = rowMeans(boldMat[, roimat[i,]==1])
    }
    output$roisignals[[nparc]] = roisignals
    
    
    # CORRELATION MATRIX
    nodeCor = cor(roisignals)
    output$nodeCor[[nparc]] = nodeCor
    
    
    
    
    # compute node damages
    if (!is.na(structuralLesion)) {
      if (!exists('boldlesion')) 
          boldlesion = antsApplyTransforms(fixed = meanbold,
                                      moving = structuralLesion,
                                      transformlist = boldmap$fwdtransforms,
                                      interpolator='NearestNeighbor')
      lesmat = imageListToMatrix( list(boldlesion), mask2bold )
      nodeLesion = roimat %*% t(lesmat)
      nodeLesion = nodeLesion / rowSums(roimat)
      
      output$nodeDamage[[nparc]] = nodeLesion
      output$boldlesion = boldlesion
    }
    

  } 
  
  if (verbose) cat(paste(format(Sys.time(), "%H:%M"), "Done!\n"))
  
  
  # add runtime info
  output$end = date()
  output$runtime = proc.time() - tic
  
  return(output)
    
}
  
  
  
  
  
  
  
  
# remove non steady img functions  
steadyRemove <- function(img, extraRuns, steadyT) {
  
  tr = antsGetSpacing(img)[4]
  steady = floor( steadyT / tr) + 1
  
  img = cropIndices(img, c(1,1,1,steady), dim(img) )
  if ( ! all( is.na( extraRuns ) ) )
  {
    if ( class( extraRuns )[[1]] != "list"  )
      stop("extraRuns must be a list of antsImages.")
    for ( i in 1:length( extraRuns ) )
    {
      timg = extraRuns[[i]]
      timg = cropIndices( timg, c(1,1,1,steady), dim(timg) )
      extraRuns[[i]] = timg
    }
  }
  return(list(img=img, extraRuns=extraRuns))
}
  
  
  
# motion correction with fusion  
mocoFuse <- function(img, extraRuns, mocoTxType="BOLDRigid", repeatMotionEst=1) {
  
  meanbold = getAverageOfTimeSeries( img )
  
  if (repeatMotionEst > 0) {
    for ( i in 1:repeatMotionEst ) {
      moco <- antsrMotionCalculation( img, fixed=meanbold, typeofTransform=mocoTxType )
      meanbold = apply.antsImage( moco$moco_img, c(1,2,3), mean)
    }
    
    # now add any additional runs and merge moco results
    if ( ! all( is.na( extraRuns ) ) ) {
      for ( i in 1:length( extraRuns ) ) {
        timg = extraRuns[[i]]
        meanbold2 = apply.antsImage(timg, c(1,2,3), mean)
        for ( i in 1:repeatMotionEst ) {  # within run moco
          mocoTemp <- antsrMotionCalculation( timg, fixed=meanbold2, typeofTransform=mocoTxType )
          meanbold2 = apply.antsImage(mocoTemp$moco_img, c(1,2,3), mean)
        }
        
        # register to 1st run? keep this idea for later
        # reg2first = antsRegistration(fixed=meanbold, moving=meanbold2, typeofTransform=mocoTxType)
        # motoTemp$moco_img = antsApplyTransforms()
        
        if ( usePkg("abind") ) {
          ttmo = as.array( moco$moco_img )
          ttmo = abind::abind( ttmo, as.array( mocoTemp$moco_img ) )
          moco$moco_img = antsCopyImageInfo( moco$moco_img, as.antsImage(ttmo) )
          rm( ttmo )
        } else stop( "need abind package for the extraRuns feature")
        moco$moco_params = rbind( moco$moco_params, mocoTemp$moco_params )
        moco$fd = rbind( moco$fd, mocoTemp$fd )
        moco$dvars = c( moco$dvars, mocoTemp$dvars )
        rm( mocoTemp )
      }
    }
    
  } else { # no motion correction performed
    moco = list()
    moco$moco_img = antsImageClone(img)
    for ( i in 1:length( extraRuns ) ) {
      ttmo = as.array( extraRuns[[i]] )
      ttmo = abind::abind( ttmo, as.array( moco$moco_img ) )
      moco$moco_img = antsCopyImageInfo( moco$moco_img, as.antsImage(ttmo) )
    }
  }
  
  return( list( moco=moco, meanbold=meanbold ) )
  
}
  
  
  
  
  

############################################################
# now use polynomial regressors to match images across runs #
#############################################################
polyRegres <- function (moco, mask, runNuis, timevals, polydegree) {
  polyNuis = NA
  if ( !is.na( polydegree ) ) # polynomial regressors
  {
    boldMat = timeseries2matrix( moco$moco_img, mask )
    meanmat = boldMat * 0
    rboldMat = boldMat * 0
    mycolmean = colMeans( boldMat )
    mycolmean = mycolmean * 1000 / mean( mycolmean ) # map to 1000
    for ( i in 1:nrow( meanmat ) ) meanmat[i,]=mycolmean
    for ( runlev in levels( runNuis ) )
    {
      polyNuis <- stats::poly( timevals[ runNuis == runlev ], degree = polydegree )
      rboldMat[ runNuis == runlev, ] <- residuals( lm( boldMat[ runNuis == runlev, ] ~ polyNuis ) )
    }
    # residualize but also keep the scale mean
    boldMat = rboldMat + meanmat
    rm( meanmat )
    rm( rboldMat )
    polyNuis <- stats::poly( timevals, degree = polydegree )
  }
  
  fusedImg = matrix2timeseries( moco$moco_img, mask, boldMat )
  
  return(list( fusedImg=fusedImg, polyNuis=polyNuis ))
}
  
  
  
# bad times function
badTimes <- function(fusedImg, mask, moco, fdthresh) {
  nVox = length(which(as.array(mask)==1))
  dvars <- computeDVARS( timeseries2matrix( fusedImg, mask ) )
  goodtimes = (1:dim(moco$moco_img)[4])
  badtimes = which(moco$fd$MeanDisplacement > fdthresh )
  badtimes = sort(c(badtimes, badtimes+1))
  haveBadTimes = FALSE
  if ( length( badtimes ) > 0 ) {
    goodtimes = goodtimes[-badtimes]
    haveBadTimes = TRUE
  } else {
    badtimes = NA
  }
  
  boldMat = timeseries2matrix( fusedImg, mask )
  nTimes = nrow( boldMat )
  if ( haveBadTimes & fdthresh != Inf ) {
    for ( v in c(1:nVox) )
    {
      boldMat[badtimes,v]=spline( c(1:nTimes)[goodtimes], boldMat[goodtimes,v],
                                  method='natural', xout=badtimes )$y
    }
  }
  
  return(
    list(
      boldMat=boldMat,
      badtimes=badtimes,
      goodtimes=goodtimes,
      dvars=dvars
    )
  )
}
  







save.ANTsR <- function(filename=file.path('.','.ANTsRsession'),
                       objects=NA,
                       env=as.environment(1),
                       overwrite=F, 
                       clonediskfiles=T,
                       ...) {
  # convert to absolute path
  filename = suppressWarnings(file.path(dirname(normalizePath(filename)),basename(filename)))
  
  # create or empty the target folder
  if (file.exists(file.path(filename,'ANTSLOAD.Rdata')) & overwrite ) {
    fnames = list.files(filename)
    assign('fremove1234567890', file.path(filename,fnames) , envir = env)
  } else {
    dir.create(filename,showWarnings = F)
    assign('fremove1234567890', '' , envir = env)
  }
  if (file.exists(file.path(filename,'ANTSLOAD.Rdata')) & ! overwrite ){
    stop(paste('Folder', filename, 'not empty and overwrite is false.'))
  }
  
  
  antslist = as.list(env)
  if(all(!is.na(objects))) {
    index = match(objects,names(antslist))
    antslist = antslist[index]
  }
  
  
  funimgS = function(x,fold=filename) {
    file = paste0(paste(sample(c(0:9, letters, LETTERS), 20, replace=T),collapse=''),'.nii.gz')
    fn = file.path(fold,file)
    antsImageWrite(x, fn)
    return(paste0('ANTSload',file))
  }
  
  funimgSf = function(x,fold=filename) {
    index = which(file.exists(x))
    if (length(index) == 0) return(x)
    
    nocopy = file.exists(file.path(fold,basename(x[index])))
    noremovef = file.path(fold,basename(x[index[nocopy]]))
    noremoveindx = match(noremovef, fremove1234567890 )
    assign('fremove1234567890', fremove1234567890[-(noremoveindx)], envir = env)
    x[index[nocopy]] = paste0('ANTSrepl', basename(x[index[nocopy]]) )
    index = index[! nocopy]
    if (length(index) == 0) return(x)
    
    for (indx in index) {      
      file = paste0(paste(sample(c(0:9, letters, LETTERS), 20, replace=T),collapse=''),
                    '_', basename(x[indx]))
      fn = file.path(fold,file)
      file.copy(x[indx],fn)
      x[indx] = paste0('ANTSrepl', file)
    }
    return(x)
  }
  
  if (clonediskfiles) antslist = rapply(antslist, funimgSf, classes='character', how='replace')
  ANTSLOAD = rapply(antslist, funimgS, classes='antsImage', how='replace')
  
  if (file.exists(file.path(filename,'ANTSLOAD.Rdata')) && overwrite && length(fremove1234567890)>0) {
    drop=file.remove(fremove1234567890) # cleanup remaining files in folder
    rm(fremove1234567890,envir = env)
  }
  
  save(ANTSLOAD,file=file.path(filename,'ANTSLOAD.Rdata'), ...)
  
}



