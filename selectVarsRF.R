# function for selecting variables with RFs and permutations
# expanded from the original pRF package in CRAN to have
# multiple runs for a more reliable estimate.
# Produces similar results to caret's recursive feature elimination.
# This implementation designed specifically for use with my qsuv
selectVarsRF <- function(x, y, ntree, nperms, nruns, 
                         memory='1G', 
                         savedir='/data/jag/dpustina/TEMP', 
                         alpha=0.05, noprint=F ) {
  
  jobid = paste(sample(LETTERS, 10, replace=T),collapse='')
  mtry = floor(ncol(x)/3)
  DATAfile = file.path(savedir,paste0('selectVarsRF.',jobid,'.Rdata'))
  save(x,y,ntree,nperms,nruns,savedir,alpha,mtry, jobid, file = DATAfile)
  
  shcmd = paste0('for run in {1..',nruns+5,'} ; do
                qsub -q all.q -l h_vmem=', memory, ',s_vmem=', memory, ' -j y -o /dev/null -e /dev/null -m e -M albnet@gmail.com \\
                 -v DATAfile=\'', DATAfile, '\' /data/jag/dpustina/Code/APHASIA/selectVarsRF.singlejob.sh \
                 done')
  rescmd = paste0('Selection = selectVarsRF.result(\'',jobid,'\')')
  if (!noprint) {
    cat('\nPLEASE RUN THE FOLLOWING COMMAND FROM Chead\n')
    cat(shcmd)
    cat('\nTO COMPUTE THE RESULTS RUN IN R:\n')
    cat(rescmd)
  } else {
    return(list(shell=shcmd,result=rescmd))
  }

}


# this function returns the results
selectVarsRF.result <- function(jobid, cleanup=F, savedir='/data/jag/dpustina/TEMP', noprint=T) {
  
  DATAfile = file.path(savedir,paste0('selectVarsRF.',jobid,'.Rdata'))
  if (!file.exists(DATAfile)) stop(paste('File',DATAfile,'does not exist'))
  load(DATAfile)
  
  resfiles = Sys.glob(file.path(savedir,paste0('selectVarsRF.RESULT.result.',jobid, '.*.csv')))
  if (length(resfiles) < nruns) {
    cat(paste0('\nProcessing is not finished: ', length(resfiles),'/', nruns, ' done.'))
    return(NULL)
  } else {
    results= lapply(resfiles,read.table)
    results = t(as.matrix(do.call(cbind, results)))
    colnames(results) = colnames(x)
    sigfiles = Sys.glob(file.path(savedir,paste0('selectVarsRF.RESULT.significant.',jobid, '.*.csv')))
    significants= lapply(sigfiles,read.table)
    significants = t(as.matrix(do.call(cbind, significants))) 
    sigs = colSums(significants)
    
    thisorder = order(sigs, decreasing= T)
    maxsig = rowSums(significants)
    avgsig = mean(maxsig)
    MYpredicts = data.frame(AverageOrder = colMeans(results)[thisorder] , SignificanceRatio=(sigs/nruns)[thisorder])
    selAvgSig = thisorder[1:round(avgsig)]
    onlytopsigs = MYpredicts$SignificanceRatio[1:avgsig] > 0.5
    selTopSig = selAvgSig[onlytopsigs]

    if (!noprint) {
      cat(paste('\nAverage significant variables:', avgsig, '(', nperms,'permutes,', nruns,'runs)\n'))
      cat(colnames(results)[selAvgSig])
    }
    
    if (cleanup) {
      rmfiles = Sys.glob(file.path(savedir,paste0('*',jobid,'*')) )
      file.remove(rmfiles)
    }

    output = list(MYpredicts=MYpredicts,
                  results=results,
                  significants=significants, 
                  sigs=sigs, 
                  maxsig=maxsig, 
                  AvgSig=avgsig,
                  selAvgSig=selAvgSig, 
                  selTopSig=selTopSig)
    return(output)
  }
    
  
}
