source('/data/jag/dpustina/Code/APHASIA/loadExcel.R', echo=F)

restppl =  SD$rest1==1 & !is.na(SD$rest1)
SD = SD[ restppl, ]
restnames = SD$NameOnDisk

iitT1 = antsImageRead('/data/jag/dpustina/IITmrtrix/IITmean_t1.nii.gz')
iitnodes = '/data/jag/dpustina/IITmrtrix/shenOnIITfromShen_GMonly.nii.gz'

iit2pennmat = '/data/jag/dpustina/IITmrtrix/ToPennTemplate/ToPennTemplate_0GenericAffine.mat'
iit2pennwarpinv = '/data/jag/dpustina/IITmrtrix/ToPennTemplate/ToPennTemplate_1InverseWarp.nii.gz'

# this part to find the native lesion mask
basedir = '/data/jag/dpustina/APHASIA/OluData/Niftii/'
datevar = read.table('/data/jag/dpustina/APHASIA/restsubs.txt')[1:56,]
rownames(datevar) = datevar$V1

# output base
outbase = '/data/jag/dpustina/APHASIA/STROKE/processing/IITconnectome/'

for (i in 1:length(restnames)) {
  cat(paste(i,''))
  
  mr = restnames[i]
  # find matrices to pennTemplate
  penn2submat = Sys.glob(paste0('/data/jag/dpustina/APHASIA/STROKE/processing/LesionMaskedRegistration/',
                       mr, '/pennTemplate_to_*_0GenericAffine.mat'))
  penn2subwarpinv = Sys.glob(paste0('/data/jag/dpustina/APHASIA/STROKE/processing/LesionMaskedRegistration/',
                                mr, '/pennTemplate_to_*_1InverseWarp.nii.gz'))
  if (any(!file.exists(c(penn2submat, penn2subwarpinv)))) stop(paste('Cannot find files for', mr))
  
  # concatenate
  conmat = c(iit2pennmat, iit2pennwarpinv, penn2submat, penn2subwarpinv) # concatenate matrices
  whichtoinvert = c(1,0,1,0)
  
  # find lesion
  cmd = paste0('ls ', basedir, mr, ' | sort ', ifelse(datevar[mr,]$V2==2,'-r ','') , '| head -n 1')
  thisdate = system(cmd, intern=T)
    lfile=Sys.glob(file.path(basedir,mr,thisdate,'rawNii','myMaskNative.nii.gz'))
    if (!file.exists(lfile)) stop(paste('Cannot find lesion for', mr))
    lesion = antsImageRead(lfile)
  t1file=Sys.glob(file.path(basedir,mr,thisdate,'rawNii','*MPRAGE.nii.gz'))
  if (!file.exists(t1file)) stop(paste('Cannot find t1 for', mr))
  mprage = antsImageRead(t1file[1])
  
  # reset the header to mprage in case...
  if ( any(antsGetOrigin(lesion) != antsGetOrigin(mprage)) ) lesion = antsCopyImageInfo(mprage,lesion)
  if ( any(antsGetDirection(lesion) != antsGetDirection(mprage)) ) lesion = antsCopyImageInfo(mprage,lesion)
  
  # bring lesion to iit
  suboutdir = file.path(outbase,mr)
  dir.create(suboutdir,showWarnings = F)
  iitlesfile = file.path(suboutdir,paste0(mr,'_Lesion_on_IIT.nii.gz'))
  iitt1file = file.path(suboutdir,paste0(mr,'_T1_on_IIT.nii.gz'))
#   if ( any(!file.exists(c(iitlesfile,iitt1file))) ) {
    t1iit = antsApplyTransforms(fixed = iitT1, moving = mprage, transformlist = conmat,
                                 whichtoinvert=whichtoinvert, interpolator = 'Linear')
    lesiit = antsApplyTransforms(fixed = iitT1, moving = lesion, transformlist = conmat,
                                 whichtoinvert=whichtoinvert, interpolator = 'NearestNeighbor')
    
    
    # save warped in subject's iit folder
    antsImageWrite(t1iit, iitt1file)
    antsImageWrite(lesiit, iitlesfile )
#   } else {
#     lesiit = antsImageRead(iitlesfile)
#     t1iit = antsImageRead(iitt1file)
#   }
  
  # save before and after overlap
  plot(mprage,lesion, outname=file.path(suboutdir, 'qa_before.png')); 
  plot(t1iit,lesiit, outname=file.path(suboutdir, 'qa_after.png')); 


  
  # tckedit to eliminate lesion streamlines
  intck = '/data/jag/dpustina/IITmrtrix/70M_10M_SIFTStep1mm-2mmHardi.tck'
  outtck = file.path(suboutdir, paste0(mr,'_lesionedTractography.tck'))
  cmd = paste('tckedit -nthread 0 -exclude ', iitlesfile, intck, outtck)
  cmd = paste('qsub -q all.q -l h_vmem=1.51G,s_vmem=1.5G -cwd -j y -o /dev/null -e /dev/null -m e -M albnet@gmail.com -b y', cmd)
  write(cmd, file='/data/jag/dpustina/APHASIA/STROKE/processing/IITconnectome/tckedit.sh', append=T)
  
  # build connectom matrix
  outmatrix = file.path(suboutdir, paste0(mr,'_FinnConnectome.csv') )
  cmd = paste('tck2connectome -nthread 0 -assignment_radial_search 15', outtck, iitnodes, outmatrix)
  cmd = paste('qsub -q all.q -l h_vmem=1.51G,s_vmem=1.5G -cwd -j y -o /dev/null -e /dev/null -m e -M albnet@gmail.com -b y', cmd)
  write(cmd, file='/data/jag/dpustina/APHASIA/STROKE/processing/IITconnectome/tck2connectome.sh', append=T)
  
}






# # write pheatmap and normalized connectomes for everyone
# outbase = '/data/jag/dpustina/APHASIA/STROKE/processing/IITconnectome/'
# fullcsv = read.table( '/data/jag/dpustina/IITmrtrix/FinnConnectome.csv', sep=' ')[,1:268]
# fullcsv=sna::symmetrize(fullcsv, 'upper')
# logfullcsv = fullcsv; logfullcsv[fullcsv>0] = log(fullcsv[fullcsv>0])
# maxlogfullcsv = max(logfullcsv[upper.tri(logfullcsv)])
# normfullcsv=logfullcsv/maxlogfullcsv
# normfullcsv[normfullcsv>1] = 1
# pheatmap(normfullcsv,  filename='/data/jag/dpustina/IITmrtrix/FinnConnectome.png', cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F, legend=F)
# 
# require(pheatmap)
# for ( i in 1:length(restnames) ) {
#   cat(paste(i,''))
#   
#   mr = restnames[i]
#   suboutdir = file.path(outbase,mr)
#   csv = read.table( file.path(suboutdir, paste0(mr,'_FinnConnectome.csv')), 
#                     sep=' ')[,1:268]
#   csv=sna::symmetrize(csv, 'upper')
#   
#   pheatmap(csv, filename=file.path(suboutdir,paste0(mr,'_rawMatrix.png')), cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F, legend=F)
#   logcsv=csv;logcsv[csv>0]=log(csv[csv>0])
#   max = max(logcsv[upper.tri(logcsv)])
#   normalcsv = logcsv/maxlogfullcsv
#   normalcsv[normalcsv>1]=1 # in case diagonal has higher values, clip to 1
#   pheatmap(normalcsv,  filename=file.path(suboutdir,paste0(mr,'_logNormalizedMatrix.png')), cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F, legend=F)
#   write.table(normalcsv, file=file.path(suboutdir,paste0(mr,'_NormalizedFinnConnectome.csv')) , row.names = F, col.names = F)
#   
#   lesioncsv = fullcsv-csv
#   lesioncsv[lesioncsv>0] = log(lesioncsv[lesioncsv>0])
#   lesioncsv = lesioncsv/maxlogfullcsv
#   lesioncsv[1,1]=1
#   lesioncsv[lesioncsv>1]=1
#   pheatmap(lesioncsv,  filename=file.path(suboutdir,paste0(mr,'_logLesionMatrix.png')), cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F, legend=F)
# 
# }





# # write tckinfo for everyone
# outbase = '/data/jag/dpustina/APHASIA/STROKE/processing/IITconnectome/'
# for (mr in restnames) {
#   suboutdir = file.path(outbase,mr)
#   cmd = paste('/share/apps/MRtrix3/2016-04-25/mrtrix3/release/bin/tckinfo', file.path(suboutdir, paste0(mr,'_lesionedTractography.tck')) )
#   info = system(cmd, intern=T)
#   write(info, file.path(suboutdir, 'trackinfo.txt') )
#   cat(paste(mr,''))
# }






# check all have tck file
# for (mr in restnames) {
#   len = length(
#     Sys.glob(file.path(outbase, mr,'*lesionedTractoghraphy.tck'))
#   )
#   if (len == 0) print(paste('Missing', mr))
# }