# code copied from densitySearch_allThreshoolds.R

library(ggplot2)
library(caret)
levorder=c(
  'DTI_deg', 'DTI_bwn', 'DTI_lot', 'DTI_eff','DTI_mat', 
  'REST_deg', 'REST_bwn', 'REST_lot', 'REST_eff','REST_mat',
  'Les_Size', 'Parc_Damage',
  'Final_All', 'Final_RFE')


behavs = c('PNTcorrect', 'WABAQ', 'WABrep', 'WABcomp')

for (i in 1:length(behavs)) {

  load(paste0('/data/tesla-data/dpustina/MOVE-ME/APHASIA/STROKE/analyses/papier/',behavs[i],'_finaldata.Rdata'))
  
  selvars.x=rfefinal$optVariables
  #selvars.y=rep(max(finaldata$cor),length(selvars.x))
  selvars.x = gsub('DTI','DTI_', selvars.x)
  selvars.x = gsub('REST','REST_', selvars.x)
  selvars.x = gsub('DMG','Parc_Damage', selvars.x)
  selvars.x = gsub('LesionSize','Les_Size', selvars.x)
  
    
  ggplot(finaldata, aes(factor(name), cor)) + 
    geom_boxplot() + # ylim(0.2,0.95) +
    scale_y_continuous(breaks=seq(0.2,0.9,0.1), limits=c(0.2,0.95), minor_breaks = NULL) +
    ylab('Correlation') +
    ggtitle(behavs[i]) +
    theme_gray() +
    guides(fill=F) +
    annotate('text',x=selvars.x, y=rep(0.95, length(selvars.x)),
             label='*', color='gray36', size=12) +
    annotate('segment', x=12.5, xend=12.5, y=0.2,yend=0.95, color='red',alpha=0.4,lwd=1.5,lty=1) + 
    annotate('segment', x=14.5, xend=14.5, y=0.2,yend=0.95, color='red4',alpha=0.4,lwd=1.5,lty=1) + # annotate("rect", xmin=12.5, xmax=14.5, ymin=0.2, ymax=0.95, alpha=0.2) +
    theme(
      axis.text.x=element_text(angle=45,hjust=1,size=16),
      axis.text.y=element_text(size=14),
      axis.title.x=element_blank(),
      axis.title.y=element_text(size=18),
      plot.title=element_text(size=20))
  
  ggsave(width=8, height=6.5, 
         filename=file.path('/data/jag/dpustina/APHASIA/STROKE/analyses/papier',paste0(behavs[i],'_SMAPplot_allThresh_BW.png') ))
}