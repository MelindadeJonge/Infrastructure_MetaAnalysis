
## ---- Variance-covariance matrix -----------------------------------------------------------------
calc.vRRDelta <- function(x) {
  v = matrix((x$ControlSD_imputed[1]^2 / (x$ControlN[1] * x$ControlA[1]^2) + 
                0.5 * (x$ControlSD_imputed[1]^4 / (x$ControlN[1]^2 * x$ControlA[1]^4))), 
             nrow=nrow(x), ncol=nrow(x))
  
  # The following are some tests to check for errors/typos in the dataset
  if (length(unique(x$ControlSD_imputed)) > 1) {
    print('Control variances not equal')
    print(x$ControlID[1])
  }
  diag(v) = x$VAR
  if (!corpcor::is.positive.definite(v)) {
    print(x$ControlID[1])
    print(v)
  }
  if (length(unique(x$SpeciesName)) > 1) {
    print(x$ControlID[1])
    print(unique(x$SpeciesName))
  }
  if (length(unique(x$CDist)) > 1) {
    print(x$ControlID[1])
    print(unique(x$CDist))
  }
  return(v)
} 

## ---- Figures ------------------------------------------------------------------------------------
createMap <- function(X, title) {
  map_theme = theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill='white', colour='White'),
                    plot.background = element_rect(fill='white', colour='White'),
                    axis.title.y = element_blank(),
                    axis.title.x = element_blank(),
                    axis.line = element_blank(),
                    axis.text.y = element_blank(),
                    axis.text.x = element_blank(),
                    axis.ticks = element_blank(),
                    legend.position="bottom",
                    legend.box = "vertical",
                    legend.title = element_blank(),
                    legend.text = element_text(size=8),
                    plot.title = element_text(hjust=0, size=8))
  
  mP = aggregate(X$SpeciesName, 
                 by=list(X$Study, X$Latitude, X$Longitude, X$InfraType, X$SpeciesCode), 
                 length)
  
  mP$Group.3 = as.numeric(mP$Group.3)
  mP$Group.2 = as.numeric(mP$Group.2)
  mP$x = log(as.numeric(mP$x))
  mP$Group.4 = factor(mP$Group.4, levels=c('Non-traffic','Paved road','Powerline','Unpaved road'))
  mapWorld = borders("world", colour=NA, fill="gray75") # create a layer of country borders
  g1 = ggplot() +
    mapWorld +
    ggtitle(paste0('       ',title)) +
    geom_point(data=mP,aes(x=Group.3, y=Group.2,
                           colour=Group.4, size=x),
               shape=21, alpha=0.5) +
    scale_size(range=c(1,2.5), guide='none') +
    scale_color_manual(breaks=c('Non-traffic','Paved road','Powerline','Unpaved road'), values=viridis(4),
                       labels=c('Non-traffic','Paved road','Power line','Unpaved road'), drop=FALSE)+
    xlab(NULL) + ylab(NULL) + ylim(c(-55,84))+
    map_theme
  
  return(g1)
}

