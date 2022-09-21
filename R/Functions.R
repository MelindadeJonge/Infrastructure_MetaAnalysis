
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


## ---- Calculate marginal and conditional R2 ------------------------------------------------------
R2.func <- function(model) {
  fixvar <-model$QM / model$k
  R2m<-fixvar / (fixvar + sum(model$sigma2))
  R2c<-(sum(model$sigma2) - model$sigma2[1] + fixvar) / (fixvar+sum(model$sigma2))
  x <- c(R2m,R2c)
  return(x)
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


distancePlots <- function(model, data, title) {
  science_theme = theme(panel.grid.major = element_blank(),
                        axis.title.y = element_text(vjust=1.25, size=8),
                        axis.title.x = element_text(vjust=-0.3, size=8),
                        axis.line = element_line(size=0.7, color="black"),
                        axis.text.y = element_text(size=8),
                        axis.text.x = element_text(size=8),
                        legend.title = element_blank(),
                        legend.position = 'right',
                        legend.direct = 'vertical',
                        plot.title = element_text(size=8))

  # Get the marginal and conditional R2 of the model
  R2 = round(R2.func(model), 2)

  # Create predictions with 95% confidence intervalls
  predDistance = seq(min(data$logD), max(data$logD), by=0.05)
  XNew = data.frame(logD=(predDistance),
                    logD2=(predDistance)^2)
  p = predict(model,
              newmods=as.matrix(XNew),
              addx=TRUE, level=95)
  predictions = data.frame(preds=p$pred,
                           ci.ub=p$ci.ub,
                           ci.lb=p$ci.lb,
                           Dist=predDistance)
  # Make plot
  xmax = max(data$logD)
  limy = max(c(abs(min(data$logRR)),abs(min(predictions$preds)),max(data$logRR),max(predictions$preds)))
  data$pointsize1 = log10(1/data$VAR) - min(log10(1/data$VAR))
  distPlot <- ggplot(data) +
    ggtitle(title) +
    geom_point(aes(logD,logRR), shape=20, color='slateblue4', alpha=0.1, size=data$pointsize1) +
    ylab(expression(LRR^Delta)) +
    xlab("Distance to infrastructure (m)") +
    scale_x_continuous(breaks=c(-1,0,1,2,3), labels=c('0.1','1','10','100','1000'), limits=c(0,xmax)) +
    geom_hline(yintercept=0, size=0.8, linetype="dashed", color="dark gray") +
    geom_ribbon(data=predictions, aes(x=Dist, ymin=ci.lb, ymax=ci.ub), fill="grey50",linetype="blank", alpha=0.6) +
    geom_line(data=predictions, aes(Dist,preds)) +
    scale_y_continuous(breaks=c(-7.5,-5,-2.5,0,2.5,5,7.5), lim=c(-limy,limy)) +
    annotate("text", x=0.60 * xmax, y=-limy + limy / 3.5,
             size=8/3, hjust=0, vjust=0.25,
             label=as.expression(bquote(italic(R[m]^2) == .(R2[1]))),
             parse=TRUE) +
    annotate("text", x=0.60 * xmax, y=-limy,
             size=8/3, hjust=0, vjust=0.25,
             label=as.expression(bquote(italic(R[c]^2) == .(R2[2]))),
             parse=TRUE) +
    theme_bw() +
    science_theme
  return(distPlot)
}

