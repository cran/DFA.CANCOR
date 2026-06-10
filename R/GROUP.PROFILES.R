

GROUP.PROFILES <- function(data, groups, variables,  
                           plot_type ='bar', bar_type = 'all',
                           rescale='standardize',
                           CI_level= 95, ylim=NULL,
                           plot_save = FALSE, plot_save_type = 'png',
                           plot_title = NULL,
                           cols_user = NULL,
                           verbose=TRUE) {
  
  donnes <- data  
  
  donnes <- donnes[,c(groups,variables)]
  donnes <- donnes[order(donnes[[groups]]),]

  # group names, in the same order as in the data matrix
  grpnames <- as.vector(as.matrix(na.omit(donnes[,groups]))) 
  grpnames <- unique(grpnames)
  grpnums  <- seq(1:length(grpnames))
  Ngroups  <- length(grpnames)
  Nvars    <- length(variables)
  
  # the critical z value (that corresponds to the specified CI) to be used in the CI computations
  zforCI <- qnorm((1 + CI_level * .01) / 2) 
  
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  if (plot_save == TRUE) {
    
    height=6; width=5
    
    if (is.null(plot_title)) plot_title = plot_type
    
    if (is.null(plot_save_type))  plot_save_type = 'png'
    
    if (plot_save_type == 'bitmap')
      bitmap(paste('Figure - ',plot_title,'.bitmap',sep=''), height=height, width=width, units='in', res=1200, pointsize=12)
    
    if (plot_save_type == 'tiff')
      tiff(paste('Figure - ',plot_title,'.tiff',sep=''), height=height, width=width, units='in', res=1200, pointsize=12)
    
    if (plot_save_type == 'png')
      png(paste('Figure - ',plot_title,'.png',sep=''), height=height, width=width, units='in', res=1200, pointsize=12)
    
    if (plot_save_type == 'jpeg')
      jpeg(paste('Figure - ',plot_title,'.jpeg',sep=''), height=height, width=width, units='in', res=1200, pointsize=12)
    
    if (plot_save_type == 'bmp')
      bmp(paste('Figure - ',plot_title,'.bmp',sep=''), height=height, width=width, units='in', res=1200, pointsize=12)
  }
  
  
  # rescale
  if (rescale == 'standardize')  donnes[, variables] <- scale(donnes[, variables])

  desdat <- DESCRIPTIVES(data=donnes, groups=groups, variables=variables, verbose=FALSE)
  
  MNs <- CIs_lb <- CIs_ub <- ymin <- ymax <- c()
  for (lupe in 1:length(desdat$DESCR_vars)) {
    MNs    <- cbind(MNs,    (desdat$DESCR_vars[[lupe]][,'Mean'] ))
    CIs_lb <- cbind(CIs_lb, (desdat$DESCR_vars[[lupe]][,'CI_lb'] ))
    CIs_ub <- cbind(CIs_ub, (desdat$DESCR_vars[[lupe]][,'CI_ub'] ))
    ymin   <- cbind(ymin,   (desdat$DESCR_vars[[lupe]][,'min'] ))
    ymax   <- cbind(ymax,   (desdat$DESCR_vars[[lupe]][,'max'] ))
  }
  colnames(MNs) <-colnames(CIs_lb) <- colnames(CIs_ub) <- variables
  rownames(MNs) <-rownames(CIs_lb) <- rownames(CIs_ub) <- grpnames
  

  # plots
  if (is.null(cols_user)) {
    cols <- c('blue', 'red', 'cyan2', 'darkviolet', 'chartreuse1', 'yellow',
              'burlywood3','darkseagreen1', 
              'mediumvioletred', 'darkgreen','bisque','cyan3', 'deeppink4')
  } else {cols <- cols_user}
  
  if (plot_type == 'bar') {	
    
    if (bar_type == 'all') {
      
      if (is.null(ylim)) ylim <- c(min(ymin), max(ymax))
      
      ylim <- c(min(MNs), max(MNs))
      
      ylim <- better_ylim(ylim, buffer = 0.15)
      
      # if (is.null(ylim)) ylim <- range( pretty(CIs_lb), pretty(CIs_ub))
      
      plot_bar <- barplot(t(MNs), beside=TRUE, legend.text=FALSE, col=cols[1:Nvars], 
                          ylim=ylim, ylab='DV', yaxt='n', main=plot_title, xpd=FALSE) #, axes=FALSE) 
      
      yyy = round( seq(ylim[1], ylim[2], by=.2), 2)
      
      axis(side=2, at=yyy, labels=yyy, las=1)
      
      if ( (min(CIs_lb) > ylim[1]) & (max(CIs_lb) < ylim[2]) ) {
        suppressWarnings(arrows(plot_bar, t(CIs_ub), plot_bar, t(CIs_lb), 
               lwd = 1.0, angle = 90, code = 3, length = 0.05))
      }
      
    graphics::legend('top', legend = variables, bty='n', lwd=2, col=cols, cex = .80)
  }
  
    if(bar_type == 'separate') {
      
      # if (nrow(activeSet_ALL) >  1)  par(mfrow=c(1,1), pty='m', mar=c(3,2,3,2) + 2.6, ask=TRUE)
      # if (nrow(activeSet_ALL) == 1)  par(mfrow=c(1,1), pty='m', mar=c(3,2,3,2) + 2.6, ask=FALSE)
      
      # oldpar <- par(no.readonly = TRUE)
      # on.exit(par(oldpar))
      
      # par(mfrow=c(2,1), pty='m',  oma = c(2,4,0,0) + 0.1, mar = c(2,2,1,1) + 2.6)   # 1.8
      
      # if (Nvars <= 4) Nplots <- Nvars
      # if (Nvars >  4) Nplots <- 4
      
      if (Nvars == 1) par(mfrow = c(1, 1))
      if (Nvars == 2) par(mfrow = c(2, 1),  oma = c(0,4,0,4) )#, mar = c(2,2,1,1))
      if (Nvars == 3) par(mfrow = c(2, 2))
      if (Nvars == 4) par(mfrow = c(2, 2))
      if (Nvars >  4) par(mfrow = c(1, 1), ask=TRUE)
      
      Nplots <- Nvars
      
      for (lupe in 1:Nplots)  {
        
        if (is.null(ylim)) ylim <- c(min(ymin), max(ymax))
        
        ylim <- c(min(MNs), max(MNs))
        
        ylim <- better_ylim(ylim, buffer = 0.15)
        
        # if (is.null(ylim)) ylim <- range( pretty(CIs_lb), pretty(CIs_ub))
        
        plot_bar <- barplot(t(MNs[,lupe]), beside=T, legend.text=FALSE, col=cols[1:Ngroups], 
                            ylim=ylim, ylab='DV', yaxt='n', xpd=FALSE) #, axes=FALSE) 
        
        yyy = round( seq(ylim[1], ylim[2], by=.2), 2)
        
        axis(side=2, at=yyy, labels=yyy, las=1)
        
        if ( (min(CIs_lb[,lupe]) > ylim[1]) & (max(CIs_lb[,lupe]) < ylim[2]) ) {
          suppressWarnings(arrows(plot_bar, t(CIs_ub[,lupe]), plot_bar, t(CIs_lb[,lupe]), 
                 lwd = 1.0, angle = 90, code = 3, length = 0.05))
        }

        graphics::title(main= paste('Group Means for',variables[lupe]))
        
        eval(parse(text=(paste('pp', lupe, ' <- recordPlot()', sep=''))))
      }
      # par(mfrow = c(1,1))
      # dev.off()
      # graphics.off()
    }
  }
  
  if (plot_type == 'profile') {	
    
    # png('Figure # - group profiles plot.tiff', width=6, height=6, units='in', res = 600, pointsize=12)
    # 
    # par(font.main=1, font.lab=1, font.axis=1, cex=1, cex.main=1, cex.lab=1, cex.axis=1,
    #     lwd=2, las=1,  pty='s',  mai=c(.8,.8,.8,.8) ) 
    
    # if (is.null(ylim))  ylim <- range( min(pretty(min(MNs))), max(pretty(max(MNs))) )
    if (is.null(ylim))  ylim <- range( min(MNs), max(MNs) )
    
    graphics::matplot(1:Ngroups, MNs, type = 'l', lty=1, lwd=3, 
                      xaxt='n', xlab='', cex.axis=1.2, cex.lab = 1.3,
                      ylab='DV Scores', ylim = ylim, cex.axis=1.2, col=cols )
    
    graphics::axis(side=1, at=grpnums, labels=grpnames, xlab='groups')   # labels=rownames(centroidsZ_2)
    
    graphics::title(main='Profiles of the Group Means')
    
    graphics::legend('top', legend = colnames(MNs), bty='n', lwd=2, col=cols)
    
    # text(seq(1, length(grpnums), by=1), par('usr')[3] - 0.2,
    #      labels = grpnames,pos = 1, xpd = TRUE)  #  srt = 90, 
    
    #	  dev.off()
  }
  
  
  if (verbose) {	
    
    if (verbose) {
      
      cat('\n\nThe confidence interval for the output: ', CI_level, '%', sep='')
      
      cat('\n\nThe rescale argument was set to:', rescale, '\n')
      
      for (lupe in 1:length(desdat$DESCR_vars)) {
        cat('\n\n', colnames(MNs)[lupe], '\n')
        print(round(desdat$DESCR_vars[[lupe]],2), print.gap=4)
      }
    }
    
    # message('\n\nGroup Means:\n') 
    # print(round(MNs,2))
    # 
    # message('\n\nGroup Standard Deviations:\n') 
    # print(round(SDs,2))
    # 
    # message('\n\nGroup Sizes:\n') 
    # print(round(Ns,2))
    # 
    # message('\n\nConfidence Intervals - lower bounds:\n') 
    # print(round(CIs_lb,2))
    # 
    # message('\n\nConfidence Intervals - upper bounds:\n') 
    # print(round(CIs_ub,2))
    
  }
  
  if (plot_save == TRUE)  dev.off()
  
  return(invisible(desdat))
}




