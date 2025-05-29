

GROUP.PROFILES <- function(data, groups, variables,  
                           plot_type ='bar', bar_type = 'all',
                           rescale='standardize',
                           CI_level= 95, ylim=NULL,
                           verbose=TRUE) {
  
  donnes <- data  
  
  donnes <- donnes[,c(groups,variables)]
  donnes <- donnes[order(donnes[,1]),]
  
  # group names, in the same order as in the data matrix
  grpnames <- as.vector(as.matrix(na.omit(donnes[,groups]))) 
  grpnames <- unique(grpnames)
  grpnums  <- seq(1:length(grpnames))
  Ngroups  <- length(grpnames)
  Nvars    <- length(variables)
  
  # the critical z value (that corresponds to the specified CI) to be used in the CI computations
  zforCI <- qnorm((1 + CI_level * .01) / 2) 
  
  # create plot data
  MNs <- SDs <- Ns <- matrix(-9999, Ngroups, Nvars)
  
  for (lupe in 1:Nvars) {
    dontemp <- data.frame(na.omit(cbind(donnes[,groups], donnes[,variables[lupe]])))
    
    # standardize DV
    if (rescale == 'standardize')  dontemp[,2] <- scale(dontemp[,2])
    
    # use specified, or default min & max range
    if (rescale == 'data') {
      
      dontemp[,2] <- RECODE(dontemp[,2], type = 'new_range',
                            real_min = min(dontemp[,2]), real_max = max(dontemp[,2]),
                            new_min  = min(donnes[,variables]), new_max = max(donnes[,variables]) )
    }
    
    MNs[,lupe] <- unlist(aggregate(dontemp$X2, list(dontemp$X1), FUN=mean)[2])
    
    SDs[,lupe] <- unlist(aggregate(dontemp$X2, list(dontemp$X1), FUN=sd)[2])
    
    Ns[,lupe]  <- unlist(aggregate(dontemp$X2, list(dontemp$X1), FUN=length)[2])
  }
  
  colnames(MNs) <-colnames(SDs) <- colnames(Ns) <- variables
  rownames(MNs) <-rownames(SDs) <- rownames(Ns) <- grpnames
  
  CIs_ub <- MNs + zforCI * (SDs / sqrt(Ns))
  CIs_lb <- MNs - zforCI * (SDs / sqrt(Ns))
  
  
  # arranging output by predictor
  output <- list()
  for (lupe in 1:ncol(MNs)) {
    temp <- cbind(MNs[,lupe], SDs[,lupe], Ns[,lupe], CIs_lb[,lupe], CIs_ub[,lupe])
    colnames(temp) <- c('Mean','SD','N','CI_lb','CI_ub')
    rownames(temp) <- grpnames
    output[[colnames(MNs)[lupe]]] <- temp
  }  
  
  # plots
  cols <- rainbow(Nvars)
  
  
  if (plot_type == 'bar') {	
    
    if (bar_type == 'all') {
      
      if (is.null(ylim)) ylim <- range( pretty(CIs_lb), pretty(CIs_ub))
      
      ylim[2] <- ylim[2] + (ylim[2] * .05)
      
      #	  barplot(1, ylim=ylim, ylab="DV")
      
      plot_bar <- barplot(t(MNs), beside=TRUE, legend.text=FALSE, col=cols, 
                          ylim=ylim, ylab="DV", yaxt="n", xpd=FALSE) #, axes=FALSE) 
      
      yyy = seq(ylim[1], ylim[2], by=.2)
      
      axis(side=2, at=yyy, labels=yyy, las=1)
      
      arrows(plot_bar, t(CIs_ub), plot_bar, t(CIs_lb), 
             lwd = 1.0, angle = 90, code = 3, length = 0.05)
      
      graphics::legend("top", legend = variables, bty="n", lwd=2, col=cols, cex = .80)
    }
    
    if(bar_type == 'separate') {
      
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
      
      # par(mfrow=c(2,1), pty="m",  oma = c(2,4,0,0) + 0.1, mar = c(2,2,1,1) + 2.6)   # 1.8
      
      if (Nvars == 1) par(mfrow = c(1, 1))
      if (Nvars == 2) par(mfrow = c(2, 1),  oma = c(0,4,0,4) )#, mar = c(2,2,1,1))
      if (Nvars == 3) par(mfrow = c(2, 2))
      if (Nvars == 4) par(mfrow = c(2, 2))
      
      if (Nvars <= 4) Nplots <- Nvars
      if (Nvars >  4) Nplots <- 4
      
      for (lupe in 1:Nplots)  {
        
        ylim <- range( pretty(CIs_lb[,lupe]), pretty(CIs_ub[,lupe]))
        
        plot_bar <- barplot(t(MNs[,lupe]), beside=F, legend.text=FALSE, col=cols[lupe], 
                            ylim=ylim, ylab="DV", yaxt="n", xpd=FALSE) #, axes=FALSE)  # , xpd=FALSE)
        
        yyy = seq(ylim[1], ylim[2], by=.2)
        
        axis(side=2, at=yyy, labels=yyy, las=1)
        
        arrows(plot_bar, t(CIs_ub[,lupe]), plot_bar, t(CIs_lb[,lupe]), 
               lwd = 1.0, angle = 90, code = 3, length = 0.05)
        
        graphics::title(main= paste('Group Means for',variables[lupe]))
        
        eval(parse(text=(paste('pp', lupe, ' <- recordPlot()', sep=''))))
      }
      par(mfrow = c(1,1))
      # dev.off()
      # graphics.off()
    }
  }
  
  if (plot_type == 'profile') {	
    
    # png("Figure # - group profiles plot.tiff", width=6, height=6, units="in", res = 600, pointsize=12)
    # 
    # par(font.main=1, font.lab=1, font.axis=1, cex=1, cex.main=1, cex.lab=1, cex.axis=1,
    #     lwd=2, las=1,  pty="s",  mai=c(.8,.8,.8,.8) ) 
    
    if (is.null(ylim))  ylim <- range( min(pretty(min(MNs))), max(pretty(max(MNs))) )
    
    graphics::matplot(1:Ngroups, MNs, type = "l", lty=1, lwd=3, 
                      xaxt='n', xlab='', cex.axis=1.2, cex.lab = 1.3,
                      ylab='DV Scores', ylim = ylim, cex.axis=1.2, col=cols )
    
    graphics::axis(side=1, at=grpnums, labels=grpnames, xlab="groups")   # labels=rownames(centroidsZ_2)
    
    graphics::title(main='Profiles of the Group Means')
    
    graphics::legend("top", legend = colnames(MNs), bty="n", lwd=2, col=cols)
    
    # text(seq(1, length(grpnums), by=1), par("usr")[3] - 0.2,
    #      labels = grpnames,pos = 1, xpd = TRUE)  #  srt = 90, 
    
    #	  dev.off()
  }
  
  
  if (verbose) {	
    
    if (verbose) {
      
      message('\nThe confidence interval for the output: ', CI_level, '%')
      
      message('\nThe rescale variables argument was set to: ', rescale)
      
      for (lupe in 1:length(output)) {
        message('\n', colnames(MNs)[lupe])
        print(round(output[[lupe]],2), print.gap=4)
      }
    }
    
    # message("\n\nGroup Means:\n") 
    # print(round(MNs,2))
    # 
    # message("\n\nGroup Standard Deviations:\n") 
    # print(round(SDs,2))
    # 
    # message("\n\nGroup Sizes:\n") 
    # print(round(Ns,2))
    # 
    # message("\n\nConfidence Intervals - lower bounds:\n") 
    # print(round(CIs_lb,2))
    # 
    # message("\n\nConfidence Intervals - upper bounds:\n") 
    # print(round(CIs_ub,2))
    
  }
}




