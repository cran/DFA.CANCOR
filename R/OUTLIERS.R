



OUTLIERS <- function(data, variables, ID=NULL, iterate=TRUE,
                        alpha_univ=.05, plot_univariates=TRUE,
                        MCD=TRUE, MCD.quantile = .75, alpha=0.025, cutoff_type = 'adjusted',
                        qqplot=TRUE, plot_iters=NULL, 
                        verbose=TRUE) {

  if (is.null(ID))  ID <- 1:nrow(data)
  
  if (is.null(variables))  donnes <- as.matrix(data)
  
  if (!is.null(variables)) {
    
    donnes <- as.matrix(data[,variables])
    
    colnames(donnes) <- variables
  }
  
  rownames(donnes) <- ID
  
  if (anyNA(donnes)) {
    
    donnes <- na.omit(donnes)
    
    N_missing <- nrow(data) - nrow(donnes)
    
    message('\n\n', N_missing, 
            ' cases with missing values were found and removed from the data matrix.')
  }
  
  if (!is.numeric(donnes))
    message('\nThe data are not numeric, expect errors.\n\n')
  
  # if (!is.data.frame(data) && !is.matrix(data))
  #   stop("Input must be one of classes \"data frame\" or \"matrix\"")
  
  #  dataframe = as.data.frame(data)
  
  if (ncol(donnes) < 2 || is.null(dim(data))) 
    stop("The number of variables must be greater than or equal to 2")
  
  if (!iterate) max_iterations <- 1
  if ( iterate) max_iterations <- 15
  
  N_vars <- ncol(donnes)
  
  
  ###################################### univariate outliers #################################
  
  # z values
  zvalues <- scale(donnes);  #summary(zvalues)
  
  z_alpha <- qnorm(alpha_univ / 2)
  
  univ_outliers <- zvalues[apply(zvalues,1,function(x) any(abs(x) > abs(z_alpha))),] 

  if (verbose) {
    message('\n\nCases with z values beyond a p value of ', alpha_univ, ':\n')
    print(round_boc(univ_outliers,3))
  }
  
  if (plot_univariates) {
    
    opar <- par(no.readonly = TRUE)      # make a copy of current settings
  
    par(mfrow=c(2,2), pty="s", ask=TRUE, mar=c(3,2,3,2) + 2.0)  # , mgp=c(2,1,0)
    
    cex_legend = 0.7; cex_lab = .9; cex_axis = .7
    
    for (lupe in 1:N_vars) {
      
      plotdon_univ <- as.matrix(donnes[,lupe])
      colnames(plotdon_univ) <- colnames(donnes)[lupe]
      
      # oldpar <- par(no.readonly = TRUE)
      # on.exit(par(oldpar))
      
      
      # histogram with a normal curve
      histogram <- hist(plotdon_univ, breaks=20, col="red", 
                        xlab=colnames(plotdon_univ), main="Histogram", xaxt='n'
      ) 
      
      xat <- seq( floor(min(plotdon_univ)) ,ceiling(max(plotdon_univ)),1)
      
      axis(side=1, at=xat, labels=xat,
           cex.lab  = cex_lab,  # font size of axis labels
           cex.axis = cex_axis #font size of axis text 
      )
      
      xfit <- seq(min(plotdon_univ), max(plotdon_univ), length=40) 
      
      yfit <- dnorm(xfit, mean=mean(plotdon_univ), sd=sd(plotdon_univ)) 
      
      yfit <- yfit * diff(histogram$mids[1:2]) * length(plotdon_univ) 
      
      lines(xfit, yfit, col="blue", lwd=2)
      
      
      
      # box plot
      outlier_values <- boxplot.stats(plotdon_univ)$out 
      # outlier values. 
      boxplot(plotdon_univ, main="Box Plot", boxwex=1.1, outwex=3,
              # ylim=c(min(plotdon_univ), max(plotdon_univ)),
              ylab = colnames(plotdon_univ),
              col = "blue",
              border = "red",
              horizontal = FALSE,
              whisklty = 1,
              staplelwd = 4, outpch = 8,
              medcol = "red", boxlty = 0,
              ylim=c( floor(min(plotdon_univ))-1 ,ceiling(max(plotdon_univ))+1),
              cex.lab = cex_lab,  # font size of axis labels
              cex.axis = cex_axis #font size of axis text 
      ) 
      
      # mtext(paste("Outliers: ", paste(round(outlier_values,2), collapse=", ")), cex=0.6)
      
      # For a given continuous variable, outliers are those observations that lie 
      # outside 1.5 * IQR, where IQR, the Inter Quartile Range is the difference 
      # between 75th and 25th quartiles. The outliers are the points outside the whiskers in below box plot.
      
      
      
      # normal probability plot
      
      sorted <- sort(plotdon_univ, decreasing = FALSE)
      
      obscumprob <- cumsum(sorted)/ sum(sorted)
      
      ranked <- rank(sorted)
      
      N_plotdon_univ <- length(plotdon_univ)

      # Score Methods  "Help Online - Origin Help - Probability Plot and Q-Q Plot.pdf"
      
      # methods for plotting percentile approximations
      
      # These formulas produce noticeable differences for short series only
      
      # Input data is ordered from smallest to largest, and then the serial number of 
      # the sorted data is scored using one of the methods listed below. In this table,   
      # is the serial number and   is the total number of the nonmissing input data.
      
      # also used in SPSS PPLOT
            
      expcumprob   <- (ranked - 0.375) / (N_plotdon_univ + 0.25)  # Blom method
      
      # expcumprob <- (ranked - 0.3)   / (N_plotdon_univ + 0.4)  # Benard method
      # 
      # expcumprob <- (ranked - 0.5)   / (N_plotdon_univ)  # Hazen method
      # 
      # expcumprob <-  ranked          / (N_plotdon_univ + 1)  # Van der Waerden method
      # 
      # expcumprob <-  ranked          / (N_plotdon_univ)  # Kaplan-Meier method
      
      
      red_black <- ifelse(abs(zvalues) > abs(z_alpha), 'red', 'black')
      
      plot(
        obscumprob,
        expcumprob,
        main = "Normal P-P plot",
        xlab = "Observed cum prob.",
        ylab = "Expected cum prob.",
        col = red_black,
        cex.lab = cex_lab,  # font size of axis labels
        cex.axis = cex_axis #font size of axis text 
      )
      abline(0,1, col='red')
      
      
      
      # Q-Q (Quantile-Quantile) plot
      
      qqnorm(sorted,
             ylab=colnames(plotdon_univ),
             xlab="Normal Scores",
             main="Normal Q-Q plot",
             cex.lab = cex_lab,  # font size of axis labels
             cex.axis = cex_axis #font size of axis text 
      )
      qqline(sorted, col = "red")
    }
    
    par(opar)          # restore original settings
    # on.exit(par(opar))
    
  }
  
  
  ############################# multivariate outliers ##############################################
  
  if (N_vars > 1) {
    
    outliers_only_CUM <- c()
    
    plot_data_list <- list()
    
    for (lupe in 1:max_iterations) {
      
      N_cases <- nrow(donnes)
      
      # if (!is.null(rownames(donnes))) { row_noms <- rownames(donnes)
      # } else { row_noms <- 1:N_cases }
      
      if (!MCD) {
        
        centers <- colMeans(donnes)
        
        covmat  <- cov(donnes)
      }
      
      if (MCD) {
        
        if ( (N_cases * MCD.quantile) < 8) {
          message('\n There are not enough cases to run the requested analyses. Expect errors.\n')
          stop("Not enough cases to proceed")
        }
        
        outp_mcd <- MASS::cov.mcd(donnes)  #, quantile.used = N_cases * MCD.quantile) 
        
        centers <- outp_mcd$center
        
        covmat  <- outp_mcd$cov
      }
      
      mahal_dist <- mahalanobis(donnes, center = centers, cov = covmat)
      
      
      # covr <- covMcd(x, alpha=quan)
       
      # dist <- mahal_dist
      
      # s <- sort(dist, index=TRUE)
      # q <- (0.5:length(dist))/length(dist)
      # qchi <- qchisq(q, df=N_vars)
       
      # plot(s$x, qchi, xlab="Ordered robust MD^2", ylab="Quantiles of Chi_p^2", main="Chi^2-Plot", col=3)
      
      
      
      
      if (cutoff_type == 'quan')  cutoff <- qchisq(p = 1 - alpha, df = N_vars)
      
      if (cutoff_type == 'adjusted') 
        cutoff <- mvoutlier::arw(x = donnes, m0 = centers, c0 = covmat, alpha = 0.025)$cn
      
      all_cases <- cbind(data.frame(rownames(donnes)), donnes, mahal_dist, NA)
      
      colnames(all_cases) <- c("ID", variables, "Mahalanobis_Distance", "Outlier")
      
      all_cases$"Outlier" <- ifelse(all_cases$"Mahalanobis_Distance" > cutoff, TRUE, FALSE)
      
      outliers_only <- subset(all_cases, Outlier == TRUE)
      
      outliers_only_CUM <- rbind(outliers_only_CUM, outliers_only)
      
      if (nrow(outliers_only) == 0) { outlier_Ns <- c(N_cases, 0)
      } else { outlier_Ns <- c( (N_cases - nrow(outliers_only)), nrow(outliers_only)) }
      
      
      # create plot data
      mahal_dist_ranked <- rank(mahal_dist)
      
      chisq_quantile <- qchisq((mahal_dist_ranked - 0.5)/N_cases, N_vars)
      
      plot_data <- data.frame(mahal_dist, chisq_quantile, NA)
      
      colnames(plot_data) <- c("Mahalanobis_Distance", "chisq_quantile", "Outlier")
      
      plot_data$Outlier <- ifelse(plot_data$Mahalanobis_Distance > cutoff, TRUE, FALSE)
      
      #    plot_data[,"Outlier"] <- ifelse(plot_data[,"Mahalanobis_Distance"] > cutoff, TRUE, FALSE)
      
      # plot_data <- data.frame(mahal_dist_ranked, chisq, NA)
      
      # plot_data_list <- c(plot_data_list, list(cbind(plot_data, all_cases$"Outlier")))
      
      plot_data_list <- c(plot_data_list, list(plot_data))
      
      
      # remove the outliers from donnes & redo
      if (nrow(outliers_only) > 0)
        donnes <- donnes[!(row.names(donnes) %in% rownames(outliers_only)),]
      # donnes <- donnes[-as.numeric(unlist(outliers_only['ID'])),]
      
      
      if (verbose) {
        
        if (nrow(outliers_only) > 0) {
          
          message('\n\nThe multivariate outliers on iteration #', lupe, ':\n')
          
          print( outliers_only[order(outliers_only$Mahalanobis_Distance, decreasing = TRUE), ], row.names=FALSE )
          
          message('\nThere are ', outlier_Ns[2], ' outliers and ', outlier_Ns[1], ' non-outliers.')
          message('\nThe robust squared Mahalanobis distance cutoff is ', round(cutoff,3), '\n')
          
          # mvn tests
          message('Tests of multivariate normality after the outliers have been removed:\n')
          mvntests <- suppressWarnings(NORMALITY(data = as.data.frame(donnes), verbose=FALSE)$multivariate_tests)
          print(round_boc(mvntests,3), row.names=FALSE, print.gap=4)
        }
        
        if (nrow(outliers_only) == 0) 
          message('\nThere were no multivariate outliers on iteration #', lupe, '.\n\n')
      }
      
      if (nrow(outliers_only) == 0) break
    }
    
    
    if (qqplot) { 
      
      if (is.null(plot_iters)) plot_iters <- c(1,2,3,4)
      
      if (length(plot_iters) > 4) {
        
        plot_iters <- plot_iters[1:4]
        
        # message('"plot_iters" has been truncated.')
      }
      
      if (max(plot_iters) > length(plot_data_list))  plot_iters <- seq(1:length(plot_iters))
      
      if (length(plot_data_list) < length(plot_iters)) 
        plot_iters <- plot_iters[1:length(plot_data_list)]
      
      # oldpar <- par(no.readonly = TRUE)
      # on.exit(par(oldpar))

      opar <- par(no.readonly = TRUE)      # make a copy of current settings
      
            
      # if (plot_save) {
      #   
      #   if (is.null(plot_save_type))  plot_save_type = 'png'
      #   
      #   if (plot_save_type == 'bitmap')
      #     bitmap(paste("Figure - MV Outliers",".bitmap",sep=""), 
      #            height=7, width=9, units='in', res=1200, pointsize=12)
      #   
      #   if (plot_save_type == 'tiff')
      #     tiff(paste("Figure - MV Outliers",".tiff",sep=""), 
      #          height=8, width=8, units='in', res=1200, pointsize=12)
      #   
      #   if (plot_save_type == 'png')
      #     png(paste("Figure - MV Outliers",".png",sep=""), 
      #         height=8, width=8, units='in', res=1200, pointsize=12)
      #   
      #   if (plot_save_type == 'jpeg')
      #     jpeg(paste("Figure - MV Outliers",".jpeg",sep=""), 
      #          height=7, width=9, units='in', res=1200, pointsize=12)
      #   
      #   if (plot_save_type == 'bmp')
      #     bmp(paste("Figure - MV Outliers",".bmp",sep=""), 
      #         height=7, width=9, units='in', res=1200, pointsize=12)
      # }
      # 
      if (length(plot_iters) == 1)  {
        
        par(mfrow=c(1,1), pty="s", ask=TRUE) #, mar=c(3,2,3,2) + 2.6)
        
        cex_legend = 1; cex_lab = 1; cex_axis = 1
      }
      
      
      # if (length(plot_iters) == 2)  par(mfrow=c(2,1), pty="s") #, mar=c(3,2,3,2) + 2.6)
      
      # if (length(plot_iters) == 2)  par(mfrow=c(2,1), pty="m", mar=c(2,6,2,6))
      
      if (length(plot_iters) >  1)  {
        
        par(mfrow=c(2,2), pty="s", ask=TRUE, mar=c(3,2,3,2) + 2.0)  # , mgp=c(2,1,0)
        
        cex_legend = 0.7; cex_lab = .9; cex_axis = .7
      }
      
      
      # nf <- layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), respect = TRUE)
      # layout.show(nf)
      
      
      for (lupe in 1:length(plot_iters)) {
        
        plotdon <- plot_data_list[[plot_iters[lupe]]]
        
        # chisq <- qchisq((rank(mahal_dist) - 0.5)/N_cases, N_vars)
        # 
        # red_black <- ifelse(all_cases$Outlier == TRUE, 'red', 'black')
        
        red_black <- ifelse(plotdon[,3] == TRUE, 'red', 'black')
        
        xlim <- ylim <- c(0, max(plotdon[,'Mahalanobis_Distance'], plotdon[,'chisq_quantile']))
        
        plot(plotdon[,'Mahalanobis_Distance'], plotdon[,'chisq_quantile'], pch = 16, 
             main = paste('Iteration #', plot_iters[lupe], sep=""),
             #xlab = "Robust Squared Mahalanobis Distance",
             xlab = "Robust Squared Distance",
             # xlab = '',
             ylab = "Chi-Square Quantile", 
             col = red_black,
             xlim = xlim,
             ylim = ylim,
             cex.lab = cex_lab,   # font size of axis labels
             cex.axis = cex_axis) #font size of axis text )
        
        # title(xlab="Robust Squared Distance", line=2, cex.lab=1.2)
        
        abline(0,1, col='blue')
        
        outlier_Ns <- as.numeric(table(plotdon[,3]))
        
        if (length(outlier_Ns) == 1) outlier_Ns <- c(outlier_Ns, 0) 
        
        legend("topleft",
               legend = c(paste("Outliers (N = ", outlier_Ns[2], ")", sep = ""),
                          paste("Non (N = ", outlier_Ns[1],")", sep = "")), 
               col = c("red", "black"), pch = 16,  bty = "n", cex=cex_legend)
      }
      # if (plot_save)  dev.off()
      
      par(opar)          # restore original settings
      # dev.off()
      # on.exit(par(opar))
    }
  }    
  
  # if (plot_save)  dev.off()
  
  output <- list(outliers_only=outliers_only_CUM)
  
  return(invisible(output))
}



# 
# plot_save = FALSE, plot_save_type = 'png',
# 
# 
# \item{plot_save}{
#   \code{}Should a plot be saved to disk? TRUE or FALSE (the default).}
# 
# \item{plot_save_type}{
#   \code{}(optional) The output format if plot_save = TRUE. The options are 'bitmap', 'tiff', 
#   'png' (the default), 'jpeg', and 'bmp'.}
# 


