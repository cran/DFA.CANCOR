


GROUP.DIFFS <-  function (data, GROUPS = NULL, DV = NULL, var.equal = FALSE, 
                          p.adjust.method = "holm", Ncomps = NULL, CI_level = 95,
                          MCMC = TRUE, Nsamples = 10000, verbose = TRUE) {
  
  # if GROUPS & DV are NULL, then the groups variable is in 1st column & the DV is in the 2nd column of data
  
  # the groups variable can be numeric or categorical
  
  # var.equal -- a logical variable indicating whether to treat the two variances as being equal. 
  # If TRUE then the pooled variance is used to estimate the variance otherwise the  
  # Welch (or Satterthwaite) approximation to the degrees of freedom is used.
  
  # p.adjust.method options are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
  
  # Ncomps is the # of pairwise comparisons; if unspecified, it will be the # of all possible GROUPS comparisons
  
  
  if (is.null(GROUPS) & is.null(DV) )    donnesT <- data
  
  if (!is.null(GROUPS) & !is.null(DV) )  donnesT <- cbind(data[,GROUPS], data[,DV])
  
  grpnames <- unique(as.vector(as.matrix(donnesT[,1]))) 
  
  Ngroups <- length(grpnames)
  
  if (is.null(Ncomps))  Ncomps <- choose(Ngroups, 2)
  
  resultsM <- noms1 <- noms2 <- c()
  
  for (lupe1 in 1:(Ngroups-1)) {
    for (lupe2 in (lupe1+1):Ngroups) {
      
      dum <- subset(donnesT, donnesT[,1] == grpnames[lupe1] | donnesT[,1] == grpnames[lupe2]  )
      
      newdon <- data.frame(dum[,1], dum[,2])
      
      MNgrp1 <- mean(subset(newdon[,2],newdon[,1]==grpnames[lupe1]))
      MNgrp2 <- mean(subset(newdon[,2],newdon[,1]==grpnames[lupe2]))
      
      SDgrp1 <- stats::sd(subset(newdon[,2],newdon[,1]==grpnames[lupe1]))
      SDgrp2 <- stats::sd(subset(newdon[,2],newdon[,1]==grpnames[lupe2]))
      
      N1 <- nrow(subset(newdon,newdon[,1]==grpnames[lupe1]))
      N2 <- nrow(subset(newdon,newdon[,1]==grpnames[lupe2]))
      
      tresults <- stats::t.test(newdon[,2] ~ newdon[,1], data=newdon, 
                                var.equal=var.equal, conf.level = CI_level * .01) 
      tgroups  <- tresults$statistic
      dfgroups <- tresults$parameter
      plevel   <- tresults$p.value
      
      plevel.adj <- stats::p.adjust(plevel, method = p.adjust.method, n = Ncomps )
      
      # confidence intervals for the mean difference
      mdiff <- MNgrp1 - MNgrp2
      ci=tresults$conf.int
      ci.lb <- unlist(ci)[1]
      ci.ub <- unlist(ci)[2]
      
      # r effect size
      reffsize =  sqrt( tgroups**2 / (tgroups**2 + dfgroups) )
      
      # d effect size -- from R & R p 303 General Formula  best, because it covers = & not = Ns
      deffsize = (tgroups * (N1+N2)) / ( sqrt(dfgroups) * sqrt(N1*N2) )
      
      # g effect size -- from R compute.es.pdf
      df <- N1 + N2 - 2			
      bigJ <- 1 - (3 / (4 * df - 1))			  
      geffsize <- bigJ * deffsize   
      
      # BESD -- 2019 Funder, Ozer - Evaluating Effect Size in Psychological Research: Sense and Nonsense  p 160
      besd <- ((reffsize * 100) / 2) + 50
      
      # Bayes factor -- http://daniellakens.blogspot.com/2014/09/bayes-factors-and-p-values-for.html
      # Bayes factors from just t & Ns
      # BF_h0h1 <- exp(suppressMessages( ttest.tstat(tgroups, N1, N2, rscale=Bayes_rscale)$bf))

      # BF_h0h1_stat <- suppressMessages( ttest.tstat(tgroups, N1, N2, simple=TRUE))
      # BF_h1h0_stat <- 1 / BF_h0h1_stat
      
      BF_h1h0_stat <- suppressMessages( ttest.tstat(tgroups, N1, N2, simple=TRUE))
      BF_h0h1_stat <- 1 / BF_h1h0_stat
      
      
      if (!MCMC) {
        results <- cbind(lupe1, N1, MNgrp1, SDgrp1, lupe2, N2, MNgrp2, SDgrp2,
                         tgroups, dfgroups, plevel, plevel.adj, 
                         mdiff, ci.lb, ci.ub, deffsize, geffsize, reffsize, 
                         besd, BF_h1h0_stat, BF_h0h1_stat)
      }
      
      if (MCMC) {
        
        # MCMC / raw data Bayesian analyses
        x <- subset(newdon, newdon[,1] == grpnames[lupe1])[2]
        y <- subset(newdon, newdon[,1] == grpnames[lupe2])[2]
        Bayes_post <- ttestBF(x = unlist(x), y = unlist(y))
        chains <- posterior(Bayes_post, iterations = 10000, progress = FALSE)
        
        # parameters estimates & CIs -- in d value metric
        Bayes_d <- mean(chains[,2])
        quant_size <- (100 - CI_level) / 2
        quant_lb <- quant_size * .01
        quant_ub <- (100 - quant_size) * .01
        Bayes_d_ci_lb <- quantile(chains[,2], probs = quant_lb)
        Bayes_d_ci_ub <- quantile(chains[,2], probs = quant_ub)
        
        # Bayes factors from raw data analyses
        
        # BF_h0h1_raw <- extractBF(Bayes_post)[1]
        # BF_h1h0_raw <- 1 / BF_h0h1_raw
        
        BF_h1h0_raw <- extractBF(Bayes_post)[1]
        BF_h0h1_raw <- 1 / BF_h1h0_raw
        
        
        results <- cbind(lupe1, N1, MNgrp1, SDgrp1, lupe2, N2, MNgrp2, SDgrp2,
                         tgroups, dfgroups, plevel, plevel.adj, 
                         mdiff, ci.lb, ci.ub, deffsize, geffsize, reffsize, besd, 
                         BF_h1h0_raw, BF_h0h1_raw,
                         Bayes_d, Bayes_d_ci_lb, Bayes_d_ci_ub)
      }
      
      resultsM <- rbind(resultsM, results)
      
      noms1 <- append(noms1, grpnames[lupe1])
      noms2 <- append(noms2, grpnames[lupe2])		                  
    }  	
  }
  
  resultsM <- data.frame(resultsM)
  rownames(resultsM) <- c()
  
  colnames(resultsM)[1:21] <- 
    c("Group_1","N_1","Mean_1","SD_1","Group_2","N_2","Mean_2","SD_2",
      "t","df","p","p adj.",
      "Diff.","ci.lb","ci.ub","d","g","r","BESD", 
      "BF_h1h0", "BF_h0h1")
  
  if (MCMC)
    colnames(resultsM)[22:24] <- c("Bayes_d", "d_ci_lb", "d_ci_ub")
  
  resultsM[,1] <- noms1
  resultsM[,5] <- noms2
  
  if (verbose) {
    message("\n\nGroup differences - descriptive statistics:\n")
    print(round_boc(resultsM[c(1:8)], round_non_p = 2), print.gap=3)
    
    message("\n\nGroup differences - significance tests:\n")
    print(round_boc(resultsM[c(1,5,9:12)], round_non_p = 2), print.gap=3)
    
    if (!MCMC) {
      message("\n\nGroup differences - Bayes Factors:\n")
      print(round_boc(resultsM[c(1,5,20,21)], round_non_p = 2), print.gap=3)
    }
    
    if (MCMC) {
      message("\nGroup differences - Bayes Factors:\n")
      print(round_boc(resultsM[c(1,5,20,21)], round_non_p = 2), print.gap=3)
    }
    
    message("\nGroup differences - confidence intervals & effect sizes:\n")
    print(round_boc(resultsM[c(1,5,13:19)], round_non_p = 2), print.gap=3)
    
    if (MCMC) {
      message("\nGroup differences - Bayesian d values:\n")
      print(round_boc(resultsM[c(1,5,22:24)], round_non_p = 2), print.gap=3)
    }
  }
  
  return(invisible(resultsM))
}



