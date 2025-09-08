

DESCRIPTIVES <- function(data, groups=NULL, variables=NULL, CI_level=95, verbose=TRUE) {
  
  if (!is.null(groups) & is.null(variables)) 
    cat('\nA "groups" variable was provided but without any "variables" to be analyzed. Expect errors.\n')
  
  if (is.null(groups) & is.null(variables)) { 
    donnes <- as.data.frame(data)
    variables <- colnames(donnes)
  } else { donnes <- as.data.frame(data[,c(groups,variables)]) }
  
  if (!all( sapply(donnes[, variables], is.numeric))) {
    non_numeric_vars <- names(donnes)[!sapply(donnes, is.numeric)]
    cat('\nThe following non-numeric columns in data were removed from the analyses:\n\n')
    print(non_numeric_vars)
    nums <- sapply(donnes, is.numeric) 

    if (is.null(groups)) {
      donnes <- donnes[, nums]
      variables <- colnames(donnes)
    }
    
    if (!is.null(groups))  {
      donnes <- cbind(donnes[,groups], donnes[, nums])
      variables <- colnames(donnes)
      colnames(donnes)[1] <- groups
    }
  }
  
  
  # whole sample / all groups combined
  DESCR_whole <- descriptives_boc(data = donnes, variables=variables, CI_level=CI_level) 
  
  if (verbose) {		
    cat('\n\nDescriptive Statistics for the whole sample:\n')
    print((round(DESCR_whole,3)), print.gap=3)
  }
  
  outp <- list(DESCR_whole)
  
  
  # when there are groups
  if (!is.null(groups)) {
    
    grpnames <- as.vector(as.matrix(donnes[,1])) # group names, in the same order as in data
    grpnames <- unique(grpnames)
    grpnums  <- seq(1:length(grpnames))
    Ngroups  <- length(grpnames)
    
    # separate stats for each group
    DESCR_grps <- list()
    for (lupeg in 1:Ngroups) { 
      dum <- subset(donnes, donnes[,1] == grpnames[lupeg])
      
      grpnom <- paste(grpnames[lupeg])
      
      DESCR_grps[[grpnom]] <- descriptives_boc(data = dum[,2:ncol(dum)], 
                                               variables=colnames(dum)[2:ncol(dum)], CI_level=CI_level) 
      
      if (verbose) {		
        cat('\n\nDescriptive Statistics for group:', grpnames[lupeg], '\n')
        print((round(DESCR_grps[[grpnom]],3)), print.gap=3)
      }
    }
    
    # arranging output by variable
    DESCR_vars <- vector(mode = "list", length = length(variables))
    
    for (lupev in 1:length(variables)) {
      for (lupeg in 1:length(DESCR_grps)) {
        
        DESCR_vars[[lupev]] <- rbind(DESCR_vars[[lupev]], DESCR_grps[[lupeg]][lupev,])
      }  
      rownames(DESCR_vars[[lupev]]) <- grpnames
    }
    names(DESCR_vars) <- variables
    
    outp <- list(DESCR_whole=DESCR_whole, 
                 DESCR_grps=DESCR_grps, 
                 DESCR_vars=DESCR_vars)
  }
  
  #print(round(DESCoutput, 2))
  
  return(invisible(outp))
  
}




# descriptives function -- N, mean, SD, min, max, CIs

descriptives_boc <- function(data, variables, CI_level= 95) {
  
  # the critical z value (that corresponds to the specified CI) to be used in the CI computations
  zforCI <- qnorm((1 + CI_level * .01) / 2) 
  
  data <- as.data.frame(data[, variables, drop=FALSE])
  
  descrout <- cbind(
    apply(data, 2, function(x) length(which(!is.na(x)))),
    sapply(data, mean, na.rm=TRUE),
    sapply(data, sd,   na.rm=TRUE),
    sapply(data, min,  na.rm=TRUE),
    sapply(data, max,  na.rm=TRUE)
  )
  
  colnames(descrout) <- c('N', 'Mean', 'SD', 'min', 'max')
  rownames(descrout) <- variables
  
  descrout <- data.frame(descrout)
  
  descrout$CI_ub <- descrout$Mean + zforCI * (descrout$SD / sqrt(descrout$N))
  descrout$CI_lb <- descrout$Mean - zforCI * (descrout$SD / sqrt(descrout$N))
  
  # print(round( descrout, 2))
  
  return(invisible(descrout))
}



