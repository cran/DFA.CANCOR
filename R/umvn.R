

umvn <- function(data) {
  
  # descriptive statistics & tests of univariate & multivariate normality -- from the MVN package
  
  # the MVN::mvn runs only 1 uni- and 1-multivariate test at a time, & the output format for Mardia is awful
  
  if (ncol(data) == 1) {
    
    datemp <- cbind(data, jitter(data)); colnames(datemp) <- c(colnames(data), 'jitd')
    
    res1 <- MVN::mvn(data=datemp, univariate_test = "SW")
    
    Shapiro_Wilk <- res1$"univariateNormality"[1,]
    
    descriptives <- res1$"Descriptives"
    
    Shapiro_Francia <- MVN::mvn(data=datemp, univariate_test = "SF")$"univariateNormality"[1,]
    
    Anderson_Darling <- MVN::mvn(data=datemp, univariate_test = "AD", desc = FALSE)$"univariateNormality"[1,]
    
    Cramer_von_Mises <- MVN::mvn(data=datemp, univariate_test = "CVM", desc = FALSE)$"univariateNormality"[1,]
    
    Lilliefors <- MVN::mvn(data=datemp, univariate_test = "Lillie", desc = FALSE)$"univariateNormality"[1,]
    
    univariate_tests <- rbind(Shapiro_Wilk, Shapiro_Francia, Anderson_Darling, Cramer_von_Mises, Lilliefors)
    
    multivariate_tests <- NA
  }
  
  if (ncol(data) > 1) {
    
    res1 <- MVN::mvn(data=data, mvn_test = "mardia", univariate_test = "SW")
    
    descriptives <- res1$"Descriptives"
    
    Shapiro_Wilk <- res1$"univariateNormality"
    
    Mardia       <- res1$"multivariateNormality"[1:2,]
    
    
    res2 <- MVN::mvn(data = data, mvn_test = "hz", univariate_test = "SF",  desc = FALSE)
    
    Shapiro_Francia <- res2$"univariateNormality"
    
    Henze_Zirkler   <- res2$"multivariateNormality"
    
    
    res3 <- MVN::mvn(data = data, mvn_test = "royston", univariate_test = "AD",  desc = FALSE)
    
    Anderson_Darling <-res3$"univariateNormality"
    
    Royston <- res3$"multivariateNormality"
    
    
    res4 <- MVN::mvn(data = data, mvn_test = "doornik_hansen", univariate_test = "CVM",  desc = FALSE)    
    
    Cramer_von_Mises <- res4$"univariateNormality"
    
    Doornik_Hansen <- res4$"multivariateNormality"
    
    
    Lilliefors <- MVN::mvn(data=data, univariate_test = "Lillie", desc = FALSE)$"univariateNormality"
    
    
    univariate_tests <- list(Shapiro_Wilk, Shapiro_Francia, Anderson_Darling, Cramer_von_Mises, Lilliefors)
    
    
    df1 <- data.frame(Test=rbind(Mardia[1], Henze_Zirkler[,1], Royston[,1], Doornik_Hansen[,1]))
    
    m2 <- as.numeric(levels(unlist(Mardia[2])))
    
    df2 <- data.frame(Statistic = rbind( m2[1], m2[2], Henze_Zirkler[,2], Royston[,2], Doornik_Hansen[,2]))
    
    m3 <- as.numeric(levels(unlist(Mardia[3])))
    
    df3 <- data.frame(p = rbind( m3[1], m3[2], Henze_Zirkler[,3], Royston[,3], Doornik_Hansen[,4]))
    
    m4 <- unlist(Mardia[4])
    
    df4 <- data.frame(Normality= rbind( m4[1], m4[2], Henze_Zirkler[,4], Royston[,4], Doornik_Hansen[,5]))
    
    multivariate_tests <- cbind(df1, df2, df3, df4); colnames(multivariate_tests)[4] <- 'Normality'
  }
  
  skewSE <- sqrt(6 / descriptives$n)
  skewZ  <- descriptives$Skew / skewSE	
  skewP  <- pnorm(abs(skewZ), lower.tail=FALSE) * 2     # 2-tailed sig test
  
  kurtosisSE <- sqrt(24 / descriptives$n)
  kurtosisZ  <- descriptives$Kurtosis / kurtosisSE	
  kurtosisP  <- pnorm(abs(kurtosisZ), lower.tail=FALSE) * 2     # 2-tailed sig test
  
  descriptives <- cbind(descriptives[,c(1:6,9)], skewZ, skewP, descriptives$Kurtosis, kurtosisZ, kurtosisP)
  colnames(descriptives)[1] <- 'N'
  colnames(descriptives)[8:12] <- c('z (Skew)','p (Skew)','Kurtosis','z (Kurtosis)','p (Kurtosis)')
  
  if (ncol(data) == 1) descriptives <- descriptives[1,]
  
  output <- list(  
    descriptives = descriptives,
    univariate_tests = univariate_tests,
    multivariate_tests = multivariate_tests
  )
  
  return(invisible(output))
  
}
