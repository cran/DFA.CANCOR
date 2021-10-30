

# rounds numeric columns in a matrix
# numeric columns named 'p' or 'plevel' are rounded to round_p places
# numeric columns not named 'p' are rounded to round_non_p places

round_boc <- function(donnes, round_non_p = 3, round_p = 5) {
	
	# identify the numeric columns
	#	numers <- apply(donnes, 2, is.numeric)  # does not work consistently 
	for (lupec in 1:ncol(donnes)) {

		if (is.numeric(donnes[,lupec]) == TRUE) 
		
			if (colnames(donnes)[lupec] == 'p' | colnames(donnes)[lupec] == 'plevel')  {
				donnes[,lupec] <- round(donnes[,lupec],round_p)
			} else {
				donnes[,lupec] <- round(donnes[,lupec],round_non_p)				
			}		
		# if (is.numeric(donnes[,lupec]) == FALSE) numers[lupec] = 'FALSE'		
		# if (colnames(donnes)[lupec] == 'p') numers[lupec] = 'FALSE'		
	}
	
	# # set the p column to FALSE
	# numers_not_p <- !names(numers) %in% "p"
	
#	donnes[,numers_not_p] = round(donnes[,numers_not_p],round_non_p) 
	
#	if (any(colnames(donnes) == 'p'))  donnes[,'p'] = round(donnes[,'p'],round_p) 

	return(invisible(donnes))
}





linquad <- function(iv, dv) {

# uses lm for the linear and quadratic association between 2 continous variables

	linquadOutput <- list()
	coefs <- matrix(-9999,2,4)

	ivsqd <- iv^2

	modlin <- lm(dv ~ iv)
	summodlin <- summary(modlin)
	coefs[1,] <- summodlin$coefficients[2,]

	modquad <- lm(dv ~ iv + ivsqd)
	summodquad <- summary(modquad)
	coefs[2,] <- summodquad$coefficients[3,]

	colnames(coefs) <- colnames(summodlin$coefficients)
	colnames(coefs)[1:2] <- c('b','SE')
	# rownames(coefs) <- c('linear','quadratic')
	rownames(coefs) <- c('idv (from dv ~ idv)','idv**2 (from dv ~ idv + idv**2)')
	beta <- betaboc(coefs[,1],iv,dv)
	coefs <- cbind(coefs[,1:2],beta,coefs[,3:4])

linquadOutput <- coefs

return(invisible(linquadOutput))

}




betaboc <- function (b,iv,dv) {
	ivsd <- sd(iv)
	dvsd <- sd(dv)
	beta <- b * ivsd / dvsd
	return(beta)
}



umvn <- function(data) {
	
	# descriptive statistics & tests of univariate & multivariate normality -- from the MVN package
		
	# # uvn1 <- uniNorm(data, type = "SW", desc = TRUE)  # Shapiro-Wilk test of univariate normality
	# # uvn2 <- uniNorm(data, type = "CVM", desc = FALSE)  # Cramer- von Mises test of univariate normality
	# # uvn3 <- uniNorm(data, type = "Lillie", desc = FALSE)  # Lilliefors (Kolmogorov-Smirnov) test of univariate normality
	# # uvn4 <- uniNorm(data, type = "SF", desc = FALSE)  # Shapiro-Francia test of univariate normality
	# # uvn5 <- uniNorm(data, type = "AD", desc = FALSE)  # Anderson-Darling test of univariate normality
			

	res1 <- MVN::mvn(data = data, mvnTest = "mardia", univariateTest = "SW") 
	res2 <- MVN::mvn(data = data, mvnTest = "hz")
	res3 <- MVN::mvn(data = data, mvnTest = "royston")
	res4 <- MVN::mvn(data = data, mvnTest = "dh",  desc = FALSE)    

	descriptives <- res1$"Descriptives"
	descriptives <- descriptives[,-c(7,8)]
	
	skewSE <- sqrt(6 / descriptives$n)
	skewZ  <- descriptives$Skew / skewSE	
	skewP  <- pnorm(abs(skewZ), lower.tail=FALSE) * 2     # 2-tailed sig test
	
	kurtosisSE <- sqrt(24 / descriptives$n)
	kurtosisZ  <- descriptives$Kurtosis / kurtosisSE	
	kurtosisP  <- pnorm(abs(kurtosisZ), lower.tail=FALSE) * 2     # 2-tailed sig test

	descriptives <- cbind(descriptives[,1:7], skewZ, skewP, descriptives$Kurtosis, kurtosisZ, kurtosisP)
	colnames(descriptives)[8:12] <- c('Skew z','Skew p','Kurtosis','Kurtosis z','Kurtosis p')

	Shapiro_Wilk   <- res1$"univariateNormality"
	
	Mardia         <- res1$"multivariateNormality"[1:2,]
	Henze_Zirkler  <- res2$"multivariateNormality"
	Royston        <- res3$"multivariateNormality"
	Doornik_Hansen <- res4$"multivariateNormality"
	
	
#	result$multivariateOutliers


umvnoutput <- list(  
   descriptives = descriptives,
   Shapiro_Wilk = Shapiro_Wilk,
   Mardia = Mardia,
   Henze_Zirkler = Henze_Zirkler,
   Royston = Royston,
   Doornik_Hansen = Doornik_Hansen
)

return(invisible(umvnoutput))

}



squareTable <- function(var1, var2, grpnames) {
    Original  <- factor(var1, levels = grpnames)
    Predicted <- factor(var2, levels = grpnames)
    table(Original, Predicted)
}



kappa.cohen <- function (var1, var2, grpnames) {

	# the data for this function (kapdon) are the category values for each of 2 columns,
	# but the analyses are performed on the contingency table	

	kapdonCT <- squareTable(var1, var2, grpnames); kapdonCT

	# based on Valiquette 1994 BRM, Table 1
	n <- sum(kapdonCT)  # Sum of Matrix elements
	kapdonP <- kapdonCT / n  # table of proportions
	po <- sum(diag(kapdonP))
	c <- rowSums(kapdonP)
	r <- colSums(kapdonP)
	pe <- r %*% c
	num <- po - pe
	den <- 1 - pe
	kappa <- num / den

	# SE and variance of kappa from Cohen 1968
	sek <- sqrt((po*(1-po))/(n*(1-pe)^2))
		
	if (n < 100) var <- pe/(n*(1-pe))  # kappa variance as reported by Cohen in his original work
	
	if (n >= 100) {
		s <- t(matrix(colSums(kapdonP))) # columns sum
		var <- (pe+pe^2-sum(diag(r%*%s) * (r+t(s)))) / (n*(1-pe)^2)
		# asymptotic kappa variance as reported by 
		# Fleiss, J. L., Lee, J. C. M., & Landis, J. R. (1979).
	    # The large sample variance of kappa in the case of different sets of raters. 
	    # Psychological Bulletin, 86, 974-977
	}

	zkappa <- kappa / sqrt(    abs(var)     )
	
	sig <- round(pnorm(abs(zkappa),lower.tail = FALSE),5) * 2 # 2-tailed test

#	print( c(kappa, sek, var, zkappa, sig))

	return(invisible(c(kappa, zkappa, sig)))

	# # based on kappa.m	
	# n <- sum(kapdon) # Sum of Matrix elements	
	# kapdonP <- kapdon/n # proportion		
	# r <- matrix(rowSums(kapdonP)) # rows sum
	# s <- t(matrix(colSums(kapdonP))) # columns sum	
	# Ex <- (r) %*% s # expected proportion for random agree
	# f <- diag(1,3)	
	# # pom <- apply(rbind(t(r),s),2,min)	
	# po <- sum(sum(kapdonP * f))  # sum(sum(x.*f))
	# pe <- sum(sum(Ex * f))
	# k <- (po-pe)/(1-pe)
	# # km <- (pom-pe)/(1-pe) # maximum possible kappa, given the observed marginal frequencies
	# # ratio <- k/km # observed as proportion of maximum possible	
	# # kappa standard error for confidence interval as reported by Cohen in his original work
	# sek <- sqrt((po*(1-po))/(n*(1-pe)^2))	
	# var <- pe/(n*(1-pe))  # kappa variance as reported by Cohen in his original work
	# # var <- (pe+pe^2-sum(diag(r*s).*(r+s')))/(n*(1-pe)^2)  # for N > 100	
	# zk  <-  k/sqrt(var)

}



# Fleiss's kappa

# source: https://en.wikipedia.org/wiki/Fleiss%27_kappa

# Fleiss's kappa is a generalisation of Scott's pi statistic, a
# statistical measure of inter-rater reliability. It is also related to
# Cohen's kappa statistic. Whereas Scott's pi and Cohen's kappa work for
# only two raters, Fleiss's kappa works for any number of raters giving
# categorical ratings (see nominal data), to a fixed number of items. It
# can be interpreted as expressing the extent to which the observed amount
# of agreement among raters exceeds what would be expected if all raters
# made their ratings completely randomly. Agreement can be thought of as
# follows, if a fixed number of people assign numerical ratings to a number
# of items then the kappa will give a measure for how consistent the
# ratings are. The scoring range is between 0 and 1. 

# Conger, A.J. (1980). Integration and generalisation of Kappas for multiple raters. Psychological Bulletin, 88, 322-328. 
# Fleiss, J.L. (1971). Measuring nominal scale agreement among many raters. Psychological Bul- letin, 76, 378-382. 
# Fleiss, J.L., Levin, B., & Paik, M.C. (2003). Statistical Methods for Rates and Proportions, 3rd Edition. New York: John Wiley & Sons. 

kappa.fleiss <- function(var1, var2) {

	kapdon <- cbind(var1, var2)
	
	# the data for this function (kapdon) are the category values for each column,
	# but the analyses are performed on a count matrix (not a contin table) = the fleissmat below	
	fleissmat <- matrix(0,nrow(kapdon),max(kapdon))
	for (luper in 1:nrow(kapdon)) {
		for (lupec in 1:ncol(kapdon)) {
			fleissmat[luper,kapdon[luper,lupec]] <- fleissmat[luper,kapdon[luper,lupec]] + 1				
		}
	}	
	n <- nrow(fleissmat) 
	m <- sum(fleissmat[1,]) 	
	a <- n * m		
	pj <- colSums(fleissmat) / a 
	b <- pj * (1-pj)
	c <- a*(m-1)
	d <- sum(b)
	kj <- 1-(colSums((fleissmat * (m-fleissmat))) / (c %*% b)) # the value of kappa for the j-th category
	# sekj <- sqrt(2/c) 
	# zkj <- kj / sekj
	# pkj <- round(pnorm(abs(zkj),lower.tail = FALSE),5) * 2  # 2-tailed test
	k <- sum(b*kj, na.rm=TRUE) / d  # Fleiss's (overall) kappa
	sek <- sqrt(2*(d^2-sum(b *(1-2 *pj))))/sum(b *sqrt(c)) 	
	#ci <- k+(c(-1,1) * (pnorm(zkj)) * sek) 	
	zk <- k / sek  # normalized kappa
	sig <- round(pnorm(abs(zk),lower.tail = FALSE),5) * 2  # 2-tailed test
	
#	print( c(k, sek, zk, sig))
	
	return(invisible(c(k, zk, sig)))
}



kappas <- function(var1, var2, grpnames) {
	
	kc <- kappa.cohen(var1, var2, grpnames)  
	kf <- kappa.fleiss(var1, var2)
	kappasOUT <- rbind(kc,kf)
	kappasOUT[,1:2] <- round(kappasOUT[,1:2],3)
	kappasOUT[,3] <- round(kappasOUT[,3],5)	
	rownames(kappasOUT) <- c( "Cohen's kappa", "Fleiss's kappa")
	colnames(kappasOUT) <- c( "    kappa", "        z", "         p" )
		
	# k2 <- irr::kappa2( grpdat )  # Unweighted Kappa for categorical data without a logical order
	# kf <- irr::kappam.fleiss( grpdat, exact = FALSE, detail = TRUE)
	# kl <- irr::kappam.light( grpdat )
	# kappasOUT <- matrix( -9999, 3, 4)
	# kappasOUT[1,] <- cbind( k2$subjects, k2$value, k2$statistic, k2$p.value)
	# kappasOUT[2,] <- cbind( kl$subjects, kl$value, kl$statistic, kl$p.value)
	# kappasOUT[3,] <- cbind( kf$subjects, kf$value, kf$statistic, kf$p.value)
	# rownames(kappasOUT) <- c( "Cohen's kappa", "Light's kappa", "Fleiss's kappa")
	# colnames(kappasOUT) <- c( "   N", "    kappa", "         z", "      p" )

	return (kappasOUT)
}






# Press' Q 

# When DFA is used to classify
# individuals in the second or holdout sample. The percentage of cases that are
# correctly classified reflects the degree to which the samples yield consistent
# information. The question, then is what proportion of cases should be correctly
# classified? This issue is more complex than many researchers acknowledge. To
# illustrate this complexity, suppose that 75% of individuals are Christian, 15%
# are Muslim, and 10% are Sikhs. Even without any information, you could thus
# correctly classify 75% of all individuals by simply designating them all as
# Christian. In other words, the percentage of correctly classified cases should
# exceed 75%.

# Nonetheless, a percentage of 76% does not indicate that classification is
# significantly better than chance. To establish this form of significance, you
# should invoke Press' Q statistic. 

# Compute the critical value, which equals the chi-square value at 1 degree of
# freedom. You should probably let alpha equal 0.05. When Q exceeds this critical
# value, classification can be regarded as significantly better than chance,
# thereby supporting cross-validation.

# The researcher can use Press's Q statistic to compare with the chi-square critical 
# value of 6.63 with 1 degree of freedom (p < .01). If Q exceeds this critical value, 
# the classification can be regarded as significantly better than chance. 

# Press'sQ = (N _ (n*K))^2 / (N * (K - 1))

# where
# N = total sample size
# n = number of observations correctly classified 
# K = number of groups_
# Given that Press's Q = 29.57 > 6.63, it can be concluded that the classification 
# results exceed the classification accuracy expected by chance at a statistically 
# significant level (p < .01). 


PressQ <- function(freqs) {
	
	N <- sum(colSums(freqs))
	n <- sum(diag(freqs))
	k <- ncol(freqs)
	PressQ <- (N - (n*k))^2 / (N * (k - 1))
	
	return(invisible(PressQ))
}






############################# T tests  ############################################################


"ttestboc" <-  function (donnesT, var.equal=FALSE) {
	
# reads raw data; the groups variable is in 1st column; the DV(s) are in the subsequent columns
# the groups variable can be categorical
# the function compares the lowest & highest values of the group variable

# var.equal -- a logical variable indicating whether to treat the two variances as being equal. 
# If TRUE then the pooled variance is used to estimate the variance otherwise the  
# Welch (or Satterthwaite) approximation to the degrees of freedom is used.

# p.adjust.methods options are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"

#cat("\n\nGroups t-test:\n")


grpnames <- as.vector(as.matrix(donnesT[,1])) # group names, in the same order as in the data matrix
grpnames <- unique(grpnames)
#grpnums <- seq(1:length(grpnames))

donnesT[,1] <- as.numeric(donnesT[,1])

resultsM <- matrix(-9999,1,13)
#ngroups <- max(donnesT[,1])
ngroups <- length(grpnames)
grpnoms <- cbind(-9999,-9999)

for (lupe1 in 1:(ngroups-1)) {
	for (lupe2 in (lupe1+1):ngroups) {

		dum <- subset(donnesT, donnesT[,1] == lupe1 | donnesT[,1] == lupe2  )
		
		newdon = data.frame(dum[,1], dum[,2])

		groupmin <- min(newdon[,1])
		groupmax <- max(newdon[,1])

		mgrp1 <- mean(subset(newdon[,2],newdon[,1]==groupmin))
		mgrp2 <- mean(subset(newdon[,2],newdon[,1]==groupmax))

		sdgrp1 <- stats::sd(subset(newdon[,2],newdon[,1]==groupmin))
		sdgrp2 <- stats::sd(subset(newdon[,2],newdon[,1]==groupmax))

		N1 <- nrow(subset(newdon,newdon[,1]==groupmin))
		N2 <- nrow(subset(newdon,newdon[,1]==groupmax))

		SE1 <- sdgrp1 / sqrt(N1)
		SE2 <- sdgrp2 / sqrt(N2)

		tresults <- stats::t.test(newdon[,2]~newdon[,1],data=newdon, var.equal=var.equal) 
		tgroups  <- tresults$statistic
		dfgroups <- tresults$parameter
		plevel   <- tresults$p.value

		# r effect size
		reffsiz =  sqrt( tgroups**2 / (tgroups**2 + dfgroups)  )

		# d effect size -- from R&R p 303 General Formula  best, because covers = & not = Ns
		deffsiz = (tgroups * (N1+N2)) / ( sqrt(dfgroups) * sqrt(N1*N2) )

		results <- cbind( groupmin, N1, mgrp1, sdgrp1, groupmax, N2, mgrp2, sdgrp2,
		                  tgroups, dfgroups, plevel, reffsiz, abs(deffsiz) )

		results <- as.matrix(cbind(results))
		resultsM <- rbind( resultsM, results)
		
		grpnoms <- rbind( grpnoms, cbind(grpnames[lupe1], grpnames[lupe2]))
	}  	
}

grpnoms <- grpnoms[-c(1),]

resultsM2 <- data.frame(resultsM[-1,,drop=FALSE])

# resultsM2[,1] <- grpnoms[,1]   # pre April 2021
# resultsM2[,5] <- grpnoms[,2]   # pre April 2021

resultsM2[,1] <- grpnoms[1]
resultsM2[,5] <- grpnoms[2]

rownames(resultsM2) <- c()
colnames(resultsM2) <- c("Group","N","Mean","SD","Group","N","Mean","SD",
                         "t","df","p","r effsize","d effsize")

return(invisible(as.data.frame(resultsM2)))
}







# # sources for the following multivariate significance tests:

# Manly, B. F. J., & Alberto, J. A. (2017). Multivariate Statistical 
# Methods: A Primer (4th Edition). Chapman & Hall/CRC, Boca Raton, FL.

# https://rpubs.com/aaronsc32/manova-test-statistics

# http://www.statpower.net/Software.html

# Gatignon, H. (2014). Statistical Analysis of Management Data. New York, NY: Springer.




Wilk <- function(rho, Ncases, p, q) {
	
	# rho = the canonical correlations (all of them)
	# Ncases = the sample size
	
	# for CCA, p = the number of variables in set 1
	# for CCA, q = the number of variables in set 2
	
	# for DFA, p = the number of discriminant functions
	# for DFA, q = the number of DVs (continuous variables)

	NCVs <- length(rho)
	eigval <- rho**2 / (1 - rho**2)			
	
	m <- Ncases - 3/2 - (p + q)/2
	
	wilks <- Fwilks <- df1 <- df2 <- vector("numeric", NCVs) 

	for (lupe in 1:NCVs) {
		
		wilks[lupe] <- prod(1 / (1 + eigval[lupe:NCVs]))	
		
		t <- sqrt((p^2 * q^2 - 4)/(p^2 + q^2 - 5))
		
		df1[lupe] <- p * q
		df2[lupe] <- m * t - p * q/2 + 1
		
		Fwilks[lupe] <- (1 - wilks[lupe]^(1/t)) / wilks[lupe]^(1/t) * df2[lupe] / df1[lupe]
		
		p <- p - 1
		q <- q - 1
	}
	
	pvalue <- pf(Fwilks, df1, df2, lower.tail = FALSE)
	
	WilksOutput <- cbind(wilks,Fwilks,df1,df2,pvalue)
	return(invisible(WilksOutput))
}






Pillai <- function(rho, Ncases, p, q) {
		
	# rho = the canonical correlations (all of them)
	# Ncases = the sample size
	
	# for CCA, p = the number of variables in set 1
	# for CCA, q = the number of variables in set 2
	
	# for DFA, p = the number of discriminant functions
	# for DFA, q = the number of DVs (continuous variables)

	pManly  <- max(p,q) # the # of variables in the larger set
	pManly2 <- p + q  # the total # of variables

	NCVs <- length(rho)
	eigval <- rho**2 / (1 - rho**2)			

	# 2017 Manly p 69 -- works 
	m = 1 # number of samples
	m2 = 1 # number of samples
	d = max(pManly,(m-1)) # the greater of pManly and m - 1
	
	Pillai <- Fpillai <- df1 <- df2 <- vector("numeric", NCVs)
	s2 <- NCVs
	
	for (lupe in 1:NCVs) {
		
		Pillai[lupe] <- sum(eigval[lupe:NCVs] / (1 + eigval[lupe:NCVs]))  
		# https://rpubs.com/aaronsc32/manova-test-statistics	
		
		s <- NCVs + 1 - lupe  # number of positive eigenvalues
		
		d = max(pManly,(m-1)) # update d
		df1[lupe] <- s * d
		pManly = pManly - 1
		
		df2[lupe] <- s2 * (Ncases - m2 - pManly2 + s2)
	
		Fpillai[lupe] <- ((Ncases - m2 - pManly2 + s2) * Pillai[lupe]) / ( d * (s - Pillai[lupe]))
		
		m2 = m2 - 2
	}
	
	pvalue <- pf(Fpillai, df1, df2, lower.tail = FALSE)
		
	PillaiOutput <- cbind(Pillai,Fpillai,df1,df2,pvalue)
	return(invisible(PillaiOutput))
}





Hotelling <- function(rho, Ncases, p, q) {
			
	# rho = the canonical correlations (all of them)
	# Ncases = the sample size
	
	# for CCA, p = the number of variables in set 1
	# for CCA, q = the number of variables in set 2
	
	# for DFA, p = the number of discriminant functions
	# for DFA, q = the number of DVs (continuous variables)

	pManly  <- abs(p - q) # p  = 2 # number of variables  8 - 6?
	pManly2 <- p + q  # the total # of variables

	NCVs <- length(rho)
	eigval <- rho**2 / (1 - rho**2)			

	# 2017 Manly p 69 
	
	Hotelling <- Fhotelling <- df1 <- df2 <- vector("numeric", NCVs)
	s2 <- NCVs
	m  <- 1 # number of samples
	m2 <- 1 # number of samples

	for (lupe in 1:NCVs) {
		
		Hotelling[lupe] <- sum(eigval[lupe:NCVs])  # https://rpubs.com/aaronsc32/manova-test-statistics
			
		s <- NCVs + 1 - lupe  # number of positive eigenvalues
				
		A <- ( abs(m - pManly - 1) - 1) / 2
		df1[lupe] <- s * (2 * A + s + 1)
				
		B <- (Ncases - m2 - pManly2 - 1) / 2
		df2[lupe] <- 2 * (s2 * B + 1)
		m2 = m2 - 2
		
		Fhotelling[lupe] <- df2[lupe] * (Hotelling[lupe] / (s2 * df1[lupe]))
		}

	pvalue <- pf(Fhotelling, df1, df2, lower.tail = FALSE)
	HotellingOutput <- cbind(Hotelling,Fhotelling,df1,df2,pvalue)	
	return(invisible(HotellingOutput))
}




RoyRoot <- function(rho, Ncases, p, q) {

	# rho = the first canonical correlation (or all of them)
	# Ncases = the sample size
	
	# for CCA, p = the number of variables in set 1
	# for CCA, q = the number of variables in set 2
	
	# for DFA, p = the number of discriminant functions
	# for DFA, q = the number of DVs (continuous variables)

	NDVs <- p + q

	eval <- rho[1]**2 / (1 - rho[1]**2)  # convert rho (canonical r) to eigenvalue
	# based on 2017 Manly p 69 
	m = 1
	d = min(p,q) # should be d <- max(p,(m-1)), but m (# IVs) is 1 for 1-way MANOVA
	df1 = d
	df2 = Ncases - m - d   # - 1
	Froy = df2 / df1 * eval
	pvalue <- pf(Froy, df1, df2, lower.tail = FALSE)

	# 2017 Manly, Alberto - Multivariate Statistical Methods  p 68:
	# It may be important to know that what some computer programs call Roys largest 
	# root statistic is 1/(11) rather than 1 itself. 	
	lambda <- rho[1]**2

	# RoyRootOutput <- cbind(rho[1]**2,Froy,df1,df2,pvalue)
	RoyRootOutput <- cbind(eval[1], lambda, Froy,df1,df2,pvalue) #, rho[1], rho[1]**2, eval[1])

	return(invisible(RoyRootOutput))
}







# Bartlett's V test (peel-down) of the significance of canonical correlations (for CCA, not DFA)

BartlettV <- function(rho, Ncases, p, q) {
	
	# rho = the canonical correlations (all of them)
	# Ncases = the sample size
	
	# for CCA, p = the number of variables in set 1
	# for CCA, q = the number of variables in set 2
		
	NCVs <- length(rho)
	eigval <- rho**2 / (1 - rho**2)			
		
	wilks <- X2 <- df <- pvalue <- vector("numeric", NCVs) 

	for (lupe in 1:NCVs) {
		
		wilks[lupe] <- prod(1 / (1 + eigval[lupe:NCVs]))	
		
		X2[lupe] <- -( (Ncases - 1) - (p + q + 1)/2) * log(wilks[lupe])
		          
		df[lupe]  <- (p - lupe +1) * (q - lupe + 1) 
		          
		pvalue[lupe] <- pchisq(X2[lupe], df[lupe], ncp=0, lower.tail = FALSE) 
	}
		
#	print(cbind(wilks,X2,df,pvalue))
	BartlettVOutput <- cbind(round(wilks,2),round(X2,2),df,pvalue)
	return(invisible(BartlettVOutput))
}





# Rao's test (peel-down) of the significance of canonical correlations (for CCA, not DFA)

Rao <- function(rho, Ncases, p, q) {
	
	# rho = the canonical correlations (all of them)
	# Ncases = the sample size
	
	# for CCA, p = the number of variables in set 1
	# for CCA, q = the number of variables in set 2
		
	NCVs <- length(rho)
	eigval <- rho**2 / (1 - rho**2)			
		
	wilks <- Frao <- df1 <- df2 <- pvalue <- vector("numeric", NCVs) 

	for (lupe in 1:NCVs) {		
		wilks[lupe] <- prod(1 / (1 + eigval[lupe:NCVs]))			
		pnew <- p - lupe + 1
		qnew <- q - lupe + 1
		t <- (Ncases - 1) - (pnew + qnew + 1) / 2       
		s <- ifelse((pnew^2 + qnew^2) <= 5, 1, sqrt((pnew^2 * qnew^2 -4) / (pnew^2 + qnew^2 -5)))              
		df1[lupe] <- pnew * qnew 		
		df2[lupe] <- (1 + t*s - pnew*qnew/2)		
		Frao[lupe] <- ((1 - wilks[lupe]^(1/s)) / wilks[lupe]^(1/s)) * df2[lupe] / df1[lupe] 
		pvalue[lupe] <- suppressWarnings(pf(Frao[lupe], df1[lupe], df2[lupe], ncp=0, lower.tail = FALSE))
    }
#	print(cbind(wilks,Frao,df1,df2,pvalue))
	RaoOutput <- cbind(round(wilks,2),round(Frao,2),df1,round(df2,2),round(pvalue,5))
	return(invisible(RaoOutput))
}
  
   
   
