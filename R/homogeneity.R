

homogeneity <- function(data, groups, variables, verbose = TRUE) {

donnes <- as.data.frame(data[,c(groups,variables)])

if (anyNA(donnes) == TRUE) {
	donnes <- na.omit(donnes)
	cat('\n\nCases with missing values were found and removed from the data matrix.\n\n')
}

grpnames <- as.vector(as.matrix(donnes[,1])) # group names, in the same order as in the data matrix
grpnames <- unique(grpnames)
ngroups  <- length(grpnames)
nDVs <- ncol(donnes) - 1
N <- nrow(donnes)

if (is.factor(donnes[,1]) == FALSE)  donnes[,1] <- factor( donnes[,1], ordered = FALSE, labels=grpnames)

grpFreqs <- as.matrix(table(donnes[,c(1)]))

logdetgrps <- 0 # Box's M test
BoxLogdets <- matrix(-9999,ngroups,1) # for Box's M test
homogeneityOutput <- list()
for (lupeg in 1:ngroups) {
	dum <- subset(donnes, donnes[,1] == grpnames[lupeg] )

	covmatgrp <- stats::cov(dum[,2:ncol(dum)])
	grpnom <- paste(grpnames[lupeg])		
	homogeneityOutput[[grpnom]]$covmatrix  <- covmatgrp

	if (verbose == TRUE) {
		cat('\nCovariance matrix for Group', paste(grpnames[lupeg]),'\n\n')
		print(round(covmatgrp,2))
	}

	logdetgrps <- logdetgrps + (nrow(dum) - 1) * log(det(covmatgrp)) # for Box's M test
	BoxLogdets[lupeg,1] <- log(det(covmatgrp)) # for Box's M test
}

# Bartlett's test of homogeneity of variances (parametric, for K samples)
Bartlett <- stats::bartlett.test(x=(donnes[,c(2:ncol(donnes))]), g=donnes[,1], data=donnes)
Bartlett <- c(Bartlett$statistic, Bartlett$parameter, Bartlett$p.value)

# Fligner-Killeen test of homogeneity of variances (non parametric, for K samples)
Fligner_Killeen <- c(stats::fligner.test(donnes[,c(2:ncol(donnes))], donnes[,1], data=donnes))
Fligner_Killeen <- c(Fligner_Killeen$statistic, Fligner_Killeen$parameter, Fligner_Killeen$p.value)
names(Fligner_Killeen)[1] <- 'Fligner_Killeen chi-squared'

# The Levene test is an alternative to the Bartlett test that is less 
# sensitive to departures from normality.
# Levene's test: Compare the variances of k samples, where k can be more than two samples. 
# The function leveneTest() [in car package] can be used.

# for the var-covar matrices from SPSS - requires Type III sums of squares
# using MANOVA to obtain the sums of squares and cross-products matrix for error, &
# the sums of squares and cross-products matrix for group/IV
# www.webpages.uidaho.edu/~kirk/calendar/R/MANOVA.doc
MV2 <-manova( as.matrix(donnes[,2:ncol(donnes)]) ~ donnes[,1], data=donnes)
sscpWithin <- (N-1)*cov(MV2$residuals) # E
sscpBetween <- (N-1)*cov(MV2$fitted.values) # H

# # using Anova from the car package to obtain the sums of squares and cross-products matrices
# outcome <- as.matrix(donnes[,c(2:ncol(donnes))])
# fit <- stats::lm(outcome ~ as.factor(donnes[,c(1)]), data = donnes)
# sscp <- car::Anova(fit, type="III") #summary(sscp)
# sscpWithin <- sscp$SSPE
# sscpBetween <- sscp$SSP$'as.factor(donnes[, c(1)])'

# # another approach
# options(contrasts = c('contr.sum','contr.poly'))
# Next, store the model:
# model <- lm(time ~ topic * sys, data=search)
# Finally, call the drop1 function on each model component:
# drop1(model, .~., test='F')
# The results give the type III SS, including the p-values from an F-test.


# the pooled within groups covariance matrix from SPSS:
PooledWithinCovarSPSS <- sscpWithin * (1/(N-ngroups))

# the pooled within groups correlation matrix from SPSS:
bigDv <- diag(suppressWarnings(sqrt(1 / PooledWithinCovarSPSS)))
scalemat <- diag(bigDv)
PooledWithinCorrelSPSS <- t(scalemat) %*% PooledWithinCovarSPSS %*% scalemat
rownames(PooledWithinCorrelSPSS) <- rownames(PooledWithinCovarSPSS)
colnames(PooledWithinCorrelSPSS) <- colnames(PooledWithinCovarSPSS)

# Box's M test
# formulas from http://www.real-statistics.com/multivariate-statistics/boxs-test-equality-covariance-matrices/boxs-test-basic-concepts/
# IBM SPSS Statistics Algorithms 20
logdetpw <- log(det(PooledWithinCovarSPSS)) 
BoxM <- ( (N - ngroups) * logdetpw ) - logdetgrps
df1 <- ( (ngroups-1) * nDVs * (nDVs+1) ) / 2
c <- (((2*(nDVs**2))+(3*nDVs)-1) / (6*(nDVs+1)*(ngroups-1))) *
     (sum( (1/(grpFreqs-1))) - (1 / (N-ngroups)))   
c2 <- (((nDVs-1)*(nDVs+2)) / (6*(ngroups-1)))  * 
      (sum( (1/((grpFreqs-1)^2)) )  - (1 / ((N-ngroups)^2))) 
df2 <- (df1 + 2) / abs(c2-c**2)
aplus <- df1 / ( 1 - c - (df1/df2))
Fplus <- BoxM / aplus
aminus <- df2 / (1 - c + (2/df2))
Fminus <- (df2 * BoxM) / (df1 * (aminus - BoxM))
if (c2 > c**2)  bigF <- Fplus
if (c2 < c**2)  bigF <- Fminus
pbigF <- stats::pf(bigF, df1, df2, lower.tail=FALSE) # p level
#rownames(BoxLogdets) <- paste("Group ", 1:ngroups, sep="") 
rownames(BoxLogdets) <- paste(grpnames) 
colnames(BoxLogdets) <- 'Log Determinant'
Pooled <- logdetpw
BoxLogdets <- rbind(BoxLogdets, Pooled)
BoxMtest <- c(BoxM, bigF, df1, df2, pbigF)


if (verbose == TRUE) {

#	cat('\n\nTests for homogeneity of variances & covariances:\n')	
	
	cat('\n\nBartlett test of homogeneity of variances (parametric):\n')
	cat('\nBartlett\'s K-squared =', round(Bartlett[1],3), '  df =', Bartlett[2],
	    '  p value =', round(Bartlett[3],5) )
	
	cat('\n\n\nFligner-Killeen test of homogeneity of variances (non parametric):\n')
	cat('\nFligner-Killeen chi-squared =', round(Fligner_Killeen[1],3), 
	    '  df =', Fligner_Killeen[2], '  p value =', round(Fligner_Killeen[3],5) )
	
	cat('\n\n\nPooled within groups covariance matrix from SPSS:\n')
	print(round(PooledWithinCovarSPSS,3))
	
	cat('\n\nPooled within groups correlation matrix from SPSS:\n')
	print(round(PooledWithinCorrelSPSS,3))
	
	cat('\n\nBox Test of equality of covariance matrices:\n')
	cat('\nLog determinants:\n')
	print(round(BoxLogdets,3))
	cat('\n\nM =', round(BoxMtest[1],3), '  F =', round(BoxMtest[2],3), '  df1 =', round(BoxMtest[3],2), 
	    '  df2 =', round(BoxMtest[4],2), '  p = ', round(BoxMtest[5],5), '\n\n\n')
}	

homogeneityOutput$Bartlett = Bartlett
homogeneityOutput$Fligner_Killeen = Fligner_Killeen
homogeneityOutput$PooledWithinCovarSPSS = PooledWithinCovarSPSS
homogeneityOutput$PooledWithinCorrelSPSS = PooledWithinCorrelSPSS
homogeneityOutput$sscpWithin = sscpWithin
homogeneityOutput$sscpBetween = sscpBetween
homogeneityOutput$BoxLogdets = BoxLogdets
homogeneityOutput$BoxMtest = BoxMtest	

return(invisible(homogeneityOutput))

}


