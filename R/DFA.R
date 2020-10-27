
DFA <- function(data, groups, variables, plot=TRUE,
                predictive=TRUE, priorprob='SIZES', verbose=TRUE) {

donnes <- as.data.frame(data[,c(groups,variables)])


if (anyNA(donnes) == TRUE) {
	donnes <- na.omit(donnes)
#	cat('\n\nCases with missing values were found and removed from the data matrix.\n\n')
	NAflag = TRUE
} else {
	NAflag = FALSE
}

grpnames <- as.vector(as.matrix(donnes[groups])) # group names, in the same order as in the data matrix
grpnames <- unique(grpnames)
grpnums  <- seq(1:length(grpnames))
Ngroups  <- length(grpnames)
Ndvs <- ncol(donnes) - 1
Ncases <- nrow(donnes)


if (is.factor(donnes[,1]) == FALSE)  donnes[,1] <- factor(donnes[,1], ordered = FALSE, labels=grpnames)

donnes <- as.data.frame(donnes)

grpFreqs <- as.matrix(table(donnes[,1]))


# from my homovarcovar function:
# for the var-covar matrices from SPSS - requires Type III sums of squares
# using MANOVA to obtain the sums of squares and cross-products matrix for error, &
# the sums of squares and cross-products matrix for group/IV
# www.webpages.uidaho.edu/~kirk/calendar/R/MANOVA.doc
MV2 <-manova( as.matrix(donnes[,2:ncol(donnes)]) ~ donnes[,1], data=donnes)
sscpwith <- (Ncases-1)*cov(MV2$residuals) # E
sscpbetw <- (Ncases-1)*cov(MV2$fitted.values) # H




# LDA from MASS package


#  the lda function from MASS produces different raw, lda coefficients when different priors are used
#  but SPSS produces the same coefficients regardless of the priors that are used
#  to produce the SPSS results, use priors based on the group sizes, as in prior=(grpFreqs/sum(grpFreqs))

#  from the Details for the lda function in MASS:
#  Specifying the prior will affect the classification unless over-ridden in predict.lda. Unlike in 
#  most statistical packages, it will also affect the rotation of the linear discriminants within their space,
#  as a weighted between-groups covariance matrix is used. Thus the first few linear discriminants emphasize 
#  the differences between groups with the weights given by the prior, which may differ from their prevalence in the dataset. 


# SPSS options for priors are "All groups equal" or "Compute from group sizes"
if (is.null(priorprob)) priorprob = 'SIZES'
if (priorprob == 'EQUAL') priors = matrix((1/Ngroups), 1, Ngroups) 
if (priorprob == 'SIZES') priors = grpFreqs/sum(grpFreqs) 

ldaoutput <- MASS::lda(x = as.matrix(donnes[,2:ncol(donnes)]), grouping=donnes[,1], prior=priors)

lda.values <- stats::predict(ldaoutput, donnes[,2:ncol(donnes)]) # obtain scores on the DFs


# eigenvalues, canonical correlations, & one-way anovas on the DFs
dfc <- data.frame(donnes[,1], lda.values$x)
evals <- matrix(-9999, ncol(ldaoutput$scaling), 4)
evals[,1] <- 1:ncol(ldaoutput$scaling)
anovaDFoutput <- matrix(-9999, ncol(ldaoutput$scaling), 5)
ttestDFoutput <- lapply(1:ncol(ldaoutput$scaling), function(x) matrix(-9999, nrow=choose(Ngroups,2), ncol=13))
names(ttestDFoutput) = c(paste("Discriminant Function ", 1:ncol(ldaoutput$scaling), sep=""))
for (luper in 1:nrow(evals)) {
	dd <- data.frame(dfc[,1], dfc[luper+1])
	colnames(dd) <- c('grp','dv')
	fit <- stats::lm(dv ~ as.factor(grp), data = dd)
	betwss <- stats::anova(fit)["as.factor(grp)", "Sum Sq"]
	withss <- stats::anova(fit)["Residuals", "Sum Sq"]
 	anovaDFoutput[luper,1] <- as.numeric(summary(fit)["r.squared"])
	anovaDFoutput[luper,2] <- stats::anova(fit)["as.factor(grp)","F value"]
	anovaDFoutput[luper,3] <- stats::anova(fit)["as.factor(grp)","Df"]
	anovaDFoutput[luper,4] <- fit$df.residual
	anovaDFoutput[luper,5] <- stats::anova(fit)["as.factor(grp)","Pr(>F)"]

	ttestDFoutput[[luper]] <- ttestboc(dd, var.equal=FALSE)			

	evals[luper,2] <- betwss / withss # eigenvalue
	evals[luper,4] <- sqrt(betwss / (betwss+withss)) # canonical correlation
}
sv <- ldaoutput$svd;  svproprotions <- sv^2/sum(sv^2) # % of variance
evals[,3] <- svproprotions
# cat('\n\n\nEigenvalues & canonical correlations:\n\n')
dimnames(evals) <- list(rep("", dim(evals)[1]))
colnames(evals) <- c('Function','eigenvalue','proportion of variance','canonical r')
#print(round(evals,3), print.gap=4)


# cat('\n\n\nMultivariate peel-down significance tests:\n\n')  # using p.asym from the CCP package

# July, 2019: for CCA, the CCP::p.asym function uses the #s of vars in set1 and set2
# as the values for p & q; but to get correct results for DFA using this function,
# p = the # of DVs and q = the # of DFs -- at least that what seems to yield correct values

rho <- evals[,4]
NCVs <- nrow(evals)

mv_Wilk <- Wilk(rho=rho, Ncases=Ncases, p = NCVs, q = Ndvs)
colnames(mv_Wilk) <- c('Wilk\'s Lambda', 'F-approx. ', 'df1', 'df2', 'p')
rownames(mv_Wilk) <- paste(1:nrow(mv_Wilk), paste("through ", nrow(mv_Wilk), sep = ""))
# print(round(mv_Wilk), print.gap=4); cat('\n\n')

mv_Pillai <- Pillai(rho=rho, Ncases=Ncases, p = NCVs, q = Ndvs)
colnames(mv_Pillai) <- c('Pillai-Bartlett Trace', 'F-approx. ', 'df1', 'df2', 'p')
rownames(mv_Pillai) <- paste(1:nrow(mv_Pillai), paste("through ", nrow(mv_Pillai), sep = ""))
# print(round(mv_Pillai), print.gap=4); cat('\n\n')

mv_Hotelling <- Hotelling(rho=rho, Ncases=Ncases, p = NCVs, q = Ndvs)
colnames(mv_Hotelling) <- c('Hotelling-Lawley Trace', 'F-approx. ', 'df1', 'df2', 'p')
rownames(mv_Hotelling) <- paste(1:nrow(mv_Hotelling), paste("through ", nrow(mv_Hotelling), sep = ""))
# print(round(mv_Hotelling), print.gap=4); cat('\n\n')

mv_Roy <- RoyRoot(rho=rho, Ncases=Ncases, p = NCVs, q = Ndvs)
colnames(mv_Roy) <- c('Roy\'s Largest Root', 'F-approx.', 'df1', 'df2', 'p')
rownames(mv_Roy) <- paste(1:nrow(mv_Roy), paste("through ", nrow(mv_Roy), sep = ""))
# print(round(mv_Roy), print.gap=4); cat('\n\n')



# cat('\n\n\n\nCanonical Discriminant Function (raw) Coefficients:\n')
colnames(ldaoutput$scaling) <-  c(paste("Function ", 1:ncol(ldaoutput$scaling), sep=""))
# print(round(ldaoutput$scaling), print.gap=4)

# centering each variable within groups
group.center <- function(var,grp) { return(var-tapply(var,grp,mean,na.rm=TRUE)[grp]) }
cdonnes <- matrix(-9999,nrow(donnes),(ncol(donnes)-1))
dfc <- cbind(donnes[,1], lda.values$x)
cdfc <- matrix(-9999,nrow(dfc),(ncol(dfc)-1))
for (lupec in 1:(ncol(donnes)-1)) { cdonnes[,lupec] <- group.center(donnes[,(lupec+1)], donnes[,1]) }

for (lupec in 1:(ncol(dfc)-1)) { cdfc[,lupec] <- group.center(dfc[,(lupec+1)], dfc[,1]) }
cdonnes <- cbind(donnes[,1], cdonnes) # placing the grouping variable back in the centered matrix
cdonnesdf <- data.frame(cdonnes)
structCoef <- stats::cor(x = cdonnes[,2:ncol(cdonnes)], y = cdfc)  # round(structCoef,2)
rownames(structCoef) <- rownames(ldaoutput$scaling)
colnames(structCoef) <- colnames(ldaoutput$scaling)
# cat('\n\nStructure Coefficients:\n')
colnames(structCoef) <-  c(paste("Function ", 1:ncol(structCoef), sep=""))
# print(round(structCoef), print.gap=4)

# standardized coefficients
pooledSDs <- as.matrix(apply(cdonnes, 2, FUN = sd)) # cdonnes contains the mean-centered data
standCoef <- (pooledSDs[2:nrow(pooledSDs),]) * ldaoutput$scaling
# cat('\n\nStandardized Coefficients:\n')
colnames(standCoef) <-  c(paste("Function ", 1:ncol(standCoef), sep=""))
# print(round(standCoef), print.gap=4)



# sscpwith <- sscps$sscpwith
# sscpbetw <- sscps$sscpbetw

# provides the standardized coefficients from SPSS:
poolwith <- sscpwith * (1/(nrow(donnes)-Ngroups))
pooledSDs <- sqrt(diag(poolwith)) # pooled SDs for SPSS results
standCoefSPSS <- pooledSDs * ldaoutput$scaling
# cat('\n\n\nStandardized Coefficients from SPSS:\n')
colnames(standCoefSPSS) <- c(paste("Function ", 1:ncol(standCoefSPSS), sep=""))
# print(round(standCoefSPSS), print.gap=4)

	
# group means & SDs on the raw ldfs
ldfscores  <- data.frame(donnes[,1],lda.values$x)
ldfscoresZ <- data.frame(donnes[,1],scale(lda.values$x))
centroids <- centroidSDs <- centroidsZ <- centroidSDsZ <- matrix(-9999,Ngroups,(ncol(ldfscores)-1))
for (lupec in 2:ncol(ldfscores)) {
	aggM  <- stats::aggregate(x = ldfscores[,lupec], by= list(ldfscores[,1]), FUN = mean)
	aggSD <- stats::aggregate(x = ldfscores[,lupec], by= list(ldfscores[,1]), FUN = sd)
	centroids[,(lupec-1)] <- aggM[,2]
	centroidSDs[,(lupec-1)] <- aggSD[,2]

	aggMZ  <- stats::aggregate(x = ldfscoresZ[,lupec], by= list(ldfscoresZ[,1]), FUN = mean)
	aggSDZ <- stats::aggregate(x = ldfscoresZ[,lupec], by= list(ldfscoresZ[,1]), FUN = sd)
	centroidsZ[,(lupec-1)] <- aggMZ[,2]
	centroidSDsZ[,(lupec-1)] <- aggSDZ[,2]
}
# cat('\n\n\nFunctions at Group Centroids\n')
# cat('\nUnstandardized canonical discriminant functions evaluated at group means:\n\n')
rownames(centroids) <- c(paste("Group ", aggM[,1], sep="")) 
colnames(centroids) <- c(paste("Function ", 1:(ncol(centroids)), sep=""))
# print(round(centroids), print.gap=4)

# cat('\n\nGroup Standard Deviations on the unstandardized functions:\n\n')
rownames(centroidSDs) <- c(paste("Group ", aggM[,1], sep="")) 
colnames(centroidSDs) <- c(paste("Function ", 1:(ncol(centroidSDs)), sep=""))
# print(round(centroidSDs), print.gap=4)

# cat('\n\nStandardized canonical discriminant functions evaluated at group means:\n\n')
rownames(centroidsZ) <- c(paste("Group ", aggM[,1], sep="")) 
colnames(centroidsZ) <- c(paste("Function ", 1:(ncol(centroidsZ)), sep=""))
# print(round(centroidsZ), print.gap=4)

# cat('\n\nGroup Standard Deviations on the standardized functions:\n\n')
rownames(centroidSDsZ) <- c(paste("Group ", aggM[,1], sep="")) 
colnames(centroidSDsZ) <- c(paste("Function ", 1:(ncol(centroidSDsZ)), sep=""))
# print(round(centroidSDsZ), print.gap=4)


# cat('\n\n\nOne-way ANOVAs using the scores on a discriminant function as the DV:\n\n')
dimnames(anovaDFoutput) <-list(rep("", dim(anovaDFoutput)[1]))
colnames(anovaDFoutput) <- c('Eta-squared','F','df','df','p')
rownames(anovaDFoutput) <- colnames(ldaoutput$scaling)
anovaDFoutput[,1:4] <- anovaDFoutput[,1:4]
anovaDFoutput[,5] <- anovaDFoutput[,5]
# print(anovaDFoutput, print.gap=4)


# one-way anovas & t-tests on the DVs
anovaDVoutput <- matrix(-9999, length(variables), 5)
ttestDVoutput <- lapply(1:length(variables), function(x) matrix(-9999, nrow=choose(Ngroups,2), ncol=13))
names(ttestDVoutput)=variables 
for (lupec in 1:length(variables)) {
	ddd <- data.frame(donnes[,1], donnes[,(lupec+1)])
	colnames(ddd) <- c('grp','dv')
	fit <- stats::lm(dv ~ as.factor(grp), data = ddd)
	betwss <- stats::anova(fit)["as.factor(grp)", "Sum Sq"]
	withss <- stats::anova(fit)["Residuals", "Sum Sq"]
 	anovaDVoutput[lupec,1] <- as.numeric(summary(fit)["r.squared"])
	anovaDVoutput[lupec,2] <- stats::anova(fit)["as.factor(grp)","F value"]
	anovaDVoutput[lupec,3] <- stats::anova(fit)["as.factor(grp)","Df"]
	anovaDVoutput[lupec,4] <- fit$df.residual
	anovaDVoutput[lupec,5] <- stats::anova(fit)["as.factor(grp)","Pr(>F)"]
	
	ttestDVoutput[[lupec]] <- ttestboc(ddd, var.equal=FALSE)			
}
# cat('\n\n\n\nOne-way ANOVAs on the original DVs\n') 
# cat('\n(provided for comparisons with the ANOVAs on the discriminant functions):\n\n')
dimnames(anovaDVoutput) <-list(rep("", dim(anovaDVoutput)[1]))
colnames(anovaDVoutput) <- c('Eta-squared','F','df','df','p')
rownames(anovaDVoutput) <- variables
anovaDVoutput[,1:4] <- anovaDVoutput[,1:4]
anovaDVoutput[,5] <- anovaDVoutput[,5]
# print(anovaDVoutput, print.gap=4)


DFAoutput <- list(  
   rawCoef=ldaoutput$scaling,
   structCoef=structCoef,
   standCoef=standCoef,
   standCoefSPSS=standCoefSPSS,
   centroids=centroids,
   centroidSDs=centroidSDs,
   centroidsZ=centroidsZ,
   centroidSDsZ=centroidSDsZ,
   DFAscores=lda.values$x,
   anovaDFoutput = anovaDFoutput,
   anovaDVoutput = anovaDVoutput
)



# 2014 Bird - Controlling the Maximum Familywise Type I Error Rate in Analyses of Multivariate Experiments

lsnoms <- names(ttestDVoutput)

if (Ngroups == 2) {
	# Two Groups 
	# Critical values for analyses controlling MFWER when k = 2 
	Ncases <- 40
	P <- Ndvs
	ALPHA <- 0.05 
	# Critical t and p values and confidence level for t* analysis 
	TCRIT <- qt(1-(1-(1-ALPHA)^(1/P))/2,Ncases-2)
	PCRIT <- 1-(1-ALPHA)^(1/P)
	CLEVEL <- 1-PCRIT 
	# Critical t and p values for Mt** analysis 
	TCRIT2 <- qt(1-(1-(1-ALPHA)**(1/(P-1)))/2,Ncases-2) 
	PCRIT2 <- 1-(1-ALPHA)^(1/(P-1))

	for (lupe in 1:length(ttestDVoutput)) {
		MFWER1 <- ttestDVoutput[[lupe]]
		MFWER1 <- MFWER1[,c(1,5,9)]
		MFWER1[,4] <- TCRIT
		MFWER1[,5] <- ifelse( abs(MFWER1[,3]) > TCRIT, 'significant', 'not significant')
		MFWER1[,6] <- CLEVEL
		colnames(MFWER1)[4:6] <- c('t1-critical','t1 decision','Confidence Level')		
		DFAoutput[[variables[lupe]]]$MFWER1.sigtest <- MFWER1		

		MFWER2 <- ttestDVoutput[[lupe]]
		MFWER2 <- MFWER2[,c(1,5,9)]
		MFWER2[,4] <- TCRIT2
		MFWER2[,5] <- ifelse( abs(MFWER2[,3]) > TCRIT2, 'significant', 'not significant')
		colnames(MFWER2)[4:5] <- c('t2-critical','t2 decision')
		DFAoutput[[variables[lupe]]]$MFWER2.sigtest <- MFWER2	
	}
}


if (Ngroups > 2) {

	# More Than Two Groups  
	# Critical values for analyses controlling MFWER when k >= 3 
	# Ncases <- Ncases
	K <- Ngroups
	P <- Ndvs 
	ALPHA <- .05
	# Critical t and p values and confidence level for t' analysis
	TCRIT <- qtukey((1-ALPHA)**(1/P),K,Ncases-K)/sqrt(2)
	PCRIT <- 2*(1-pt(TCRIT,Ncases-K))
	CLEVEL <- 1-PCRIT
	# Critical F and critical t and p values for F't'' analysis
	FCRIT <- qf((1-ALPHA)^(1/P),K-1,Ncases-K)
	TCRIT2 <- qtukey((1-ALPHA)^(1/P),K-1,Ncases-K)/sqrt(2)
	PCRIT2 <- 2*(1-pt(TCRIT2,Ncases-K))
	decisionF <- c(1:length(ttestDVoutput))
	
	for (lupe in 1:length(ttestDVoutput)) {
		MFWER1 <- ttestDVoutput[[lupe]]
		MFWER1 <- MFWER1[,c(1,5,9)]
		MFWER1[,4] <- TCRIT
		MFWER1[,5] <- ifelse( abs(MFWER1[,3]) > TCRIT, 'significant', 'not significant')
		MFWER1[,6] <- CLEVEL
		colnames(MFWER1)[4:6] <- c('t1-critical','t1 decision','Confidence Level')		
		DFAoutput[[variables[lupe]]]$MFWER1.sigtest <- MFWER1		

		decisionF[lupe] <- ifelse( anovaDVoutput[lupe,2] > FCRIT, 'significant', 'not significant')
		MFWER2 <- ttestDVoutput[[lupe]]
		MFWER2 <- MFWER2[,c(1,5,9)]
		MFWER2[,4] <- TCRIT2
		MFWER2[,5] <- ifelse( abs(MFWER2[,3]) > TCRIT2, 'significant', 'not significant')
		colnames(MFWER2)[4:5] <- c('t2-critical','t2 decision')
		DFAoutput[[variables[lupe]]]$MFWER2.sigtest <- MFWER2		
	}
}
	
	

if (predictive == TRUE | is.null(predictive)) {

	freqs_OR <- squareTable(var1=donnes[,1], var2=lda.values$class, grpnames)
	
	PropOrigCorrect <- round((sum(diag(freqs_OR)) / sum(freqs_OR)),3)
	
	chi_square_OR <- summary(freqs_OR)
	PressQ_OR <- PressQ(freqs_OR)
	
	rowfreqs_OR <- margin.table(freqs_OR, 1)
	colfreqs_OR <- margin.table(freqs_OR, 2)
	
	cellprops_OR <- prop.table(freqs_OR)
	rowprops_OR  <- prop.table(freqs_OR, 1)
	colprops_OR  <- prop.table(freqs_OR, 2)
	
	# kappas_cvo <- kappas(stats::na.omit(cbind(lda.values$class, donnes[,1])))
	# kappas_cvo <- kappas(freqs)
	  kappas_cvo_OR <- kappas(var1=donnes[,1], var2=lda.values$class, grpnames)

	
	# Frequencies: Original vs Cross-Validated (leave-one-out cross-validation)
	
	# classifications from leave-one-out cross-validation
	ldaoutput_CV <- MASS::lda(x = as.matrix(donnes[,c(2:ncol(donnes))]),
	   grouping=donnes[,1], prior=priors, CV = TRUE)
	
	# freqs_cvp <- data.frame(cbind(ldaoutput_CV$class, lda.values$class))
	# colnames(freqs_cvp) <-  c("Cross-Validated", "Predicted") 
	freqs_CV <-  squareTable(ldaoutput_CV$class, lda.values$class, grpnames)
	colnames(freqs_CV) <- paste(grpnames)
	rownames(freqs_CV) <- paste(grpnames)

	PropCrossValCorrect <- sum(diag(freqs_CV)) / sum(freqs_CV)
	
	chi_square_CV <- summary(freqs_CV)
	PressQ_CV <- PressQ(freqs_CV)

	rowfreqs_CV <- margin.table(freqs_CV, 1)
	colfreqs_CV <- margin.table(freqs_CV, 2)
	
	cellprops_CV <- prop.table(freqs_CV)
	rowprops_CV  <- prop.table(freqs_CV, 1)
	colprops_CV  <- prop.table(freqs_CV, 2)
	
	# Agreement (kappas) between the Cross-Validated and Original Group Memberships
	# kappas_CVO <- kappas(stats::na.omit(cbind(ldaoutput_CV$class, donnes[,1])))
	# freqs_CV <-  squareTable(ldaoutput_CV$class, donnes[,1])
	# colnames(freqs_CV) <- paste(grpnames)
	# rownames(freqs_CV) <- paste(grpnames)
	# kappas_CVO <- kappas(freqs_CV)
	kappas_CVO <- kappas(ldaoutput_CV$class, donnes[,1], grpnames)

	# Agreement (kappas) between the Cross-Validated and Predicted Group Memberships
	# kappas_CVP <- kappas(stats::na.omit(cbind(ldaoutput_CV$class, lda.values$class)))
	# kappas_CVP <- kappas(freqs_CV)
	kappas_CVP <- kappas(ldaoutput_CV$class, lda.values$class, grpnames)
}



# plot
if (plot == TRUE) {

	colnames(centroidsZ) <-  c(paste("DF ", 1:ncol(centroidSDsZ), sep=""))
	
	graphics::matplot(1:length(grpnames), centroidsZ, type = "l", lty=1, lwd=3, 
	        xaxt='n', xlab=groups, cex.axis=1.2, cex.lab = 1.3,
	        ylab='Discriminant Function z Scores', ylim = c(-3, 3), cex.axis=1.2        )
	graphics::axis(side=1, at=grpnums, labels=rownames(centroidsZ), xlab="groups")
	graphics::title(main='Mean Standardized Discriminant Function Scores for the Groups')
	graphics::legend("topright", legend = colnames(centroidsZ), bty="n", lwd=2, col=1:ncol(centroidsZ))
}



if (verbose == TRUE) {
	
	cat('\n\n\nLinear Discriminant Function Analysis:\n')
	
	# eigenvalues, canonical correlations, & one-way anovas on the DFs
	cat('\n\n\nEigenvalues & canonical correlations:\n\n')
	print(round_boc(evals), print.gap=4)
	
	cat('\n\n\nMultivariate peel-down significance tests:\n\n')  # using p.asym from the CCP package
	
	print(round_boc(mv_Wilk), print.gap=4); cat('\n\n')

	print(round_boc(mv_Pillai), print.gap=4); cat('\n\n')
	
	print(round_boc(mv_Hotelling), print.gap=4); cat('\n\n')
	
	print(round_boc(mv_Roy), print.gap=4); cat('\n\n')
	
	cat('\n\n\nCanonical Discriminant Function (raw) Coefficients:\n')
	print(round_boc(ldaoutput$scaling), print.gap=4)
	
	cat('\n\nStructure Coefficients:\n')
	print(round_boc(structCoef), print.gap=4)
	
	cat('\n\nStandardized Coefficients:\n')
	print(round_boc(standCoef), print.gap=4)
	
	cat('\n\n\nStandardized Coefficients from SPSS:\n')
	print(round_boc(standCoefSPSS), print.gap=4)
	
	cat('\n\n\nFunctions at Group Centroids\n')
	cat('\nUnstandardized canonical discriminant functions evaluated at group means:\n\n')
	print(round_boc(centroids), print.gap=4)
	
	cat('\n\nGroup Standard Deviations on the unstandardized functions:\n\n')
	print(round_boc(centroidSDs), print.gap=4)
	
	cat('\n\nStandardized canonical discriminant functions evaluated at group means:\n\n')
	print(round_boc(centroidsZ), print.gap=4)
	
	cat('\n\nGroup Standard Deviations on the standardized functions:\n\n')
	print(round_boc(centroidSDsZ), print.gap=4)
	
	cat('\n\n\nOne-way ANOVAs using the scores on a discriminant function as the DV:\n\n')
	print(round_boc(anovaDFoutput), print.gap=4)
	
	cat('\n\n\n\nOne-way ANOVAs on the original DVs\n') 
	cat('\n(provided for comparisons with the ANOVAs on the discriminant functions):\n\n')
	print(round_boc(anovaDVoutput), print.gap=4)
	
	cat('\n\n\nt-tests and effect sizes for group differences on the discriminant functions:\n\n\n')
	DFnames <- names(ttestDFoutput)
	for (lupe in 1:length(ttestDFoutput)) {
		cat(DFnames[lupe],'\n\n')
		print(round_boc(ttestDFoutput[[lupe]]), row.names = FALSE, print.gap=3); cat('\n\n')
	}
	
	cat('\n\n\n\nt-tests and effect sizes for group differences on the original DVs\n')
	cat('\n(provided for comparisons with the t-tests on the discriminant functions):\n\n\n')
	for (lupe in 1:length(ttestDVoutput)) {
		cat(lsnoms[lupe],'\n\n')
		print(round_boc(ttestDVoutput[[lupe]]), row.names = FALSE, print.gap=3); cat('\n\n')
	}
	

			
	# 2014 Bird - Controlling the Maximum Familywise Type I Error Rate in Analyses of Multivariate Experiments
	
	cat('\n\nAppropriate univariate significance tests (ignoring the discriminant functions)')
	cat('\nfor controlling the maximum familywise Type I error rate (MFWER), based on:')
	cat('\n\nBird, K. D., & Hadzi-Pavlovic, D. (2014). Controlling the maximum familywise Type I ')
	cat('\nerror rate in analyses of multivariate experiments. Psychological Methods, 19(2), 265-280.\n')
	
	if (Ngroups == 2) {
		
		cat('\nWhen there are only two groups and the researcher decides to control the MFWER') 
		cat('\nat .05 by (only) carrying out multiple t tests, then the MANOVA test is')
		cat('\nirrelevant and the critical t value is the t1-critical value below.')
		cat('\nThe t value for the data must be greater than t1-critical value for a pairwise')
		cat('\ncomparison to be significant.\n\n\n')
						
		for (lupe in 1:length(ttestDVoutput)) {
			cat(lsnoms[lupe],'\n\n')
			print(round_boc(DFAoutput[[variables[lupe]]]$MFWER1.sigtest), row.names = FALSE, print.gap=4); cat('\n\n')
		}		
	
		cat('\n\nWhen there are only two groups and the researcher prefers the two-stage')
		cat('\napproach to controling the MFWER at .05, then the MANOVA F test must be')
		cat('\nsignificant and the t value for each comparison must be greater than') 
		cat('\nt2-critical value below for a pairwise comparison to be significant.\n\n\n')
	
		for (lupe in 1:length(ttestDVoutput)) {
			cat(lsnoms[lupe],'\n\n')
			print(round_boc(DFAoutput[[variables[lupe]]]$MFWER2.sigtest), row.names = FALSE, print.gap=4); cat('\n\n')
		}		
	}
		
	if (Ngroups > 2) {
			
		cat('\n\nFor designs with more than two groups, Bird et al. (2014) could find no')
		cat('\njustification for the inclusion of an initial MANOVA test in MCPs designed to') 
		cat('\ncontrol the MFWER with protected t tests.\n')
		
		cat('\nWhen a researcher decides to control the MFWER at .05 by solely carrying out') 
		cat('\nmultiple t tests, then the MANOVA and ANOVA tests are irrelevant and the')
		cat('\ncritical t value is the t1-critical value below. The t value for the data must')
		cat('\nbe greater than t1-critical value for a pairwise comparison to be significant.\n\n\n')
			
		for (lupe in 1:length(ttestDVoutput)) {
			cat(lsnoms[lupe],'\n\n')
			print(round_boc(DFAoutput[[variables[lupe]]]$MFWER1.sigtest), row.names = FALSE, print.gap=4); cat('\n\n')
		}		
	
		cat('\n\nWhen a researcher prefers the two-stage approach to controlling the MFWER at .05,')
		cat('\nthen the ANOVA F value for the data must be greater than the F-critical value below,')
		cat('\nand the t value for each comparison must be greater than t2-critical value') 
		cat('\nbelow for a pairwise comparison to be significant.\n\n\n')
	
		for (lupe in 1:length(ttestDVoutput)) {
			cat(lsnoms[lupe],'\n\n')
			cat('\nF =', anovaDVoutput[lupe,2],'F-critical =', FCRIT,'decision =', decisionF[lupe],'\n\n')
			print(round_boc(DFAoutput[[variables[lupe]]]$MFWER2.sigtest), row.names = FALSE, print.gap=4); cat('\n\n')
		}		
	}

		
	if (predictive == TRUE | is.null(predictive)) {
	
		cat('\n\n\nPREDICTIVE DISCRIMINANT ANALYSIS\n')
		
		cat('\n\nPrior Probabilities for Groups:\n'); print(round(ldaoutput$prior,3), print.gap=4)
		
		cat('\n\nCross-Tabulation of the Original and Predicted Group Memberships:\n\n'); print(freqs_OR, print.gap=4)	
		
		cat('\n\nProportion of original grouped cases correctly classified:  ', round(PropOrigCorrect,3))
		
		cat('\n\n\nChi-square test of independence:\n'); print(chi_square_OR, print.gap=4) # chi-square test of indepedence

		cat('\n\nPress\'s Q significance test of classifiation accuracy:\n')
		
		if (PressQ_OR < 3.8415) { 
			cat('\nPress\'s Q =', PressQ_OR, ', which is < 3.8415, indicating non-significance')
		} else if (PressQ_OR > 6.63501) { 
			cat('\nPress\'s Q =', PressQ_OR, 'which is > 6.63501, indicating p < .01')
		} else if (PressQ_OR < 6.63501 & PressQ_OR > 3.8415) { 
			cat('\nPress\'s Q =', PressQ_OR, 'which is > 3.8415, indicating p < .05')
		}
				
		cat('\n\n\nRow Frequencies:\n\n'); print(rowfreqs_OR, print.gap=4) # A frequencies (summed over B) 
	
		cat('\n\nColumn Frequencies:\n\n'); print(colfreqs_OR, print.gap=4) # B frequencies (summed over A)
		
		cat('\n\nCell Proportions:\n\n'); print(round(cellprops_OR,2), print.gap=4)
	
		cat('\n\nRow-Based Proportions:\n\n'); print(round(rowprops_OR,2), print.gap=4) 
	
		cat('\n\nColumn-Based Proportions:\n\n'); print(round(colprops_OR,2), print.gap=4) 
		
		cat('\n\nAgreement (kappas) between the Predicted and Original Group Memberships:\n\n'); print(kappas_cvo_OR,3, print.gap=4)
			
		# Frequencies: Original vs Cross-Validated (leave-one-out cross-validation)
		
		cat('\n\n\n\nCross-Tabulation of the Cross-Validated and Predicted Group Memberships:\n\n'); print(freqs_CV, print.gap=4)	
		
		cat('\n\nProportion of cross-validated grouped cases correctly classified:  ', round(PropCrossValCorrect,3), print.gap=4)
		
		cat('\n\nChi-square test of indepedence:\n'); print(chi_square_CV) # chi-square test of indepedence
		
		cat('\n\nPress\'s Q significance test of classifiation accuracy:\n')
		
		if (PressQ_CV < 3.8415) { 
			cat('\nPress\'s Q =', PressQ_CV, ', which is < 3.8415, indicating non-significance')
		} else if (PressQ_CV > 6.63501) { 
			cat('\nPress\'s Q =', PressQ_CV, 'which is > 6.63501, indicating p < .01')
		} else if (PressQ_CV < 6.63501 & PressQ_CV > 3.8415) { 
			cat('\nPress\'s Q =', PressQ_CV, 'which is > 3.8415, indicating p < .05')
		}
				
		cat('\n\n\nRow Frequencies:\n\n'); print(rowfreqs_CV, print.gap=4) # A frequencies (summed over B) 
	
		cat('\n\nColumn Frequencies:\n\n'); print(colfreqs_CV, print.gap=4) # B frequencies (summed over A)
		
		cat('\n\nCell Proportions:\n\n'); print(round(cellprops_CV,2), print.gap=4)
	
		cat('\n\nRow-Based Proportions:\n\n'); print(round(rowprops_CV,2), print.gap=4)
		
		cat('\n\nColumn-Based Proportions:\n\n'); print(round(colprops_CV,2), print.gap=4) 
		
		cat('\n\nAgreement (kappas) between the Cross-Validated and Original Group Memberships:\n\n')
		print(kappas_CVO, print.gap=4)
		
		cat('\n\nAgreement (kappas) between the Cross-Validated and Predicted Group Memberships:\n\n')
		print(kappas_CVP, print.gap=4)
		
		cat('\n\n\n')	
		
		DFAoutput$ldaoutput_CV <- ldaoutput_CV
		DFAoutput$freqs_OR <- freqs_OR
		DFAoutput$PropOrigCorrect <- PropOrigCorrect
		DFAoutput$chi_square_OR <- chi_square_OR
		DFAoutput$PressQ_OR <- PressQ_OR
		DFAoutput$rowfreqs_OR <- rowfreqs_OR
		DFAoutput$colfreqs_OR <- colfreqs_OR
		DFAoutput$cellprops_OR <- cellprops_OR
		DFAoutput$rowprops_OR <- rowprops_OR
		DFAoutput$colprops_OR <- colprops_OR
		DFAoutput$kappas_cvo_OR <- kappas_cvo_OR
		DFAoutput$freqs_CV <- freqs_CV
		DFAoutput$PropCrossValCorrect <- PropCrossValCorrect
		DFAoutput$chi_square_CV <- chi_square_CV
		DFAoutput$PressQ_CV <- PressQ_CV
		DFAoutput$rowfreqs_CV <- rowfreqs_CV
		DFAoutput$colfreqs_CV <- colfreqs_CV
		DFAoutput$cellprops_CV <- cellprops_CV
		DFAoutput$rowprops_CV <- rowprops_CV
		DFAoutput$colprops_CV <- colprops_CV
		DFAoutput$kappas_CVO <- kappas_CVO
		DFAoutput$kappas_CVP <- kappas_CVP		
	}
}

return(invisible(DFAoutput))

}


