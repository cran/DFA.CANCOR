
DFA <- function(data, groups, variables, plot=TRUE,
                predictive=TRUE, priorprob='SIZES', covmat_type='within', CV=TRUE, verbose=TRUE) {

donnes <- as.data.frame(data[,c(groups,variables)])


if (anyNA(donnes) == TRUE) {
	donnes <- na.omit(donnes)
#	message('\n\nCases with missing values were found and removed from the data matrix.\n')
	NAflag = TRUE
} else {
	NAflag = FALSE
}


if (is.factor(donnes[,1])) {
	grp_levels <- levels(donnes[,1])
	# grp_labels <- labels(donnes[,1])	
}

if (!is.factor(donnes[,1])) {
	donnes[,1] <- factor(donnes[,1], ordered = FALSE)
	grp_levels <- levels(donnes[,1])
}


grpnames <- grp_levels


# grpnames <- as.vector(as.matrix(donnes[groups])) # group names, in the same order as in the data matrix
# grpnames <- unique(grpnames)
grpnums  <- seq(1:length(grpnames))

Ngroups  <- length(grpnames)
Ndvs <- ncol(donnes) - 1
Ndfs <- min( (Ngroups-1), Ndvs)
Ncases <- nrow(donnes)


# if (!is.factor(donnes[,1]))  donnes[,1] <- factor(donnes[,1], ordered = FALSE, labels=grpnames, levels=grpnames)

donnes <- as.data.frame(donnes)

DVnames <- colnames(donnes[,2:(Ndvs+1)])


groupNs <- table(donnes[,1])
groupNs <- groupNs[grpnames] # because table gives results in alphabetical order

# groupNs <- as.matrix(table(donnes[,1]))
	# # Group Sizes
	# # groupNs <- aggregate(x = DFAposteriors[,postvarnoms[1]], by = list(DFAposteriors$Group), FUN=length)[,2]
	# # names(groupNs) <- grpnames
	# groupNs <- table(DFAposteriors$Group)
# groupNs[order(match(groupNs,grpnames))]


grpmeans <- sapply(2:ncol(donnes), function(x) tapply(donnes[,x], INDEX = donnes[,1], FUN = mean))
grpmeans <- grpmeans[grpnames,] # resorting rows because mean gives results in alphabetical order





# from my homovarcovar function:
# for the var-covar matrices from SPSS - requires Type III sums of squares
# using MANOVA to obtain the sums of squares and cross-products matrix for error, &
# the sums of squares and cross-products matrix for group/IV
# www.webpages.uidaho.edu/~kirk/calendar/R/MANOVA.doc
MV2 <- manova( as.matrix(donnes[,2:ncol(donnes)]) ~ donnes[,1], data=donnes)
sscpwith <- (Ncases-1) * cov(MV2$residuals) # E
sscpbetw <- (Ncases-1) * cov(MV2$fitted.values) # H


W <- sscpwith * (1/(Ncases-Ngroups))   # vcv 

B <- 1/(3-1) * sscpbetw   # vcv 

W.inv <- solve(W)

A <- W.inv %*% B   # the canonical matrix 

A.vectors <- eigen(A)$vectors[,1:Ndfs]  # the # of cols matters

# Complex numbers  -- http://blog.phytools.org/2013/07/complex-numbers-in-r-and-new-version-of.html
# It turns out that part spectral decomposition (eigenanalysis) can result in very slightly negative 
# eigenvalues, which result in complex numbers during some steps of the calculations of CCA. In later 
# steps, the imaginary parts go to zero - but R keeps the class of a number as "complex" even if the 
# imaginary part is zero.
# To fix this, I merely check to see if all imaginary parts of the canonical correlations are zero, 
# and (if so) set the canonical correlation to the real part of the number only:
if(all(Im(A.vectors)==0)) A.vectors <- Re(A.vectors)


# *** DA_lab  Utente Discriminant Analysis - 3rd TUTORIAL
# The A matrix (here indicated as 'Coefficients of linear discriminants') yielded by the lda R function 
# is a bit different from what we obtained with the single value decomposition of matrix: here it is 
# normalized so that the within groups covariance matrix is spherical.
# By multiplying A vectors for psi we find again the coefficients of the linear discriminants yield by lda:
# MASS scaling:
# a matrix which transforms observations to discriminant functions, normalized so that within groups covariance matrix is spherical

psi <- t(A.vectors) %*% W %*% A.vectors

psiINV <- solve(psi)
psiINV <- ifelse (psiINV < .000000001, .0000001, psiINV) # setting vv small or neg values to a small #
coefs_raw <- A.vectors %*% (psiINV^(1/2))

colnames(coefs_raw) <-  c(paste("Function ", 1:Ndfs, sep=""))
rownames(coefs_raw) <-  DVnames


# scores on the LDFs
# https://stackoverflow.com/questions/68307682/r-lda-linear-discriminant-analysis-how-to-get-compute-lda-scores-from-lda-co

# 1 -- Calculate the group means for each variable
# grpmeans <- sapply(2:(Ndvs+1), function(x) tapply(donnes[,x], INDEX = donnes[,1], FUN = mean))

# 2 -- Calculate the mean of group means for each variable -- BOC = must be weighted grp means, by grp Ns
# center <- colMeans(grpmeans)
# groupNs <- table(donnes[,1])
center <- apply(grpmeans, 2, function(m, w) { sum(m * w) / sum(w) }, w=groupNs)

# 3 -- Center the data at the mean of group means
# The x argument in scale can be the original data, or any new data one wants to project (predict) 
# into the fitted discriminant space. However, one always has to use the centering vector defined 
# by the original data (used for LDA model fitting, center in our example) to center new data accordingly.
DVs_centered <- scale(x = donnes[,2:ncol(donnes)], center = center, scale = FALSE)

# 4 -- multiply the centered data by the LD coefficients to get the actual scores
dfa_scores <- DVs_centered %*% coefs_raw

# add the group values to the first column 
dfa_scores <- data.frame(donnes[,1], dfa_scores); colnames(dfa_scores)[1] <- 'group'



# eigenvalues, canonical correlations, & one-way anovas on the DFs
evals <- matrix(-9999, Ndfs, 4)
evals[,1] <- 1:Ndfs
anovaDFoutput <- matrix(-9999, Ndfs, 8)
ttestDFoutput <- lapply(1:Ndfs, function(x) matrix(-9999, nrow=choose(Ngroups,2), ncol=13))
names(ttestDFoutput) = c(paste("Discriminant Function ", 1:Ndfs, sep=""))
for (luper in 1:Ndfs) {
	dd <- data.frame(dfa_scores[,1], dfa_scores[luper+1])
	colnames(dd) <- c('grp','dv')
	fit <- stats::lm(dv ~ as.factor(grp), data = dd)
	betwss <- stats::anova(fit)["as.factor(grp)", "Sum Sq"]
	withss <- stats::anova(fit)["Residuals", "Sum Sq"]
 	anovaDFoutput[luper,1] <- as.numeric(summary(fit)["r.squared"])
 	anovaDFoutput[luper,2] <- 1 - as.numeric(summary(fit)["r.squared"])
	anovaDFoutput[luper,3] <- stats::anova(fit)["as.factor(grp)","F value"]
	anovaDFoutput[luper,4] <- stats::anova(fit)["as.factor(grp)","Df"]
	anovaDFoutput[luper,5] <- fit$df.residual
	anovaDFoutput[luper,6] <- stats::anova(fit)["as.factor(grp)","Pr(>F)"]

	anBF <- suppressMessages(anovaBF(dv ~ grp, data = dd, progress=FALSE))
	anovaDFoutput[luper,7] <- as.numeric(suppressMessages(extractBF(anBF)[1]))

	anovaDFoutput[luper,8] <- 1 / as.numeric(suppressMessages(extractBF(anBF)[1]))

	ttestDFoutput[[luper]] <- GROUP.DIFFS(dd, var.equal=FALSE, verbose=FALSE)			

	evals[luper,2] <- betwss / withss # eigenvalue
	
	evals[luper,4] <- sqrt(betwss / (betwss+withss)) # canonical correlation
}

evals[,3] <- evals[,2] / sum(evals[,2])   # proportions of variance
	
dimnames(evals) <- list(rep("", dim(evals)[1]))
colnames(evals) <- c('Function','eigenvalue','proportion of variance','canonical r')


# multivariate tests

# July, 2019: for CCA, the CCP::p.asym function uses the #s of vars in set1 and set2
# as the values for p & q; but to get correct results for DFA using this function,
# p = the # of DVs and q = the # of DFs 

mv_Wilks <- Wilks(rho=evals[,4], Ncases=Ncases, p = Ndfs, q = Ndvs)
colnames(mv_Wilks) <- c('Wilk\'s Lambda', 'F-approx.', 'df1', 'df2', 'p')
rownames(mv_Wilks) <- paste(1:nrow(mv_Wilks), paste("through ", nrow(mv_Wilks), sep = ""))
# print(round_boc(mv_Wilks), print.gap=4); message('\n')

mv_Pillai <- Pillai(rho=evals[,4], Ncases=Ncases, p = Ndfs, q = Ndvs)
colnames(mv_Pillai) <- c('Pillai-Bartlett Trace', 'F-approx.', 'df1', 'df2', 'p')
rownames(mv_Pillai) <- paste(1:nrow(mv_Pillai), paste("through ", nrow(mv_Pillai), sep = ""))
# print(round_boc(mv_Pillai), print.gap=4); message('\n')

mv_Hotelling <- Hotelling(rho=evals[,4], Ncases=Ncases, p = Ndfs, q = Ndvs)
colnames(mv_Hotelling) <- c('Hotelling-Lawley Trace', 'F-approx.', 'df1', 'df2', 'p')
rownames(mv_Hotelling) <- paste(1:nrow(mv_Hotelling), paste("through ", nrow(mv_Hotelling), sep = ""))
# print(round_boc(mv_Hotelling), print.gap=4); message('\n')

mv_Roy <- RoyRoot(rho=evals[,4], Ncases=Ncases, p = Ndfs, q = Ndvs)
colnames(mv_Roy) <- c('Roy\'s Largest Root', 'lambda ', 'F-approx. ', 'df1', 'df2', 'p')
rownames(mv_Roy) <- paste(1:nrow(mv_Roy), paste("through ", nrow(mv_Roy), sep = ""))
# print(round_boc(mv_Roy), print.gap=4); message('\n')



# centering each DV within groups
group.center <- function(var,grp) { return(var - tapply(var,grp,mean,na.rm=TRUE)[grp]) }
centd_donnes <- matrix(-9999,nrow(donnes),(ncol(donnes)-1))
for (lupec in 1:(Ndvs)) { centd_donnes[,lupec] <- group.center(donnes[,(lupec+1)], donnes[,1]) }

# centering each DF within groups
centd_dfa_scores <- matrix(-9999, Ncases, Ndfs)
for (lupec in 2:(Ndfs+1)) { centd_dfa_scores[,(lupec-1)] <- group.center(dfa_scores[,lupec], dfa_scores[,1]) }

# the structure coefficients
coefs_structure <- stats::cor(x = centd_donnes, y = centd_dfa_scores) 
rownames(coefs_structure) <- rownames(coefs_raw)
colnames(coefs_structure) <- c(paste("Function ", 1:ncol(coefs_structure), sep=""))

# the standardized coefficients
centd_donnes <- data.frame(cbind(donnes[,1], centd_donnes)) # placing the grouping variable back in the centered matrix
pooledSDs <- as.matrix(apply(centd_donnes, 2, FUN = sd)) # centd_donnes contains the group mean-centered data
coefs_standardized <- (pooledSDs[2:nrow(pooledSDs),]) * coefs_raw
colnames(coefs_standardized) <-  c(paste("Function ", 1:ncol(coefs_standardized), sep=""))

# the SPSS standardized coefficients 
poolwith <- sscpwith * (1/(nrow(donnes)-Ngroups))
pooledSDs <- sqrt(diag(poolwith)) # pooled SDs for SPSS results
coefs_standardizedSPSS <- pooledSDs * coefs_raw
colnames(coefs_standardizedSPSS) <- c(paste("Function ", 1:ncol(coefs_standardizedSPSS), sep=""))


# Centroids	
# group means & SDs on the raw ldfs
centroids   <- sapply(2:(Ndfs+1), function(x) tapply(dfa_scores[,x], INDEX = dfa_scores[,1], FUN = mean))
centroidsSDs <- sapply(2:(Ndfs+1), function(x) tapply(dfa_scores[,x], INDEX = dfa_scores[,1], FUN = sd))
# group means & SDs on the standardized ldfs
dfa_scoresZ <- data.frame(dfa_scores[,1], scale(dfa_scores[,2:(Ndfs+1)]))
centroidsZ   <- sapply(2:(Ndfs+1), function(x) tapply(dfa_scoresZ[,x], INDEX = dfa_scoresZ[,1], FUN = mean))
centroidsSDsZ <- sapply(2:(Ndfs+1), function(x) tapply(dfa_scoresZ[,x], INDEX = dfa_scoresZ[,1], FUN = sd))
colnames(centroids) <- colnames(centroidsSDs) <- colnames(centroidsZ) <- colnames(centroidsSDsZ) <- colnames(coefs_raw)

# rownames(centroids) <- rownames(centroidsSDs) <- rownames(centroidsZ) <- rownames(centroidsSDsZ) <- grpnames  -- prob

# resorting rows to be consistent with the groups appearing in the same order as they do in donnes
centroids     <- centroids[grpnames,] 
centroidsSDs  <- centroidsSDs[grpnames,] 
centroidsZ    <- centroidsZ[grpnames,] 
centroidsSDsZ <- centroidsSDsZ[grpnames,] 

if (Ndfs == 1) {
	centroids     <- matrix(centroids,     length(centroids),     1)  
	centroidsSDs  <- matrix(centroidsSDs,  length(centroidsSDs),  1)  
	centroidsZ    <- matrix(centroidsZ,    length(centroidsZ),    1)
	centroidsSDsZ <- matrix(centroidsSDsZ, length(centroidsSDsZ), 1)  
	
	colnames(centroids) <- colnames(centroidsSDs) <- colnames(centroidsZ) <- colnames(centroidsSDsZ) <- "Function 1"
	rownames(centroids) <- rownames(centroidsSDs) <- rownames(centroidsZ) <- rownames(centroidsSDsZ) <- grpnames
}



# one-way ANOVAs using the scores on a discriminant function as the DV
dimnames(anovaDFoutput) <-list(rep("", dim(anovaDFoutput)[1]))
colnames(anovaDFoutput) <- c('Eta-squared','Wilks_Lambda','F','df','df_res','p','Bayes_Factor_alt_vs_null','Bayes_Factor_null_vs_alt')
rownames(anovaDFoutput) <- colnames(coefs_raw)
# anovaDFoutput[,1:4] <- anovaDFoutput[,1:4]
# anovaDFoutput[,5] <- anovaDFoutput[,5]


# one-way anovas & t-tests on the DVs
anovaDVoutput <- matrix(-9999, length(variables), 8)
ttestDVoutput <- lapply(1:length(variables), function(x) matrix(-9999, nrow=choose(Ngroups,2), ncol=13))
names(ttestDVoutput)=variables 
for (lupec in 1:length(variables)) {
	ddd <- data.frame(donnes[,1], donnes[,(lupec+1)])
	colnames(ddd) <- c('grp','dv')
	fit <- stats::lm(dv ~ as.factor(grp), data = ddd)
	betwss <- stats::anova(fit)["as.factor(grp)", "Sum Sq"]
	withss <- stats::anova(fit)["Residuals", "Sum Sq"]
 	anovaDVoutput[lupec,1] <- as.numeric(summary(fit)["r.squared"])
 	anovaDVoutput[lupec,2] <- 1 - as.numeric(summary(fit)["r.squared"])
	anovaDVoutput[lupec,3] <- stats::anova(fit)["as.factor(grp)","F value"]
	anovaDVoutput[lupec,4] <- stats::anova(fit)["as.factor(grp)","Df"]
	anovaDVoutput[lupec,5] <- fit$df.residual
	anovaDVoutput[lupec,6] <- stats::anova(fit)["as.factor(grp)","Pr(>F)"]

	anBF <- suppressMessages(anovaBF(dv ~ grp, data = ddd, progress=FALSE))
	anovaDVoutput[lupec,7] <- as.numeric(extractBF(anBF)[1])

	anovaDVoutput[lupec,8] <- 1 / as.numeric(extractBF(anBF)[1])
	
	ttestDVoutput[[lupec]] <- GROUP.DIFFS(ddd, var.equal=FALSE, verbose=FALSE)			
}
dimnames(anovaDVoutput) <-list(rep("", dim(anovaDVoutput)[1]))
colnames(anovaDVoutput) <- c('Eta-squared','Wilks_Lambda','F','df','df_res','p','Bayes_Factor_alt_vs_null','Bayes_Factor_null_vs_alt')
rownames(anovaDVoutput) <- variables
# anovaDVoutput[,1:4] <- anovaDVoutput[,1:4]
# anovaDVoutput[,5] <- anovaDVoutput[,5]


DFAoutput <- list(  
   evals = evals,
   mv_Wilks =mv_Wilks,
   mv_Pillai = mv_Pillai,
   mv_Hotelling = mv_Hotelling,
   mv_Roy = mv_Roy,
   coefs_raw=coefs_raw,
   coefs_structure=coefs_structure,
   coefs_standardized=coefs_standardized,
   coefs_standardizedSPSS=coefs_standardizedSPSS,
   centroids=centroids,
   centroidsSDs=centroidsSDs,
   centroidsZ=centroidsZ,
   centroidsSDsZ=centroidsSDsZ,
   dfa_scores=dfa_scores,
   anovaDFoutput = anovaDFoutput,
   anovaDVoutput = anovaDVoutput
)



# 2014 Bird - Controlling the Maximum Familywise Type I Error Rate in Analyses of Multivariate Experiments

lsnoms <- names(ttestDVoutput)

if (Ngroups == 2) {
	# Two Groups 
	# Critical values for analyses controlling MFWER when k = 2 
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
	
	

if (predictive | is.null(predictive)) {

	classes_PRED <- DFA_classes(donnes=donnes, grpmeans=grpmeans, Ngroups=Ngroups, groupNs=groupNs, grpnames=grpnames,
	                            Ncases=Ncases, W=W, priorprob=priorprob, covmat_type=covmat_type )

	freqs_ORIG_PRED <- squareTable(var1=donnes[,1], var2=classes_PRED$dfa_class, faclevels=grpnames, tabdimnames=c('Original','Predicted'))
	
	chi_square_ORIG_PRED <- summary(freqs_ORIG_PRED)

	PressQ_ORIG_PRED <- PressQ(freqs_ORIG_PRED)

	kappas_ORIG_PRED <- kappas(var1=donnes[,1], var2=classes_PRED$dfa_class, grpnames)
	
	PropOrigCorrect <- round((sum(diag(freqs_ORIG_PRED)) / sum(freqs_ORIG_PRED)),3)	
	
	# rowfreqs_ORIG_PRED <- margin.table(freqs_ORIG_PRED, 1)
	# colfreqs_ORIG_PRED <- margin.table(freqs_ORIG_PRED, 2)
	# cellprops_ORIG_PRED <- prop.table(freqs_ORIG_PRED)
	# rowprops_ORIG_PRED  <- prop.table(freqs_ORIG_PRED, 1)
	# colprops_ORIG_PRED  <- prop.table(freqs_ORIG_PRED, 2)
	
	# posterior probability stats
	grp_post_stats <- DFC_post_class_stats(classes_PRED$posteriors, grpnames=grpnames, 
	                                       Ngroups=Ngroups, groupNs=groupNs, verbose=FALSE) 


	# faclevels <- levels(donnes[,1])
		
	if (CV) {

		# Original vs Cross-Validated (leave-one-out cross-validation)
		
		# classifications from leave-one-out cross-validation
		classes_CV <- DFA_classes_CV(donnes=donnes, priorprob=priorprob, covmat_type=covmat_type, grpnames=grpnames )
	
		freqs_ORIG_CV <- squareTable(var1=donnes[,1], var2=classes_CV, 
		                             faclevels=grpnames, tabdimnames=c('Original','Cross-Validated'))
		# freqs_ORIG_CV <- squareTable(var1=classes_PRED$dfa_class, var2=classes_CV, tabdimnames=c('Original','Cross-Validated'))
	
		colnames(freqs_ORIG_CV) <- paste(grpnames)
		rownames(freqs_ORIG_CV) <- paste(grpnames)

		chi_square_ORIG_CV <- summary(freqs_ORIG_CV)

		PressQ_ORIG_CV <- PressQ(freqs_ORIG_CV)
	
		kappas_ORIG_CV <- kappas(donnes[,1], classes_CV, grpnames)

		PropCrossValCorrect <- sum(diag(freqs_ORIG_CV)) / sum(freqs_ORIG_CV)			
		
		# rowfreqs_ORIG_CV <- margin.table(freqs_ORIG_CV, 1)
		# colfreqs_ORIG_CV <- margin.table(freqs_ORIG_CV, 2)		
		# cellprops_ORIG_CV <- prop.table(freqs_ORIG_CV)
		# rowprops_ORIG_CV  <- prop.table(freqs_ORIG_CV, 1)
		# colprops_ORIG_CV  <- prop.table(freqs_ORIG_CV, 2)

		
		# Predicted vs Cross-Validated (leave-one-out cross-validation)


		# attributes(var2) <- attributes(var1) 


# attributes(freqs_ORIG_CV)


		# freqs_PRED_CV <- squareTable(var1=classes_PRED$dfa_class, var2=classes_CV, 
		                             # faclevels=grpnames, tabdimnames=c('Predicted','Cross-Validated'))

		# colnames(freqs_PRED_CV) <- paste(grpnames)
		# rownames(freqs_PRED_CV) <- paste(grpnames)

		# chi_square_PRED_CV <- summary(freqs_PRED_CV)

		# PressQ_PRED_CV <- PressQ(freqs_PRED_CV)

		# kappas_PRED_CV <- kappas(var1 = classes_PRED$dfa_class, var2 = classes_CV, grpnames)

	}

	DFAoutput$classes_PRED <- classes_PRED$dfa_class
	if (CV) DFAoutput$classes_CV <- classes_CV
	DFAoutput$posteriors <- classes_PRED$posteriors 
	DFAoutput$grp_post_stats <- grp_post_stats
	DFAoutput$freqs_ORIG_PRED <- freqs_ORIG_PRED
	DFAoutput$chi_square_ORIG_PRED <- chi_square_ORIG_PRED
	DFAoutput$PressQ_ORIG_PRED <- PressQ_ORIG_PRED
	DFAoutput$kappas_ORIG_PRED <- kappas_ORIG_PRED
	DFAoutput$PropOrigCorrect <- PropOrigCorrect

	# DFAoutput$rowfreqs_ORIG_PRED <- rowfreqs_ORIG_PRED
	# DFAoutput$colfreqs_ORIG_PRED <- colfreqs_ORIG_PRED
	# DFAoutput$cellprops_ORIG_PRED <- cellprops_ORIG_PRED
	# DFAoutput$rowprops_ORIG_PRED <- rowprops_ORIG_PRED
	# DFAoutput$colprops_ORIG_PRED <- colprops_ORIG_PRED


	DFAoutput$freqs_ORIG_CV <- freqs_ORIG_CV
	DFAoutput$chi_square_ORIG_CV <- chi_square_ORIG_CV
	DFAoutput$PressQ_ORIG_CV <- PressQ_ORIG_CV
	DFAoutput$kappas_ORIG_CV <- kappas_ORIG_CV
	DFAoutput$PropCrossValCorrect <- PropCrossValCorrect

	# DFAoutput$rowfreqs_ORIG_CV <- rowfreqs_ORIG_CV
	# DFAoutput$colfreqs_ORIG_CV <- colfreqs_ORIG_CV
	# DFAoutput$cellprops_ORIG_CV <- cellprops_ORIG_CV
	# DFAoutput$rowprops_ORIG_CV <- rowprops_ORIG_CV
	# DFAoutput$colprops_ORIG_CV <- colprops_ORIG_CV

	# DFAoutput$freqs_PRED_CV <- freqs_PRED_CV
	# DFAoutput$chi_square_PRED_CV <- chi_square_PRED_CV
	# DFAoutput$PressQ_PRED_CV <- PressQ_PRED_CV
	# DFAoutput$kappas_PRED_CV <- kappas_PRED_CV
	
}



# plot
if (plot == TRUE) {

	# if (Ndfs == 1) {
		# centroidsZ_2 <- matrix(centroidsZ, length(centroidsZ), 1)
		# colnames(centroidsZ_2) <- "DF 1"
	# }
	# if (Ndfs > 1) {
		# centroidsZ_2 <- centroidsZ		
		# colnames(centroidsZ_2) <-  c(paste("DF ", 1:ncol(centroidsSDsZ), sep=""))		
	# }


	centroidsZ_2 <- centroidsZ		
	colnames(centroidsZ_2) <-  c(paste("DF ", 1:ncol(centroidsSDsZ), sep=""))		

	graphics::matplot(1:length(grpnames), centroidsZ_2, type = "l", lty=1, lwd=3, 
	        xaxt='n', xlab=groups, cex.axis=1.2, cex.lab = 1.3,
	        ylab='Discriminant Function z Scores', ylim = c(-3, 3), cex.axis=1.2        )
	graphics::axis(side=1, at=grpnums, labels=grpnames, xlab="groups")   # labels=rownames(centroidsZ_2)
	graphics::title(main='Mean Standardized Discriminant Function Scores for the Groups')
	graphics::legend("topright", legend = colnames(centroidsZ_2), bty="n", lwd=2, col=1:ncol(centroidsZ_2))
}



if (verbose == TRUE) {
	
	message('\n\nLinear Discriminant Function Analysis')

	if (NAflag) message('\n\nCases with missing values were found and removed from the data matrix.\n')
			
	# eigenvalues, canonical correlations, & one-way anovas on the DFs
	message('\n\nEigenvalues & canonical correlations:\n')
	print(round_boc(evals), print.gap=4)
	
	message('\n\nMultivariate peel-down significance tests:\n')  # using p.asym from the CCP package
	
	print(round_boc(mv_Wilks), print.gap=4); message('\n')

	print(round_boc(mv_Pillai), print.gap=4); message('\n')
	
	print(round_boc(mv_Hotelling), print.gap=4); message('\n')
	
	print(round_boc(mv_Roy), print.gap=4)
	
	message('\n\nCanonical Discriminant Function (raw) Coefficients:\n')
	print(round_boc(coefs_raw), print.gap=4)
	
	message('\n\nStructure Coefficients:\n')
	print(round_boc(coefs_structure), print.gap=4)
	
	message('\n\nStandardized Coefficients:\n')
	print(round_boc(coefs_standardized), print.gap=4)
	
	message('\n\nStandardized Coefficients from SPSS:\n')
	print(round_boc(coefs_standardizedSPSS), print.gap=4)
	
	message('\n\nFunctions at Group Centroids\n')
	message('\nUnstandardized canonical discriminant functions evaluated at group means:\n')
	print(round_boc(centroids), print.gap=4)
	
	message('\n\nGroup Standard Deviations on the unstandardized functions:\n')
	print(round_boc(centroidsSDs), print.gap=4)
	
	message('\n\nStandardized canonical discriminant functions evaluated at group means:\n')
	print(round_boc(centroidsZ), print.gap=4)
	
	message('\n\nGroup Standard Deviations on the standardized functions:\n')
	print(round_boc(centroidsSDsZ), print.gap=4)
	
	message('\n\nOne-way ANOVAs using the scores on a discriminant function as the DV:\n')
	print(round_boc(anovaDFoutput), print.gap=4)
	
	message('\n\nOne-way ANOVAs on the original DVs') 
	message('\n(provided for comparisons with the ANOVAs on the discriminant functions;')
	message(' In SPSS, this is labelled, "Tests of Equality of Group Means":)\n')
	print(round_boc(anovaDVoutput), print.gap=4)
	
	message('\n\n\nt-tests and effect sizes for group differences on the discriminant functions\n\n')
	DFnames <- names(ttestDFoutput)
	for (lupe in 1:length(ttestDFoutput)) {
		# message(DFnames[lupe],'\n')
		# print(round_boc(ttestDFoutput[[lupe]]), row.names = FALSE, print.gap=3); message('\n')

		message(lsnoms[lupe],':',sep="")

		resultsM <- ttestDFoutput[[lupe]]

		message("\nGroup comparisons - significance tests\n")
		print(round_boc(resultsM[1:12], round_non_p = 2), print.gap=3)
		
		message("\nGroup comparisons - confidence intervals, effect sizes, and Bayes Factors\n")
		print(round_boc(resultsM[c(1,5,16:21)], round_non_p = 2), print.gap=3); message('\n')
	}
	
	message('\nt-tests and effect sizes for group differences on the original DVs\n')
	message('(provided for comparisons with the t-tests on the discriminant functions)\n\n')
	for (lupe in 1:length(ttestDVoutput)) {
		message(lsnoms[lupe],':',sep="")
		# print(round_boc(ttestDVoutput[[lupe]]), row.names = FALSE, print.gap=3); message('\n')

		resultsM <- ttestDVoutput[[lupe]]

		message("\nGroup comparisons - significance tests:\n")
		print(round_boc(resultsM[1:12], round_non_p = 2), print.gap=3)
		
		message("\nGroup comparisons - confidence intervals, effect sizes, and Bayes Factors:\n")
		print(round_boc(resultsM[c(1,5,16:21)], round_non_p = 2), print.gap=3); message('\n')
	}
	
	
			
	# 2014 Bird - Controlling the Maximum Familywise Type I Error Rate in Analyses of Multivariate Experiments
	
	message('\nAppropriate univariate significance tests (ignoring the discriminant functions)')
	message('for controlling the maximum familywise Type I error rate (MFWER), based on:')
	message('\nBird, K. D., & Hadzi-Pavlovic, D. (2014). Controlling the maximum familywise Type I ')
	message('error rate in analyses of multivariate experiments. Psychological Methods, 19(2), 265-280.\n')
	
	if (Ngroups == 2) {
		
		message('\nWhen there are only two groups and the researcher decides to control the MFWER') 
		message('at .05 by (only) carrying out multiple t tests, then the MANOVA test is')
		message('irrelevant and the critical t value is the t1-critical value below.')
		message('The t value for the data must be greater than t1-critical value for a pairwise')
		message('comparison to be significant.\n\n')
						
		for (lupe in 1:length(ttestDVoutput)) {
			message(lsnoms[lupe],'\n')
			print(round_boc(DFAoutput[[variables[lupe]]]$MFWER1.sigtest), row.names = FALSE, print.gap=4); message('\n')
		}		
	
		message('\nWhen there are only two groups and the researcher prefers the two-stage')
		message('approach to controling the MFWER at .05, then the MANOVA F test must be')
		message('significant and the t value for each comparison must be greater than') 
		message('t2-critical value below for a pairwise comparison to be significant.\n\n')
	
		for (lupe in 1:length(ttestDVoutput)) {
			message(lsnoms[lupe],'\n')
			print(round_boc(DFAoutput[[variables[lupe]]]$MFWER2.sigtest), row.names = FALSE, print.gap=4)
		}		
	}
		
	if (Ngroups > 2) {
			
		message('\nFor designs with more than two groups, Bird et al. (2014) could find no')
		message('justification for the inclusion of an initial MANOVA test in MCPs designed to') 
		message('control the MFWER with protected t tests.')
		
		message('\nWhen a researcher decides to control the MFWER at .05 by solely carrying out') 
		message('multiple t tests, then the MANOVA and ANOVA tests are irrelevant and the')
		message('critical t value is the t1-critical value below. The t value for the data must')
		message('be greater than t1-critical value for a pairwise comparison to be significant.\n\n')
			
		for (lupe in 1:length(ttestDVoutput)) {
			message(lsnoms[lupe],'\n')
			print(round_boc(DFAoutput[[variables[lupe]]]$MFWER1.sigtest), row.names = FALSE, print.gap=4)
		}		
	
		message('\nWhen a researcher prefers the two-stage approach to controlling the MFWER at .05,')
		message('then the ANOVA F value for the data must be greater than the F-critical value below,')
		message('and the t value for each comparison must be greater than t2-critical value') 
		message('below for a pairwise comparison to be significant.\n\n')
	
		for (lupe in 1:length(ttestDVoutput)) {
			message(lsnoms[lupe],'\n')
			message('\nF = ', round(anovaDVoutput[lupe,2],2),'   F-critical = ', 
			        round(FCRIT,2),'   decision =', decisionF[lupe],'\n')
			print(round_boc(DFAoutput[[variables[lupe]]]$MFWER2.sigtest), row.names = FALSE, print.gap=4); message('\n')
		}		
	}

		
	if (predictive == TRUE | is.null(predictive)) {
	
		message('\n\nPREDICTIVE DISCRIMINANT ANALYSIS')
		
		message('\n\nPrior Probabilities for Groups:\n')
		print(round(classes_PRED$prior,3), print.gap=4, row.names=FALSE)
				
		if (covmat_type == 'within') {
			message('\n\nClassification Function Coefficients:\n')
			print(round(classes_PRED$classifcoefs,3), print.gap=4, row.names=FALSE)

			message('\n\nClassification Function Intercepts:\n')
			print(round(classes_PRED$classifints,3), print.gap=4, row.names=FALSE)
		}
				
		message('\n\nCross-Tabulation of the Original and Predicted Group Memberships:\n')
		print(freqs_ORIG_PRED, print.gap=4)	
		
		message('\n\nProportion of original grouped cases correctly classified:  ', round(PropOrigCorrect,3))
		
		message('\n\nChi-square test of independence:\n')
		print(chi_square_ORIG_PRED, print.gap=4) 

		message('\n\nPress\'s Q significance test of classifiation accuracy:')
		
		if (PressQ_ORIG_PRED < 3.8415) { 
			message('\nPress\'s Q = ', round(PressQ_ORIG_PRED,2), ', which is < 3.8415, indicating non-significance')
		} else if (PressQ_ORIG_PRED > 6.63501) { 
			message('\nPress\'s Q = ', round(PressQ_ORIG_PRED,2), ' which is > 6.63501, indicating p < .01')
		} else if (PressQ_ORIG_PRED < 6.63501 & PressQ_ORIG_PRED > 3.8415) { 
			message('\nPress\'s Q = ', round(PressQ_ORIG_PRED,2), ' which is > 3.8415, indicating p < .05')
		}
				
		# message('\n\nRow Frequencies:\n'); print(rowfreqs_ORIG_PRED, print.gap=4) # A frequencies (summed over B) 
	
		# message('\n\nColumn Frequencies:\n'); print(colfreqs_ORIG_PRED, print.gap=4) # B frequencies (summed over A)
		
		# message('\n\nCell Proportions:\n'); print(round(cellprops_ORIG_PRED,2), print.gap=4)
	
		# message('\n\nRow-Based Proportions:\n'); print(round(rowprops_ORIG_PRED,2), print.gap=4) 
	
		# message('\n\nColumn-Based Proportions:\n'); print(round(colprops_ORIG_PRED,2), print.gap=4) 
		
		message('\n\nAgreement (kappas) between the Original and Predicted Group Memberships:\n')
		print(kappas_ORIG_PRED,3, print.gap=4)

		message("\n\nGroup mean posterior classification probabilities: \n"); 
		print(round(grp_post_stats$grpMNprobs,2), print.gap=4)
			
		message("\n\nNumber of cases per level of posterior classification probability:\n") 
		print(grp_post_stats$grp_prob_Ns, print.gap=4)
			
		message("\n\nProportions of cases per level of posterior classification probability:\n") 
		print(round(grp_post_stats$grp_prob_proports,2), print.gap=4)		
			
		if (CV) {

			message('\n\nCross-Tabulation of the Original and Cross-Validated Group Memberships:\n')
			print(freqs_ORIG_CV, print.gap=4)	
			
			message('\n\nProportion of cross-validated grouped cases correctly classified:  ', round(PropCrossValCorrect,3), print.gap=4)
			
			message('\n\nChi-square test of independence:\n')
			print(chi_square_ORIG_CV)
			
			message('\n\nPress\'s Q significance test of classifiation accuracy:')
			
			if (PressQ_ORIG_CV < 3.8415) { 
				message('\nPress\'s Q = ', round(PressQ_ORIG_CV,2), ', which is < 3.8415, indicating non-significance')
			} else if (PressQ_ORIG_CV > 6.63501) { 
				message('\nPress\'s Q = ', round(PressQ_ORIG_CV,2), ' which is > 6.63501, indicating p < .01')
			} else if (PressQ_ORIG_CV < 6.63501 & PressQ_ORIG_CV > 3.8415) { 
				message('\nPress\'s Q = ', round(PressQ_ORIG_CV,2), ' which is > 3.8415, indicating p < .05')
			}
					
			# message('\n\nRow Frequencies:\n'); print(rowfreqs_ORIG_CV, print.gap=4) # A frequencies (summed over B) 
		
			# message('\n\nColumn Frequencies:\n'); print(colfreqs_ORIG_CV, print.gap=4) # B frequencies (summed over A)
			
			# message('\n\nCell Proportions:\n'); print(round(cellprops_ORIG_CV,2), print.gap=4)
		
			# message('\n\nRow-Based Proportions:\n'); print(round(rowprops_ORIG_CV,2), print.gap=4)
			
			# message('\n\nColumn-Based Proportions:\n'); print(round(colprops_ORIG_CV,2), print.gap=4) 
			
			message('\n\nAgreement (kappas) between the Original and Cross-Validated Group Memberships:\n')
			print(kappas_ORIG_CV, print.gap=4)


			# message('\n\nCross-Tabulation of the Predicted and Cross-Validated Group Memberships:\n')
			# print(freqs_PRED_CV, print.gap=4)	
			
			# message('\n\nChi-square test of independence:\n')
			# print(chi_square_ORIG_CV) 
			
			# message('\n\nPress\'s Q significance test of classifiation accuracy:')
			
			# if (PressQ_PRED_CV < 3.8415) { 
				# message('\nPress\'s Q = ', round(PressQ_PRED_CV,2), ', which is < 3.8415, indicating non-significance')
			# } else if (PressQ_PRED_CV > 6.63501) { 
				# message('\nPress\'s Q = ', round(PressQ_PRED_CV,2), ' which is > 6.63501, indicating p < .01')
			# } else if (PressQ_PRED_CV < 6.63501 & PressQ_PRED_CV > 3.8415) { 
				# message('\nPress\'s Q = ', round(PressQ_PRED_CV,2), ' which is > 3.8415, indicating p < .05')
			# }
					
			# message('\n\nAgreement (kappas) between the Predicted and Cross-Validated Group Memberships:\n')
			# print(kappas_PRED_CV, print.gap=4)
			
		}
				
		message('\n\n')	
	}
}

return(invisible(DFAoutput))

}


