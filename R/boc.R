
	




DFC_post_class_stats <- function (DFAposteriors, grpnames, Ngroups, groupNs, verbose=FALSE) {

	Group <- NULL # to avoid a "no visible binding for global variable 'Group'" warning

	postvarnoms <- names(DFAposteriors[,-grep("Group",colnames(DFAposteriors))])	

	# # Group mean classification probabilities -- aggregate produces a list that is too long/tiresome to work with
	# grpMNprobs <- aggregate(x = DFAposteriors[,postvarnoms], by = list(DFAposteriors$Group), FUN=mean)
	# grpMNprobs <- grpMNprobs[,-1]
	# grpMNprobs <- diag(as.matrix(t(sapply(grpMNprobs,c))))		
	# names(grpMNprobs) <- grpnames
	# message("\n\nGroup mean classification probabilities:\n")
	# print(round_boc(grpMNprobs,2), print.gap=4)

	grpMNprobs <- matrix(NA,1,Ngroups)
	classifprobs <- cbind(0, .1, .2, .3, .4, .5, .6, .7, .75, .80, .85, .90, .95, 1)
	grp_prob_Ns <- grp_prob_proports <- matrix(NA,ncol(classifprobs),Ngroups)

	for (lupec in 1:Ngroups) {	
		dumm <- subset(DFAposteriors, Group == grpnames[lupec], select=postvarnoms)

		# grpMNprobs[1,lupec] <- mean(dumm[,lupec])
		postcolnum <- which( substring(names(DFAposteriors), 11) == grpnames[lupec])
		grpMNprobs[1,lupec] <- mean(dumm[,postcolnum])
	
		for (luper in 1:nrow(grp_prob_Ns) ) 			
			grp_prob_Ns[luper,lupec] <- length(dumm[,lupec][dumm[,lupec] > classifprobs[1,luper]]) 
			# grp_prob_proports[luper,lupec] <- grp_prob_Ns[luper,lupec] / groupNs[lupec] 
	}	
	grpMNprobs <- matrix(grpMNprobs, 1, length(grpMNprobs))
	dimnames(grpMNprobs) <- list(rep("", dim(grpMNprobs)[1])); colnames(grpMNprobs) <- grpnames

	grp_prob_proports <- mapply('/', data.frame(grp_prob_Ns), groupNs)


	# # message("\n\nAverage of the highest classification probabilities: ", 
	# message("\n\nAverage of the (original) group classification probabilities: ", 
	    # round(mean(apply( DFAposteriors[,-grep("Group",colnames(DFAposteriors))], 1, max)),2)  
		
	# # the uncertainties of the classifications
	# uncertanties <- 1 - (apply(DFAposteriors[,-grep("Group",colnames(DFAposteriors))], 1, max))
	# message("\n\nAverage of the classification uncertainties: ", round(mean(uncertanties),2))
	
	grp_prob_Ns <- cbind( t(classifprobs), grp_prob_Ns)
	dimnames(grp_prob_Ns) <-list(rep("", dim(grp_prob_Ns)[1]))
	colnames(grp_prob_Ns) <- c("Classif. prob.", grpnames )
	
	grp_prob_proports <- cbind( t(classifprobs), grp_prob_proports)
	dimnames(grp_prob_proports) <-list(rep("", dim(grp_prob_proports)[1]))
	colnames(grp_prob_proports) <- c("Classif. prob.", grpnames )

	
	if (verbose) {
		message("\nGroup Sizes: \n"); print(groupNs)
			
		message("\n\nGroup mean posterior classification probabilities: \n"); 
		print(round(grpMNprobs,2), print.gap=4)
			
		message("\n\nNumber of cases per level of posterior classification probability:\n") 
		print(grp_prob_Ns, print.gap=4)
			
		message("\n\nProportions of cases per level of posterior classification probability:\n") 
		print(grp_prob_proports, print.gap=4); message("\n")		
	}	

	output <- list(groupNs=groupNs, grpMNprobs=grpMNprobs, grp_prob_Ns=grp_prob_Ns, grp_prob_proports=grp_prob_proports) 

	return(invisible(output))
} 








DFA_classes <- function(donnes, grpmeans, Ngroups, groupNs, grpnames, Ncases, W, priorprob='SIZES', covmat_type=covmat_type) {

	Ndvs <- ncol(donnes) - 1

	# the group prior probabilities
	if (priorprob == 'EQUAL')  prior <- matrix( (1 / Ngroups), 1, Ngroups)
	if (priorprob == 'SIZES')  prior <- matrix( (groupNs / Ncases), 1, Ngroups)
	names(prior) <- grpnames	


	if (covmat_type == 'within')  {   
		
		# # APPROACH 1 = T & F 2013, p 389 - 391

		# # the weights
		# Cj_TabFid <- t(apply(grpmeans, 1, function(grpmeans, W) { solve(W) %*% grpmeans }, W=W))
		# colnames(Cj_TabFid) <- colnames(donnes[,2:(Ndvs+1)])
	
		# # # the intercepts
		# # # T & F 2013 p 389 = -.5 * t(Cj_TabFid[x,]) %*% grpmeans[x,];   but SPSS = log(prior[x]) -.5 * t(Cj_TabFid[x,]) %*% grpmeans[x,]
		# # Cj0_TabFid <- sapply(1:Ngroups, function(x, grpmeans, Cj_TabFid) { log(prior[x]) -.5 * t(Cj_TabFid[x,]) %*% grpmeans[x,] }, grpmeans=grpmeans, Cj_TabFid=Cj_TabFid)
		# # names(Cj0_TabFid) <- grpnames

		# # to get the same intercepts as Rencher:
		# Cj0_TabFid <- sapply(1:Ngroups, function(x, grpmeans, Cj_TabFid) { -.5 * t(Cj_TabFid[x,]) %*% grpmeans[x,] - 1 }, grpmeans=grpmeans, Cj_TabFid=Cj_TabFid)
		# names(Cj0_TabFid) <- grpnames

		# # # Cj0_TabFid <- Cj0_Rencher ; 		names(Cj0_TabFid) <- grpnames
		
		# # the TabFid intercepts (i.e., computed with log(prior) are = the SPSS intercepts, but not = Rencher's, slight diff), 
		# # but the classes are NOT the same as the SPSS classes if the TabFid intercepts are used
		# # the classes are = the SPSS classes & the MASS classes when the Rencher intercepts are used
		# # my guess: SPSS uses log(prior) to display the intercepts, but not to compute the classes

		# # scores on the classification functions  -- T & F 2013 p 391 (altho they do not use log(prior) for equal grp sizes)
		# clsfxnvals <- matrix(NA, Ncases, Ngroups)
		# for (ng in 1:length(grpnames)) { # for each observation get the scores on the group classification functions
			# for (luper in 1:Ncases) {
				# clsfxnvals[luper,ng] <- Cj0_TabFid[grpnames[ng]] + sum(Cj_TabFid[grpnames[ng],] * donnes[luper,(2:(Ndvs+1))]) + log(prior[grpnames[ng]])  
			# }
		# }

		# dfa_class_TabFid <- c()
		# for (luper in 1:Ncases)  dfa_class_TabFid[luper] <- grpnames[which.max(clsfxnvals[luper,])]

 		# dfa_class_TabFid <- factor(dfa_class_TabFid, levels=grpnames)
			

		
		# APPROACH 2 = Linear Discriminant Analysis of Several Groups 
		# Rencher (2002, p. 304), but R code adpated from 
		# Schlegel  https://rpubs.com/aaronsc32/lda-classification-several-groups	
		# Schlegel - Linear Classification Functions for Several Groups.pdf
			
		# split the data into groups
		donnes.groups <- split(donnes[,2:ncol(donnes)], donnes[,1])
		
		# get the group mean vectors
		donnes.grpmeans <- lapply(donnes.groups, function(x) { c(apply(x, 2, mean)) })		

		# get the groups covariance matrices
		Si <- lapply(donnes.groups, function(x) cov(x))		
		
		# pooling the group covariance matrices
		summat <- matrix(0, Ndvs, Ndvs)
		for (Ngs in 1:Ngroups) { summat = summat + ( (groupNs[Ngs] - 1) * matrix(unlist(Si[Ngs]), Ndvs, byrow = TRUE) ) }
		sp1 <- (1 / (Ncases - Ngroups)) * summat  # the pooled covmat

		# ISSUE: SPSS uses Type III sums of squares for the var-covar matrices, Rencher does not 
		# print(sp1);	print(W)		
		sp1 <- W  # using Rencher's commands but with the SPSS vcv instead of Rencher's sp1

		sp1_inv <- solve(sp1)

		Li.y <- matrix(-9999, Ncases, Ngroups)
		dfa_class_Rencher <- c()    #matrix(-9999, Ncases, 1) 
				
		for (luper in 1:Ncases) {			  
			y <- as.numeric(matrix( donnes[luper,2:ncol(donnes)], Ndvs, 1) )			  
													  
			for (ng in 1:length(grpnames)) { # for each observation get the scores on the group classification functions

				y.bar <- matrix( unlist(donnes.grpmeans[grpnames[ng]]), Ndvs, 1)
					
				Li.y[luper,ng] <- log(prior[grpnames[ng]]) + t(y.bar) %*% sp1_inv %*% (y) - .5 * t(y.bar) %*% sp1_inv %*% (y.bar)
			}			  
			dfa_class_Rencher[luper] <- grpnames[which.max(Li.y[luper,])]   # the group number which maximizes the function
		}
		dfa_class_Rencher <- factor(dfa_class_Rencher, levels=grpnames)


		# Rencher Cj0 & Cj   -- p 305
		Cj0_Rencher <- c()
		Cj_Rencher  <- matrix(NA, Ngroups, Ndvs)
		for (ng in 1:length(grpnames)) { 
			y.bar <- matrix( unlist(donnes.grpmeans[grpnames[ng]]), Ndvs, 1)
			Cj0_Rencher[ng] <- -.5 * t(y.bar) %*% sp1_inv %*% (y.bar) - 1	# to get Schleger's intercepts   -1 correct?????					
			# # add log(prior[grpnames[ng]]) to get the SPSS intercepts
			# Cj0_Rencher[ng] <- log(prior[grpnames[ng]]) -.5 * t(y.bar) %*% sp1_inv %*% (y.bar)		   					
			Cj_Rencher[ng,] <-  t(y.bar) %*% sp1_inv 			
		}			  
		colnames(Cj_Rencher) <- colnames(donnes[,2:(Ndvs+1)])
		rownames(Cj_Rencher) <- grpnames
		names(Cj0_Rencher) <- grpnames



		# # APPROACH 3 = MASS::lda 
		# ldaoutput <- MASS::lda(x = as.matrix(donnes[,2:ncol(donnes)]), grouping=donnes[,1], prior = prior) 		
		# lda.values <- stats::predict(ldaoutput, donnes[,2:ncol(donnes)]) # obtain scores on the DFs
		# dfa_class_MASS <- lda.values$class
		# # Cj_MASS  <- 
		# # Cj0_MASS <- 
				

		# # comparing TabFid with Rencher   INTERCEPTS
		# message('\n\ncomparing TabFid with Rencher:\n')
		# print(Cj0_TabFid)
		# print(Cj0_Rencher)
		# print(round((Cj0_TabFid - Cj0_Rencher),3))

		# # comparing TabFid with Rencher   WEIGHTS
		# message('\n\ncomparing TabFid with Rencher:\n')
		# print(Cj_TabFid)
		# print(Cj_Rencher)
		# print(round((Cj_TabFid - Cj_Rencher),3))


		# # comparing ORIG with TabFid
		# message('\ncomparing ORIG with TabFid:\n')
		# print(squareTable(donnes[,1], dfa_class_TabFid, faclevels=grpnames, tabdimnames=c('variable 1','variable 2')))

		# # comparing ORIG with Rencher
		# message('\ncomparing ORIG with Rencher:\n')
		# print(squareTable(donnes[,1], dfa_class_Rencher, faclevels=grpnames, tabdimnames=c('variable 1','variable 2')))

		# # comparing ORIG with MASS
		# message('\ncomparing ORIG with MASS:\n')
		# print(squareTable(donnes[,1], dfa_class_MASS, faclevels=grpnames, tabdimnames=c('variable 1','variable 2')))


		# # comparing TabFid with Rencher
		# message('\ncomparing TabFid with Rencher:\n')
		# print(squareTable(dfa_class_TabFid, dfa_class_Rencher, faclevels=grpnames, tabdimnames=c('variable 1','variable 2')))

		# # comparing TabFid with MASS::lda
		# message('\ncomparing TabFid with MASS::lda:\n')
		# print(squareTable(dfa_class_TabFid, dfa_class_MASS, faclevels=grpnames, tabdimnames=c('variable 1','variable 2')))

		# # comparing Rencher with MASS::lda
		# message('\ncomparing Rencher with MASS::lda:\n')
		# print(squareTable(dfa_class_Rencher, dfa_class_MASS, faclevels=grpnames, tabdimnames=c('variable 1','variable 2')))


	
		# posterior probabilities -- adapted from 
		# 2019 Boedeker - LDA for Prediction of Group Membership - A User-Friendly Primer - Appendix A
		# the code below produces the same posteriors values as MASS::lda
		detW <- det(W)   # determinant of the pooled variance-covariance matrix
		invW <- solve(W)		
		fs_formula_pt1 <- (1/(((2*pi)^(Ndvs/2))*(detW^.5)))
		posteriors <- c()
		for (luper in 1:Ncases) {
			fs <- c()
			for (lupeg in 1:Ngroups) {
				fs   <- append(fs, fs_formula_pt1 * 
				                   exp(-.5*(t(t(donnes[luper,2:ncol(donnes)]) - grpmeans[lupeg,])) 
				                   %*% invW %*% (t(donnes[luper,2:ncol(donnes)]) - grpmeans[lupeg,])) )	
			}
			posteriors <- rbind(posteriors, ((fs * prior) / sum(fs * prior)) )
		}
		colnames(posteriors) <- paste("posterior_", grpnames, sep="")
		posteriors <- as.data.frame(posteriors)
		posteriors$Group <- donnes[,1]
		rownames(posteriors) <- 1:Ncases
	

		# producing a version of prior that does not display attr separately
		prior.print <- data.frame(prior)
		dimnames(prior.print)[[1]] <- '';  colnames(prior.print) <- grpnames	
		
		DFAclass_output <- list(dfa_class=dfa_class_Rencher, prior=prior.print, 
		                        classifcoefs=t(Cj_Rencher), classifints=Cj0_Rencher, posteriors=posteriors) 
	}


	if (covmat_type == 'separate') {

		# use of separate-groups covariance matrices for classification is called Quadratic Discriminant Analysis
		# SPSS has a separate-groups covariance matrices option, but it uses group cov matrices based on the DFs - not great
		# the MASS package has a qda function (which my code replicates) -- prob = no clear way to get/display the classif function coeffs
		
		# Quadratic Discriminant Analysis of Several Groups 
		# Rencher (2002, p. 306), but R code adpated from Schlegel  https://rpubs.com/aaronsc32/qda-several-groups		
		# split the data into groups & get the groups covariance matrices and group mean vectors
		donnes.groups <- split(donnes[,2:ncol(donnes)], donnes[,1])
		Si <- lapply(donnes.groups, function(x) cov(x))		
		donnes.grpmeans <- lapply(donnes.groups, function(x) { c(apply(x, 2, mean)) })		
		
		# # Rencher Cj0 & Cj   -- p 305
		# Cj0_Rencher <- c()
		# Cj_Rencher  <- matrix(NA, Ngroups, Ndvs)
	
		dfa_class_Rencher <- c()     

		for (luper in 1:Ncases) {			  
			y <- donnes[luper,2:ncol(donnes)] 			  
			l2i <- c()		  
			for (j in 1:Ngroups) { # For each group, calculate the QDA function 
				y.bar <- unlist(donnes.grpmeans[j])
				Si.j <- matrix(unlist(Si[j]), Ndvs, byrow = TRUE)
				Si.j_inv <- solve(Si.j)
				l2i <- append(l2i, -.5 * log(det(Si.j)) - .5 * as.numeric(y - y.bar) %*% Si.j_inv %*% as.numeric(y - y.bar) + log(prior[grpnames[j]]) ) 

				# Cj0_Rencher[j] <- -.5 * t(y.bar) %*% Si.j_inv %*% (y.bar)		               - 1	#  -1 correct?????					
				# Cj_Rencher[j,] <-  t(y.bar) %*% Si.j_inv 			
			}			  
			dfa_class_Rencher <- append(dfa_class_Rencher, grpnames[which.max(l2i)])   # the group which maximizes the function
		}

		# colnames(Cj_Rencher) <- colnames(donnes[,2:(Ndvs+1)])
		# rownames(Cj_Rencher) <- grpnames
		# names(Cj0_Rencher) <- grpnames

		dfa_class_Rencher <- factor(dfa_class_Rencher, levels=grpnames)


		# # comparing with MASS::qda
		# # library(MASS)
		# qda_MASS <- qda(x = as.matrix(donnes[,2:ncol(donnes)]), grouping=donnes[,1], prior=prior, CV=FALSE)
		# dfa_class_MASS <- predict(qda_MASS)$class 
		# # print(table(dfa_class_Rencher, dfa_class_MASS))
		# # 1 - sum(dfa_class_Rencher == dfa_class_MASS) / Ncases   # error rate
 		# message('\ncomparing Rencher with MASS::qda:\n')
		# print(squareTable(dfa_class_Rencher, dfa_class_MASS, faclevels=grpnames, tabdimnames=c('variable 1','variable 2')))



		# posterior probabilities -- adapted from 
		# 2019 Boedeker - LDA for Prediction of Group Membership - A User-Friendly Primer - Appendix A
		# the code below produces the same posteriors values as MASS::qda
		detWgrp <- c()
		invWgrp <- replicate(Ngroups, matrix(NA, Ndvs, Ndvs), simplify = F)
		for (lupeg in 1:Ngroups) {
			dum <- subset(donnes, donnes[,1] == grpnames[lupeg] )		
			Wgrp <- stats::cov(dum[,2:ncol(dum)])
			invWgrp[[lupeg]] <- solve(Wgrp)
			detWgrp <- append(detWgrp, det(Wgrp))   
		}			
		posteriors <- c()
		for (luper in 1:Ncases) {
			fs <- c()
			for (lupeg in 1:Ngroups) {
				fs   <- append(fs, (1/(((2*pi)^(Ndvs/2))*(detWgrp[lupeg]^.5))) * 
				                   exp(-.5*(t(t(donnes[luper,2:ncol(donnes)]) - grpmeans[lupeg,])) 
				                   %*% invWgrp[[lupeg]] %*% (t(donnes[luper,2:ncol(donnes)]) - grpmeans[lupeg,])) )	
			}
			posteriors <- rbind(posteriors,  ((fs * prior) / sum(fs * prior)) )
		}
		colnames(posteriors) <- paste("posterior_", grpnames, sep="")
		posteriors <- as.data.frame(posteriors)
		posteriors$Group <- donnes[,1]
		rownames(posteriors) <- 1:Ncases


		# producing a version of prior that does not display attr separately
		prior.print <- data.frame(prior)
		dimnames(prior.print)[[1]] <- '';  colnames(prior.print) <- grpnames	
		
		DFAclass_output <- list(dfa_class=dfa_class_Rencher, prior=prior.print, 
		                        classifcoefs=NULL, classifints=NULL, posteriors=posteriors)
	}

	return(invisible(DFAclass_output))
}




 





DFA_classes_CV <- function(donnes, priorprob='SIZES', covmat_type, grpnames) {

	# Cross-Validation of Predicted Groups
	# In cases with small sample sizes, prediction error rates can tend to 
	# be optimistic. Cross-validation is a technique used to estimate how accurate a predictive 
	# model may be in actual practice. When larger sample sizes are available, the more common 
	# approach of splitting the data into test and training sets may still be employed. There 
	# are many different approaches to cross-validation, including leave-p-out and k-fold 
	# cross-validation. One particular case of leave-p-out cross-validation is the leave-one-out 
	# approach, also known as the holdout method.
	
	# Leave-one-out cross-validation is performed by using all but one of the sample 
	# observation vectors to determine the classification function and then using that 
	# classification function to predict the omitted observations group membership. 
	# The procedure is repeated for each observation so that each is classified by a 
	# function of the other observations.

	Ncases <- nrow(donnes)
	Ncasesm1 <- Ncases - 1
	Ndvs <- ncol(donnes) - 1
	Ngroups  <- length(unique(donnes[,1]))
		
	if (covmat_type == 'within')  {   # classification -- T & F 2013, p 389 - 391

		clsfxnvals <- matrix(NA, Ncases, Ngroups)

		for (luper in 1:nrow(donnes)) {
			
			donnesm1 <- donnes[-luper,]
	
			MV2 <- manova( as.matrix(donnesm1[,2:ncol(donnesm1)]) ~ donnesm1[,1], data=donnesm1)
			sscpwith <- (Ncasesm1 - 1) * cov(MV2$residuals) # E	
			W <- sscpwith * (1 / (Ncasesm1 - Ngroups))   # vcv 
			
			# classification -- T & F 2013, p 389 - 391
	
			grpmeans_CV <- sapply(2:ncol(donnesm1), function(x) tapply(donnesm1[,x], INDEX = donnesm1[,1], FUN = mean))
	
			groupNs_CV <- as.matrix(table(donnesm1[,1]))
						
			# the weights
			Cj <- t(apply(grpmeans_CV, 1, function(grpmeans_CV, W) { solve(W) %*% grpmeans_CV }, W=W))
			colnames(Cj) <- colnames(donnes[,2:(Ndvs+1)])
		
			# the constants
			Cj0 <- sapply(1:Ngroups, function(x, grpmeans_CV, Cj) { -.5 * t(Cj[x,]) %*% grpmeans_CV[x,] }, grpmeans_CV=grpmeans_CV, Cj=Cj)
			names(Cj0) <- grpnames

			if (priorprob == 'EQUAL') {
				for (ng in 1:Ngroups) clsfxnvals[luper,ng] <- Cj0[grpnames[ng]] + sum(Cj[grpnames[ng],] * donnes[luper,(2:(Ndvs+1))])
			}
			if (priorprob == 'SIZES') {  
				prior <- matrix( (groupNs_CV / Ncasesm1), 1, Ngroups)
				names(prior) <- grpnames	
				for (ng in 1:Ngroups) 
					 clsfxnvals[luper,ng] <- Cj0[grpnames[ng]] + sum(Cj[grpnames[ng],] * donnes[luper,(2:(Ndvs+1))]) + log(prior[grpnames[ng]]) 
			}
		}				

		dfa_class_TabFid_CV <- c()
		for (luper in 1:Ncases)  dfa_class_TabFid_CV[luper] <- grpnames[which.max(clsfxnvals[luper,])]

 		dfa_class_TabFid_CV <- factor(dfa_class_TabFid_CV, levels=grpnames) 	

		dfa_class_CV <- dfa_class_TabFid_CV
	}


	if (covmat_type == 'separate') {

		dfa_class_Rencher_CV <- c()
		for (luper in 1:nrow(donnes)) {
			
			donnesm1 <- donnes[-luper,]
	
			# use of separate-groups covariance matrices for classification is called Quadratic Discriminant Analysis
			# SPSS has a separate-groups covariance matrices option, but it uses group cov matrices based on the DFs - not great
			# the MASS package has a qda function (which my code replicates) -- prob = no clear way to get/display the classif function coeffs
			
			# Quadratic Discriminant Analysis of Several Groups 
			# Rencher (2002, p. 306), but R code adpated from Schlegel  https://rpubs.com/aaronsc32/qda-several-groups
			
			# split the data into groups & get the groups covariance matrices and group mean vectors
			donnesm1.groups <- split(donnesm1[,2:ncol(donnesm1)], donnesm1[,1])
			Si <- lapply(donnesm1.groups, function(x) cov(x))		
			donnesm1.grpmeans_CV <- lapply(donnesm1.groups, function(x) { c(apply(x, 2, mean)) })
			groupNs_CV <- as.matrix(table(donnesm1[,1]))
			
			# the group prior probabilities
			if (priorprob == 'EQUAL')  prior <- matrix( (1 / Ngroups), 1, Ngroups)
			if (priorprob == 'SIZES')  prior <- matrix( (groupNs_CV / Ncasesm1), 1, Ngroups)
			names(prior) <- grpnames	
		
				y <- donnes[luper,2:ncol(donnesm1)] 
				  
				l2i <- c()		  
				for (j in 1:Ngroups) { # For each group, calculate the QDA function 
					y.bar <- unlist(donnesm1.grpmeans_CV[j])
					Si.j <- matrix(unlist(Si[j]), Ndvs, byrow = TRUE)
					Si.j_inv <- solve(Si.j)
					l2i <- append(l2i, -.5 * log(det(Si.j)) - .5 * as.numeric(y - y.bar) %*% Si.j_inv %*% as.numeric(y - y.bar) + log(prior[grpnames[j]]) )
				}
				  
				dfa_class_Rencher_CV <- append(dfa_class_Rencher_CV, grpnames[which.max(l2i)])   # the group which maximizes the function
		}
		dfa_class_Rencher_CV <- factor(dfa_class_Rencher_CV, levels=grpnames)


		# # comparing with MASS::qda
		# # library(MASS)
		# dfa_class_MASS_CV <- qda(x = as.matrix(donnes[,2:ncol(donnes)]), grouping=donnes[,1], prior=prior, CV=TRUE)$class
		# # print(table(dfa_class_Rencher_CV, dfa_class_MASS_CV))
		# # 1 - sum(dfa_class_Rencher_CV == dfa_class_MASS_CV) / Ncases   # error rate
 		# message('\nCV: comparing Rencher with MASS::qda:\n')
		# print(squareTable(dfa_class_Rencher_CV, dfa_class_MASS_CV, faclevels=grpnames, tabdimnames=c('variable 1','variable 2')))


		dfa_class_CV <- dfa_class_Rencher_CV
	}
	
return(invisible(dfa_class_CV))
}		
		








# rounds numeric columns in a matrix
# numeric columns named 'p' or 'plevel' or 'plevel.adj' are rounded to round_p places
# numeric columns not named 'p' are rounded to round_non_p places

round_boc <- function(donnes, round_non_p = 3, round_p = 5) {
	
	# identify the numeric columns
	#	numers <- apply(donnes, 2, is.numeric)  # does not work consistently 
	for (lupec in 1:ncol(donnes)) {

		if (is.numeric(donnes[,lupec]) == TRUE) 
		
			if (colnames(donnes)[lupec] == 'p' | colnames(donnes)[lupec] == 'plevel'  | 
			    colnames(donnes)[lupec] == 'p adj.')  {
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





squareTable <- function(var1, var2, faclevels=NULL, tabdimnames=c('variable 1','variable 2')) {

	# # # my original commands
	# # # prob = non-square contin table when the # of factor levels are not the same for var1 and var2
    # var1fact <- factor(var1, labels = grpnames)
    # var2fact <- factor(var2, labels = grpnames)
    # table(var1fact, var2fact, dnn=tabdimnames)

	# length(levels(var1fact)) == length(levels(var2fact))
	# var2fact <- factor(var2)  # var1 is already a factor, but var2 isn't



	if (!is.factor(var1)) { var1fact <- factor(var1, levels=faclevels) } else { var1fact <- var1; levels(var1fact) <- faclevels }
           
    var2fact <- factor(var2, levels=levels(var1fact) )   #, labels = sort(levels1[unique(var1)]) )
     
    contintab <- table(var1fact, var2fact, dnn=tabdimnames)


# # 	length(levels(var1fact)) == length(levels(var2fact))
	# var2fact <- factor(var2)  # var1 is already a factor, but var2 isn't


   # var2fact <- var2
        
		# attributes(var2fact) <- attributes(var1fact) 




	# if (is.factor(var1)) {    # var1 <- as.integer(var1)
	
		# levels1 <- levels(var1)
		
# #		var2 <- factor(var2, labels = sort(levels1[unique(var2)]) )

		# # R provides a simple and transparent way to apply one object's attributes to another object. 
		# # For factors, the attributes include the factor levels/labels.
		# # https://stackoverflow.com/questions/27453416/r-apply-one-factors-levels-to-another-object
		# attributes(var2) <- attributes(var1) 

		# contintab <- as.table(matrix(0, length(levels(var1)), length(levels(var1))))
		
		# rownames(contintab) <- colnames(contintab) <- levels(var1)

	# } else {

		# contintab <- as.table(matrix(0, length(unique(var1)), length(unique(var1))))
		
		# rownames(contintab) <- colnames(contintab) <- unique(var1)		
	# }
		
    
	# for (lupe in 1:length(var1)) {		
		# contintab[which(rownames(contintab) == var1[lupe]),  which(rownames(contintab) == var2[lupe])   ] =
		# contintab[which(rownames(contintab) == var1[lupe]),  which(rownames(contintab) == var2[lupe])   ] + 1
	# }

	# names(dimnames(contintab)) <- tabdimnames

    return(invisible(contintab))
}








kappa.cohen <- function (var1, var2, grpnames) {

	# the data for this function (kapdon) are the category values for each of 2 columns,
	# but the analyses are performed on the contingency table	

	kapdonCT <- squareTable(var1, var2, faclevels=grpnames); kapdonCT

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









# # sources for the following multivariate significance tests:

# Manly, B. F. J., & Alberto, J. A. (2017). Multivariate Statistical 
# Methods: A Primer (4th Edition). Chapman & Hall/CRC, Boca Raton, FL.

# https://rpubs.com/aaronsc32/manova-test-statistics

# http://www.statpower.net/Software.html

# Gatignon, H. (2014). Statistical Analysis of Management Data. New York, NY: Springer.




Wilks <- function(rho, Ncases, p, q) {
	
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
		
		if ( (p^2 * q^2 - 4) > 0)  { t <- sqrt((p^2 * q^2 - 4)/(p^2 + q^2 - 5))
		} else { t <- 1}
			
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
	RaoOutput <- cbind(round(wilks,2),round(Frao,2),df1,round(df2,2),round(pvalue,5))
	return(invisible(RaoOutput))
}
  
   
   
