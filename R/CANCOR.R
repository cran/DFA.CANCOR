

CANCOR <- function(data, set1, set2, plot=TRUE, plotCV=1, plotcoefs='structure', verbose=TRUE ) {

data <- as.data.frame(data[,c(set1,set2)])
 
if (anyNA(data) == TRUE) {
	data <- na.omit(data)
#	message('\n\nCases with missing values were found and removed from the data matrix.\n')
	NAflag = TRUE
} else {
	NAflag = FALSE
}

set1data <- as.data.frame(data[,set1])
set2data <- as.data.frame(data[,set2])

Ncases <- nrow(set1data)
NVset1 <- ncol(set1data)
NVset2 <- ncol(set2data)

# Pearson correlations
CorrelSet1   <- stats::cor(set1data)
CorrelSet2   <- stats::cor(set2data)
CorrelSet1n2 <- stats::cor(set2data,set1data)

# the CCA
output <- canonical.cor(set1data, set2data) # source("http://www.statpower.net/R312/CanCorr.r")

cancorrels <- output$cancorrels

mv_Wilks <- Wilks(rho= cancorrels[,2], Ncases=Ncases, p = NVset1, q = NVset2)
colnames(mv_Wilks) <- c('Wilk\'s Lambda', 'F-approx. ', 'df1', 'df2', 'p')
rownames(mv_Wilks) <- paste(1:nrow(mv_Wilks), paste("through ", nrow(mv_Wilks), sep = ""))
# print(round(mv_Wilks,4)); message('\n')

mv_Pillai <- Pillai(rho= cancorrels[,2], Ncases=Ncases, p = NVset1, q = NVset2)
colnames(mv_Pillai) <- c('Pillai-Bartlett Trace', 'F-approx. ', 'df1', 'df2', 'p')
rownames(mv_Pillai) <- paste(1:nrow(mv_Pillai), paste("through ", nrow(mv_Pillai), sep = ""))
# print(round(mv_Pillai,4)); message('\n')

mv_Hotelling <- Hotelling(rho= cancorrels[,2], Ncases=Ncases, p = NVset1, q = NVset2)
colnames(mv_Hotelling) <- c('Hotelling-Lawley Trace', 'F-approx. ', 'df1', 'df2', 'p')
rownames(mv_Hotelling) <- paste(1:nrow(mv_Hotelling), paste("through ", nrow(mv_Hotelling), sep = ""))
# print(round(mv_Hotelling,4)); message('\n')

mv_Roy <- RoyRoot(rho = cancorrels[,2], Ncases=Ncases, p = NVset1, q = NVset2)
colnames(mv_Roy) <- c('Roy\'s Largest Root', 'lambda ', 'F-approx. ', 'df1', 'df2', 'p')
rownames(mv_Roy) <- paste(1:nrow(mv_Roy), paste("through ", nrow(mv_Roy), sep = ""))
# print(round(mv_Roy,4)); message('\n')

mv_BartlettV <- BartlettV(rho= cancorrels[,2], Ncases=Ncases, p = NVset1, q = NVset2)
# message('\n\nBartlett\'s V test:\n')
colnames(mv_BartlettV) <- c('Wilk\'s Lambda', 'F-approx. ', 'df', 'p')
rownames(mv_BartlettV) <- paste(1:nrow(mv_BartlettV), paste("through ", nrow(mv_BartlettV), sep = ""))
# print(round(mv_BartlettV,4)); message('\n')

mv_Rao <- Rao(rho= cancorrels[,2], Ncases=Ncases, p = NVset1, q = NVset2)
# message('\n\nRao\'s V test:\n')
colnames(mv_Rao) <- c('Wilk\'s Lambda', 'F-approx. ', 'df1', 'df2', 'p')
rownames(mv_Rao) <- paste(1:nrow(mv_Rao), paste("through ", nrow(mv_Rao), sep = ""))
# print(round(mv_Rao,4)); message('\n')



# bivariate correlations for the canonical variates
# message('\n\nCanonical correlations:\n')
colnames(output$cancorrels) <- c('Eigenvalue', 'Canonical r', 'Canonical r sq.','t','df','p value')
rownames(output$cancorrels) <- paste(" Canonical function ", 1:nrow(output$cancorrels),sep = "")
# print(round(output$cancorrels,3))
# message('\nThe above t tests are for the Pearson correlations between the canonical variate scores')
# message('\nfor each function, i.e., they are not the multivariate significance tests.\n')

  
# raw canonical coefficients
# message('\n\nRaw canonical coefficients for Set 1:\n')
colnames(output$raw1) <- paste('CV', 1:ncol(output$raw1),sep = "")
# print(round(output$raw1,2))

# message('\n\nRaw canonical coefficients for Set 2:\n')
colnames(output$raw2) <- paste('CV', 1:ncol(output$raw2),sep = "")
# print(round(output$raw2,2))


# structure coefficients (canonical loadings)
# message('\n\nStructure coefficients (Pearson correlations) for Set 1 variables with the Set 1 variates:\n')
colnames(output$struct11) <- paste('CV', 1:ncol(output$struct11),sep = "")
# print(round(output$struct11,2))

# message('\n\nStructure coefficients (Pearson correlations) for Set 2 variables with the Set 1 variates:\n')
colnames(output$struct21) <- paste('CV', 1:ncol(output$struct21),sep = "")
# print(round(output$struct21,2))

# message('\n\nStructure coefficients (Pearson correlations) for Set 1 variables with the Set 2 variates:\n')
colnames(output$struct12) <- paste('CV', 1:ncol(output$struct12),sep = "")
# print(round(output$struct12,2))

# message('\n\nStructure coefficients (Pearson correlations) for Set 2 variables with the Set 2 variates:\n')
colnames(output$struct22) <- paste('CV', 1:ncol(output$struct22),sep = "")
# print(round(output$struct22,2))


# standardized canonical coefficients 
# message('\n\nStandardized coefficients for Set 1 variables:\n')
colnames(output$stand1) <- paste('CV', 1:ncol(output$stand1),sep = "")
# print(round(output$stand1,2))

# message('\n\nStandardized coefficients for Set 2 variables:\n')
colnames(output$stand2) <- paste('CV', 1:ncol(output$stand2),sep = "")
# print(round(output$stand2,2))


# scores on the canonical functions
# standardize the variables
Z.set1 <- scale(set1data, center=TRUE, scale=FALSE)
Z.set2 <- scale(set2data, center=TRUE, scale=FALSE)
set1_scores <- Z.set1 %*% output$raw1
set2_scores <- Z.set2 %*% output$raw2



if (plot == TRUE | is.null(plot)) {

	oldpar <- par(no.readonly = TRUE)
	on.exit(par(oldpar))

	if (is.null(plot)) plotCV = 1

	# bar plots - structure coefficients
	if (plotcoefs == 'structure'| is.null(plotcoefs)) {
		layout(matrix(c(1, 2, 1, 2), nrow=2, byrow=T));  # layout.show(n=2)		
		barplot(output$struct11[,plotCV], ylim=c(-1,1), ylab='CV1', 
		        main='Set 1 Structure Coefficients', col='blue', las=2)
		box()	
		barplot(output$struct22[,plotCV], ylim=c(-1,1), ylab='CV1', 
		        main='Set 2 Structure Coefficients', col='blue', las=2)
		box()

		# par(mfrow=c(1,2), mar = c(9,5,3,3)) # set the margin on all sides
		# barplot(output$struct11[,plotCV], ylim=c(-1,1), ylab='CV1', 
		        # main='Set 1 Structure Coefficients', col='blue', las=2)
		# box()	
		# barplot(output$struct22[,plotCV], ylim=c(-1,1), ylab='CV1', 
		        # main='Set 2 Structure Coefficients', col='blue', las=2)
		# box()
	}
		
	# bar plots - standardized coefficients
	if (plotcoefs == 'standardized') {
		layout(matrix(c(1, 2, 1, 2), nrow=2, byrow=T));  # layout.show(n=2)		
		barplot(output$stand1[,plotCV], ylim=c(-1.3,1.3), ylab='CV1', 
		        main='Set 1 Standardized Coefficients', col='blue', las=2)
		box()	
		barplot(output$stand2[,plotCV], ylim=c(-1.3,1.3), ylab='CV1', 
		        main='Set 2 Standardized Coefficients', col='blue', las=2)
		box()

		# par(mfrow=c(1,2), mar = c(9,5,3,3)) # set the margin on all sides
		# barplot(output$stand1[,plotCV], ylim=c(-1.3,1.3), ylab='CV1', 
		        # main='Set 1 Standardized Coefficients', col='blue', las=2)
		# box()	
		# barplot(output$stand2[,plotCV], ylim=c(-1.3,1.3), ylab='CV1', 
		        # main='Set 2 Standardized Coefficients', col='blue', las=2)
		# box()
	}
	#dev.off()  

	# # helio plot -- from the yacca package
	# cca.fit <- yacca::cca(set1data, set2data)
	# yacca::helio.plot(cca.fit, x.name="Set 1", y.name="Set 2", cv=plotCV)
	# boc.fit <- list( xstructcorr=struct11, ystructcorr=struct22, xstructcorrsq=struct11**2,
                     # ystructcorrsq=struct22**2, xlab=rownames(struct11), ylab=rownames(struct22) )
	# yacca::helio.plot(boc.fit, x.name="Set 1", y.name="Set 2", cv=plotCV)
	# message('\n\nThe plot is provided by the yacca package. Helio plots display data in')
	# message('\nradial bars, with larger values pointing outward from a base reference circle')
	# message('\nand smaller (more negative) values pointing inward.')
}



if (verbose) {
	
	message('\n\nCanonical Correlation Analysis')
	 
	if (NAflag) message('\n\nCases with missing values were found and removed from the data matrix.\n')
		
	message('\n\nPearson correlations for Set 1:\n'); print(round(CorrelSet1,2), print.gap=4)
	message('\n\nPearson correlations for Set 2:\n'); print(round(CorrelSet2,2), print.gap=4)
	message('\n\nPearson correlations between Set 1 & Set 2:\n'); print(round(CorrelSet1n2,2), print.gap=4)
	
	message('\n\nMultivariate peel-down significance tests:\n')
	
	print(round(mv_Wilks,4), print.gap=4); message('\n')
	
	print(round(mv_Pillai,4), print.gap=4); message('\n')
	
	print(round(mv_Hotelling,4), print.gap=4); message('\n')
	
	print(round(mv_Roy,4), print.gap=4); message('\n')
	
	message('\nBartlett\'s V test:\n')
	print(round(mv_BartlettV,4), print.gap=4); message('\n')
	
	message('\nRao\'s V test:\n')
	print(round(mv_Rao,4), print.gap=4); message('\n')
		
	# bivariate correlations for the canonical variates
	message('\nCanonical correlations:\n')
	print(round(output$cancorrels,3), print.gap=4)
	message('\nThe above t tests are for the Pearson correlations between the canonical variate scores')
	message('for each function, i.e., they are not the multivariate significance tests.\n')
	
	# raw canonical coefficients
	message('\n\nRaw canonical coefficients for Set 1:\n')
	print(round(output$raw1,2), print.gap=4)
	
	message('\n\nRaw canonical coefficients for Set 2:\n')
	print(round(output$raw2,2), print.gap=4)
	
	# structure coefficients (canonical loadings)
	message('\n\nStructure coefficients (Pearson correlations) for Set 1 variables with the Set 1 variates:\n')
	print(round(output$struct11,2), print.gap=4)
	
	message('\n\nStructure coefficients (Pearson correlations) for Set 2 variables with the Set 1 variates:\n')
	print(round(output$struct21,2), print.gap=4)
	
	message('\n\nStructure coefficients (Pearson correlations) for Set 1 variables with the Set 2 variates:\n')
	print(round(output$struct12,2), print.gap=4)
	
	message('\n\nStructure coefficients (Pearson correlations) for Set 2 variables with the Set 2 variates:\n')
	print(round(output$struct22,2), print.gap=4)
	
	# standardized canonical coefficients 
	message('\n\nStandardized coefficients for Set 1 variables:\n')
	print(round(output$stand1,2), print.gap=4)
	
	message('\n\nStandardized coefficients for Set 2 variables:\n')
	print(round(output$stand2,2), print.gap=4)
	
	message('\n\n')
}


CANCORoutput <- list(  
   cancorrels = cancorrels,
   mv_Wilks = mv_Wilks,
   mv_Pillai = mv_Pillai,
   mv_Hotelling = mv_Hotelling,
   mv_Roy = mv_Roy,
   mv_BartlettV = mv_BartlettV,
   mv_Rao = mv_Rao,
   CoefRawSet1 = output$raw1,
   CoefRawSet2 = output$raw2,
   CoefStruct11 = output$struct11,
   CoefStruct21 = output$struct21,
   CoefStruct12 = output$struct12,
   CoefStruct22 = output$struct22,
   CoefStandSet1 = output$stand1,
   CoefStandSet2 = output$stand2,   
   CorrelSet1 = CorrelSet1,
   CorrelSet2 = CorrelSet2,
   CorrelSet1n2 = CorrelSet1n2,
   set1_scores = set1_scores,
   set2_scores = set2_scores
)

return(invisible(CANCORoutput))

class(CANCORoutput) <- "CANCORoutput"

}



