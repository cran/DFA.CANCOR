

CANCOR <- function(data, set1, set2, plot=TRUE, plotCV=1, plotcoefs='structure', verbose=TRUE ) {

data <- as.data.frame(data[,c(set1,set2)])
 
if (anyNA(data) == TRUE) {
	data <- na.omit(data)
#	cat('\n\nCases with missing values were found and removed from the data matrix.\n\n')
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

mv_Wilk <- Wilk(rho= cancorrels[,2], Ncases=Ncases, p = NVset1, q = NVset2)
colnames(mv_Wilk) <- c('Wilk\'s Lambda', 'F-approx. ', 'df1', 'df2', 'p')
rownames(mv_Wilk) <- paste(1:nrow(mv_Wilk), paste("through ", nrow(mv_Wilk), sep = ""))
# print(round(mv_Wilk,4)); cat('\n\n')

mv_Pillai <- Pillai(rho= cancorrels[,2], Ncases=Ncases, p = NVset1, q = NVset2)
colnames(mv_Pillai) <- c('Pillai-Bartlett Trace', 'F-approx. ', 'df1', 'df2', 'p')
rownames(mv_Pillai) <- paste(1:nrow(mv_Pillai), paste("through ", nrow(mv_Pillai), sep = ""))
# print(round(mv_Pillai,4)); cat('\n\n')

mv_Hotelling <- Hotelling(rho= cancorrels[,2], Ncases=Ncases, p = NVset1, q = NVset2)
colnames(mv_Hotelling) <- c('Hotelling-Lawley Trace', 'F-approx. ', 'df1', 'df2', 'p')
rownames(mv_Hotelling) <- paste(1:nrow(mv_Hotelling), paste("through ", nrow(mv_Hotelling), sep = ""))
# print(round(mv_Hotelling,4)); cat('\n\n')

mv_Roy <- RoyRoot(rho = cancorrels[,2], Ncases=Ncases, p = NVset1, q = NVset2)
colnames(mv_Roy) <- c('Roy\'s Largest Root', 'lambda ', 'F-approx. ', 'df1', 'df2', 'p')
rownames(mv_Roy) <- paste(1:nrow(mv_Roy), paste("through ", nrow(mv_Roy), sep = ""))
# print(round(mv_Roy,4)); cat('\n\n')

mv_BartlettV <- BartlettV(rho= cancorrels[,2], Ncases=Ncases, p = NVset1, q = NVset2)
# cat('\n\n\nBartlett\'s V test:\n')
colnames(mv_BartlettV) <- c('Wilk\'s Lambda', 'F-approx. ', 'df', 'p')
rownames(mv_BartlettV) <- paste(1:nrow(mv_BartlettV), paste("through ", nrow(mv_BartlettV), sep = ""))
# print(round(mv_BartlettV,4)); cat('\n\n')

mv_Rao <- Rao(rho= cancorrels[,2], Ncases=Ncases, p = NVset1, q = NVset2)
# cat('\n\n\nRao\'s V test:\n')
colnames(mv_Rao) <- c('Wilk\'s Lambda', 'F-approx. ', 'df1', 'df2', 'p')
rownames(mv_Rao) <- paste(1:nrow(mv_Rao), paste("through ", nrow(mv_Rao), sep = ""))
# print(round(mv_Rao,4)); cat('\n\n')



# bivariate correlations for the canonical variates
# cat('\n\n\nCanonical correlations:\n\n')
colnames(output$cancorrels) <- c('Eigenvalue', 'Canonical r', 'Canonical r sq.','t','df','p value')
rownames(output$cancorrels) <- paste(" Canonical function ", 1:nrow(output$cancorrels),sep = "")
# print(round(output$cancorrels,3))
# cat('\nThe above t tests are for the Pearson correlations between the canonical variate scores')
# cat('\nfor each function, i.e., they are not the multivariate significance tests.\n\n')

  
# raw canonical coefficients
# cat('\n\n\n\nRaw canonical coefficients for Set 1:\n\n')
colnames(output$raw1) <- paste('CV', 1:ncol(output$raw1),sep = "")
# print(round(output$raw1,2))

# cat('\n\nRaw canonical coefficients for Set 2:\n\n')
colnames(output$raw2) <- paste('CV', 1:ncol(output$raw2),sep = "")
# print(round(output$raw2,2))


# structure coefficients (canonical loadings)
# cat('\n\n\nStructure coefficients (Pearson correlations) for Set 1 variables with the Set 1 variates:\n\n')
colnames(output$struct11) <- paste('CV', 1:ncol(output$struct11),sep = "")
# print(round(output$struct11,2))

# cat('\n\nStructure coefficients (Pearson correlations) for Set 2 variables with the Set 1 variates:\n\n')
colnames(output$struct21) <- paste('CV', 1:ncol(output$struct21),sep = "")
# print(round(output$struct21,2))

# cat('\n\nStructure coefficients (Pearson correlations) for Set 1 variables with the Set 2 variates:\n\n')
colnames(output$struct12) <- paste('CV', 1:ncol(output$struct12),sep = "")
# print(round(output$struct12,2))

# cat('\n\nStructure coefficients (Pearson correlations) for Set 2 variables with the Set 2 variates:\n\n')
colnames(output$struct22) <- paste('CV', 1:ncol(output$struct22),sep = "")
# print(round(output$struct22,2))


# standardized canonical coefficients 
# cat('\n\n\nStandardized coefficients for Set 1 variables:\n\n')
colnames(output$stand1) <- paste('CV', 1:ncol(output$stand1),sep = "")
# print(round(output$stand1,2))

# cat('\n\nStandardized coefficients for Set 2 variables:\n\n')
colnames(output$stand2) <- paste('CV', 1:ncol(output$stand2),sep = "")
# print(round(output$stand2,2))



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
	# cat('\n\n\nThe plot is provided by the yacca package. Helio plots display data in')
	# cat('\nradial bars, with larger values pointing outward from a base reference circle')
	# cat('\nand smaller (more negative) values pointing inward.')
}



if (verbose) {
	
	cat('\n\n\nCanonical Correlation Analysis\n')
	 
	if (NAflag) cat('\n\nCases with missing values were found and removed from the data matrix.\n\n')
		
	cat('\n\nPearson correlations for Set 1:\n\n'); print(round(CorrelSet1,2), print.gap=4)
	cat('\n\nPearson correlations for Set 2:\n\n'); print(round(CorrelSet2,2), print.gap=4)
	cat('\n\nPearson correlations between Set 1 & Set 2:\n\n'); print(round(CorrelSet1n2,2), print.gap=4)
	
	cat('\n\n\nMultivariate peel-down significance tests:\n\n')
	
	print(round(mv_Wilk,4), print.gap=4); cat('\n\n')
	
	print(round(mv_Pillai,4), print.gap=4); cat('\n\n')
	
	print(round(mv_Hotelling,4), print.gap=4); cat('\n\n')
	
	print(round(mv_Roy,4), print.gap=4); cat('\n\n')
	
	cat('\n\n\nBartlett\'s V test:\n')
	print(round(mv_BartlettV,4), print.gap=4); cat('\n\n')
	
	cat('\n\n\nRao\'s V test:\n')
	print(round(mv_Rao,4), print.gap=4); cat('\n\n')
		
	# bivariate correlations for the canonical variates
	cat('\n\n\nCanonical correlations:\n\n')
	print(round(output$cancorrels,3), print.gap=4)
	cat('\nThe above t tests are for the Pearson correlations between the canonical variate scores')
	cat('\nfor each function, i.e., they are not the multivariate significance tests.\n\n')
	
	# raw canonical coefficients
	cat('\n\n\n\nRaw canonical coefficients for Set 1:\n\n')
	print(round(output$raw1,2), print.gap=4)
	
	cat('\n\nRaw canonical coefficients for Set 2:\n\n')
	print(round(output$raw2,2), print.gap=4)
	
	# structure coefficients (canonical loadings)
	cat('\n\n\nStructure coefficients (Pearson correlations) for Set 1 variables with the Set 1 variates:\n\n')
	print(round(output$struct11,2), print.gap=4)
	
	cat('\n\nStructure coefficients (Pearson correlations) for Set 2 variables with the Set 1 variates:\n\n')
	print(round(output$struct21,2), print.gap=4)
	
	cat('\n\nStructure coefficients (Pearson correlations) for Set 1 variables with the Set 2 variates:\n\n')
	print(round(output$struct12,2), print.gap=4)
	
	cat('\n\nStructure coefficients (Pearson correlations) for Set 2 variables with the Set 2 variates:\n\n')
	print(round(output$struct22,2), print.gap=4)
	
	# standardized canonical coefficients 
	cat('\n\n\nStandardized coefficients for Set 1 variables:\n\n')
	print(round(output$stand1,2), print.gap=4)
	
	cat('\n\nStandardized coefficients for Set 2 variables:\n\n')
	print(round(output$stand2,2), print.gap=4)
	
	cat('\n\n\n\n')
}


CANCORoutput <- list(  
   cancorrels = cancorrels,
   CoefRawSet1 = output$raw1,
   CoefRawSet2 = output$raw2,
   CoefStruct11 = output$struct11,
   CoefStruct21 = output$struct21,
   CoefStruct12 = output$struct12,
   CoefStruct22 = output$struct22,
   CoefStandSet1 = output$stand1,
   CoefStandSet2 = output$stand2,   
   mv_Wilk = mv_Wilk,
   mv_Pillai = mv_Pillai,
   mv_Hotelling = mv_Hotelling,
   mv_Roy = mv_Roy,
   mv_BartlettV = mv_BartlettV,
   mv_Rao = mv_Rao,
   CorrelSet1 = CorrelSet1,
   CorrelSet2 = CorrelSet2,
   CorrelSet1n2 = CorrelSet1n2
)

#return(invisible(CANCORoutput))

class(CANCORoutput) <- "CANCORout"

}



