
LINEARITY <- function(data, variables=NULL, groups=NULL, idvs=NULL, dv=NULL, verbose = TRUE) {

# uses lm for the linear and quadratic associations between 2 continous variables


# when there are no groups
if (is.null(groups) == TRUE) {

	# idvs & dv are not specified
	if (is.null(idvs) == TRUE & is.null(dv) == TRUE) {
	
		if (is.null(variables) == TRUE)  donnes <- as.data.frame(data)
		
		if (is.null(variables) == FALSE) donnes <- as.data.frame(data[,variables])
	
		if (anyNA(donnes) == TRUE) {
			donnes <- na.omit(donnes)
			message('\n\nCases with missing values were found and removed from the data matrix.\n')
		}
	
		varnames <- colnames(donnes)
		LINEARITYOutput <- list()
		
		if (verbose == TRUE) message('\n\nTests for linear and quadratic effects')
	
		for (lupevars1 in 1:(length(varnames)-1)) {
			for (lupevars2 in (lupevars1 + 1):(length(varnames))) {
	
				iv <- donnes[,lupevars1]
				ivname <- varnames[lupevars1]
				dv <- donnes[,lupevars2]
				dvname <- varnames[lupevars2]
	
				linquadOut <- linquad(iv,dv)
				
				lname <- paste(paste(varnames[lupevars1]),'_(X)_&_', 
				               paste(varnames[lupevars2]),'_(Y))',sep='')
				LINEARITYOutput[[lname]]$coefs <- linquadOut
						
				if (verbose == TRUE) {		
					linquadOut <- cbind(round(linquadOut[,1:4],2),round(linquadOut[,5],6))
					colnames(linquadOut)[5] <- 'Pr(>|t|)'
					message('\n',paste(paste(varnames[lupevars1]),' (idv) & ', 
					           paste(varnames[lupevars2]),' (dv)',sep=''),'\n')
					print(round(linquadOut,5))		
				}
			}		
		}
	}


	# when idvs & dv are specified
	if (is.null(idvs) == FALSE & is.null(dv) == FALSE) {
	
		# idvs <- as.data.frame(data[,idvs])
		# idvnames <- colnames(donnes)
		
		if (verbose == TRUE) message('\n\nTests for linear and quadratic effects:\n')
	
		LINEARITYOutput <- list()
		
		for (lupevars in 1:length(idvs)) {

			dontemp <- data.frame(data[,idvs[lupevars]], data[,dv])
			
			if (anyNA(dontemp) == TRUE) {
				dontemp <- na.omit(dontemp)
				message('\n\nCases with missing values were found and removed from the data matrix.\n')
			}
								
			linquadOut <- linquad(dontemp[,1],dontemp[,2])
			
			lname <- paste(paste(idvs[lupevars]),'_(X)_&_', paste(dv),'_(Y))',sep='')
#			LINEARITYOutput[[lname]]$coefs <- linquadOut
					
			if (verbose == TRUE) {		
				linquadOut <- cbind(round(linquadOut[,1:4],2),round(linquadOut[,5],6))
				colnames(linquadOut)[5] <- 'Pr(>|t|)'
				message('\n',paste(paste(idvs[lupevars]),' (idv) & ', paste(dv),' (dv)',sep=''),'\n')
				print(round(linquadOut,5))		
			}		
		}
	}
}




# when there are groups
if (is.null(groups) == FALSE) {

	grpnames <- as.vector(as.matrix(data[,groups])) # group names, in the same order as in the data matrix
	grpnames <- unique(grpnames)
	grpnums  <- seq(1:length(grpnames))
	Ngroups  <- length(grpnums)

	# when idvs & dv are not specified
	if (is.null(idvs) == TRUE & is.null(dv) == TRUE) {
	
		if (is.null(variables) == TRUE)  {
			variables <- colnames(data)		
			variables[ variables != groups]	 
		}
		
		LINEARITYOutput <- list()
		
		if (verbose == TRUE) message('\n\nTests for linear and quadratic effects:')
	
		for (lupevars1 in 1:(length(variables)-1)) {
			for (lupevars2 in (lupevars1 + 1):(length(variables))) {
				for (lupeg in 0:Ngroups) {
					if (lupeg == 0) {					
						# whole sample / all groups combined
						message('\n\n\nWhole-sample statistics\n')		 										
						dontemp <- data.frame(data[,variables[lupevars1]], data[,variables[lupevars2]])			
						if (anyNA(dontemp) == TRUE) {
							dontemp <- na.omit(dontemp)
							message('\n\nCases with missing values were found and removed from the data matrix.\n')
						}											
						linquadOut <- linquad(dontemp[,1],dontemp[,2])
						lname <- paste(paste(variables[lupevars1]),'_(X)_&_', 
				                       paste(variables[lupevars2]),'_(Y))',sep='')
						LINEARITYOutput[[lname]]$coefs <- linquadOut
					}else{						
						dontemp <- data.frame(data[,groups],data[,variables[lupevars1]], data[,variables[lupevars2]])			
						dontemp <- subset(dontemp, dontemp[,1] == grpnames[lupeg])						
						message('\n\n\nGroup', paste(grpnames[lupeg]),'\n')
						linquadOut <- linquad(dontemp[,2],dontemp[,3])
						lname <- paste(paste(variables[lupevars1]),'_(X)_&_', 
				                       paste(variables[lupevars2]),'_(Y))',sep='')
						grpnom <- paste(grpnames[lupeg])
						LINEARITYOutput[[lname]][[grpnom]]$coefs <- linquadOut
					}						
					if (verbose == TRUE) {		
						linquadOut <- cbind(round(linquadOut[,1:4],2),round(linquadOut[,5],6))
						colnames(linquadOut)[5] <- 'Pr(>|t|)'
						message('\n',paste(paste(variables[lupevars1]),' (idv) & ', 
						           paste(variables[lupevars2]),' (dv)',sep=''),'\n')
						print(round(linquadOut,5))		
					}
				}
			}		
		}
	}

	# when idvs & dv are specified
	if (is.null(idvs) == FALSE & is.null(dv) == FALSE) {
	
		if (is.null(variables) == TRUE)  {
			variables <- colnames(data)		
			variables[ variables != groups]	 
		}
		
		LINEARITYOutput <- list()
		
		if (verbose == TRUE) message('\n\nTests for linear and quadratic effects')
	
		for (lupevars in 1:length(idvs)) {
			for (lupeg in 0:Ngroups) {
				if (lupeg == 0) {					
					# whole sample / all groups combined
					message('\n\n\nWhole-sample statistics\n')		 										
					dontemp <- data.frame(data[,variables[lupevars]], data[,dv])			
					if (anyNA(dontemp) == TRUE) {
						dontemp <- na.omit(dontemp)
						message('\n\nCases with missing values were found and removed from the data matrix.\n')
					}											
					linquadOut <- linquad(dontemp[,1],dontemp[,2])
					lname <- paste(paste(idvs[lupevars]),'_(X)_&_', paste(dv),'_(Y))',sep='')
					LINEARITYOutput[[lname]]$coefs <- linquadOut
				}else{						
					dontemp <- data.frame(data[,groups],data[,idvs[lupevars]], data[,dv])			
					dontemp <- subset(dontemp, dontemp[,1] == grpnames[lupeg])						
					message('\n\nGroup', paste(grpnames[lupeg]),'\n')
					linquadOut <- linquad(dontemp[,2],dontemp[,3])
					lname <- paste(paste(idvs[lupevars]),'_(X)_&_', paste(dv),'_(Y))',sep='')
					grpnom <- paste(grpnames[lupeg])
					LINEARITYOutput[[lname]][[grpnom]]$coefs <- linquadOut
				}						
				if (verbose == TRUE) {		
					linquadOut <- cbind(round(linquadOut[,1:4],2),round(linquadOut[,5],6))
					colnames(linquadOut)[5] <- 'Pr(>|t|)'
					message('\n',paste(paste(idvs[lupevars]),' (idv) & ', paste(dv),' (dv)',sep=''),'\n')
					print(round(linquadOut,5))		
				}
			}		
		}
	}
	
}

return(invisible(LINEARITYOutput))

}

