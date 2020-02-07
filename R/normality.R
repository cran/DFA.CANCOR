
normality <- function(data, groups=NULL, variables=NULL, verbose = TRUE) {

# descriptive statistics & tests of univariate & multivariate normality -- from the MVN package
# uses the normality.R function; runs it for all variables and for variables within groups (optional)

# no groups
if (is.null(groups) == TRUE) {

	if (is.null(variables) == TRUE)  donnes <- as.data.frame(data)
	
	if (is.null(variables) == FALSE) donnes <- as.data.frame(data[,variables])

	if (anyNA(donnes) == TRUE) {
		donnes <- na.omit(donnes)
		cat('\n\nCases with missing values were found and removed from the data matrix.\n\n')
	}

	normres <- umvn(donnes)

	normalityoutput <- list(  
	   descriptives = normres$descriptives,
	   Shapiro_Wilk = normres$Shapiro_Wilk,
	   Mardia = normres$Mardia,
	   Henze_Zirkler = normres$Henze_Zirkler,
	   Royston = normres$Royston,
	   Doornik_Hansen = normres$Doornik_Hansen
	)

	if (verbose == TRUE) {		
		cat('\n\nDescriptive Statistics:\n\n')
		print(round(normres$descriptives,3))		
		cat('\n\n\nShapiro-Wilk tests of univariate normality:\n\n')
		print(normres$Shapiro_Wilk)		
		cat('\n\n\nTests of multivariate normality:\n\n')
		print(normres$Mardia, row.names = FALSE); cat('\n\n')
		print(normres$Henze_Zirkler, row.names = FALSE); cat('\n\n')
		print(normres$Royston, row.names = FALSE); cat('\n\n')
		print(normres$Doornik_Hansen, row.names = FALSE); cat('\n\n')
	}
}

# when there are groups
if (is.null(groups) == FALSE) {
	
	donnes <- as.data.frame(data[,c(groups,variables)])

	if (anyNA(donnes) == TRUE) {
		donnes <- na.omit(donnes)
		cat('\n\nCases with missing values were found and removed from the data matrix.\n\n')
	}

	grpnames <- as.vector(as.matrix(donnes[,1])) # group names, in the same order as in the data matrix
	grpnames <- unique(grpnames)
	grpnums  <- seq(1:length(grpnames))
	ngroups  <- length(grpnames)
	
	# whole sample / all groups combined
	cat('\n\nWhole-sample statistics\n')	
	normres <- umvn(donnes[,2:ncol(donnes)])

	normalityoutput <- list(  
	   descriptives = normres$descriptives,
	   Shapiro_Wilk = normres$Shapiro_Wilk,
	   Mardia = normres$Mardia,
	   Henze_Zirkler = normres$Henze_Zirkler,
	   Royston = normres$Royston,
	   Doornik_Hansen = normres$Doornik_Hansen
	)

	if (verbose == TRUE) {
		cat('\n\nDescriptive Statistics:\n\n')
		print(round(normres$descriptives,3))
		cat('\n\n\nShapiro-Wilk tests of univariate normality:\n\n')
		print(normres$Shapiro_Wilk)		
		cat('\n\n\nTests of multivariate normality:\n\n')
		print(normres$Mardia, row.names = FALSE); cat('\n\n')
		print(normres$Henze_Zirkler, row.names = FALSE); cat('\n\n')
		print(normres$Royston, row.names = FALSE); cat('\n\n')
		print(normres$Doornik_Hansen, row.names = FALSE); cat('\n\n')
	}
	
	# separate stats for each group
	cat('\n\nGroup statistics\n\n')
	for (lupeg in 1:ngroups) { 
		dum <- subset(donnes, donnes[,1] == grpnames[lupeg])
		cat('\n\nGroup', paste(grpnames[lupeg]),'\n')
		normres <- umvn(dum[,2:ncol(dum)])

		grpnom <- paste(grpnames[lupeg])
		
		normalityoutput[[grpnom]]$descriptives   <- normres$descriptives
		normalityoutput[[grpnom]]$Shapiro_Wilk   <- normres$Shapiro_Wilk
		normalityoutput[[grpnom]]$Mardia         <- normres$Mardia
		normalityoutput[[grpnom]]$Henze_Zirkler  <- normres$Henze_Zirkler
		normalityoutput[[grpnom]]$Royston        <- normres$Royston
		normalityoutput[[grpnom]]$Doornik_Hansen <- normres$Doornik_Hansen

		if (verbose == TRUE) {
			cat('\n\nDescriptive Statistics:\n\n')
			print(round(normres$descriptives,3))
			cat('\n\n\nShapiro-Wilk tests of univariate normality:\n\n')
			print(normres$Shapiro_Wilk)		
			cat('\n\n\nTests of multivariate normality:\n\n')
			print(normres$Mardia, row.names = FALSE); cat('\n\n')
			print(normres$Henze_Zirkler, row.names = FALSE); cat('\n\n')
			print(normres$Royston, row.names = FALSE); cat('\n\n')
			print(normres$Doornik_Hansen, row.names = FALSE); cat('\n\n')
		}
	}
}	

return(invisible(normalityoutput))

}

