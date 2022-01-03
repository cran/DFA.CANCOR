
NORMALITY <- function(data, groups=NULL, variables=NULL, verbose = TRUE) {

# descriptive statistics & tests of univariate & multivariate normality -- from the MVN package
# uses the normality.R function; runs it for all variables and for variables within groups (optional)

# no groups
if (is.null(groups) == TRUE) {

	if (is.null(variables) == TRUE)  donnes <- as.data.frame(data)
	
	if (is.null(variables) == FALSE) donnes <- as.data.frame(data[,variables])

	if (anyNA(donnes) == TRUE) {
		donnes <- na.omit(donnes)
		message('\n\nCases with missing values were found and removed from the data matrix.\n\n')
	}

	normres <- umvn(donnes)

	NORMALITYoutput <- list(  
	   descriptives = normres$descriptives,
	   Shapiro_Wilk = normres$Shapiro_Wilk,
	   Mardia = normres$Mardia,
	   Henze_Zirkler = normres$Henze_Zirkler,
	   Royston = normres$Royston,
	   Doornik_Hansen = normres$Doornik_Hansen
	)

	if (verbose == TRUE) {		
		message('\n\nDescriptive Statistics:\n\n')
		print(round(normres$descriptives,3))		
		message('\n\nShapiro-Wilk tests of univariate normality:\n\n')
		print(normres$Shapiro_Wilk)		
		message('\n\nTests of multivariate normality:\n\n')
		print(normres$Mardia, row.names = FALSE); message('\n\n')
		print(normres$Henze_Zirkler, row.names = FALSE); message('\n\n')
		print(normres$Royston, row.names = FALSE); message('\n\n')
		print(normres$Doornik_Hansen, row.names = FALSE); message('\n\n')
	}
}

# when there are groups
if (is.null(groups) == FALSE) {
	
	donnes <- as.data.frame(data[,c(groups,variables)])

	if (anyNA(donnes) == TRUE) {
		donnes <- na.omit(donnes)
		message('\n\nCases with missing values were found and removed from the data matrix.\n\n')
	}

	grpnames <- as.vector(as.matrix(donnes[,1])) # group names, in the same order as in the data matrix
	grpnames <- unique(grpnames)
	grpnums  <- seq(1:length(grpnames))
	ngroups  <- length(grpnames)
	
	# whole sample / all groups combined
	message('\n\nWhole-sample statistics')	
	normres <- umvn(donnes[,2:ncol(donnes)])

	NORMALITYoutput <- list(  
	   descriptives = normres$descriptives,
	   Shapiro_Wilk = normres$Shapiro_Wilk,
	   Mardia = normres$Mardia,
	   Henze_Zirkler = normres$Henze_Zirkler,
	   Royston = normres$Royston,
	   Doornik_Hansen = normres$Doornik_Hansen
	)

	if (verbose == TRUE) {
		message('\n\nDescriptive Statistics:\n')
		print(round(normres$descriptives,3))
		message('\n\nShapiro-Wilk tests of univariate normality:\n')
		print(normres$Shapiro_Wilk)		
		message('\n\nTests of multivariate normality:\n')
		print(normres$Mardia, row.names = FALSE); message('\n')
		print(normres$Henze_Zirkler, row.names = FALSE); message('\n')
		print(normres$Royston, row.names = FALSE); message('\n')
		print(normres$Doornik_Hansen, row.names = FALSE)
	}
	
	# separate stats for each group
	message('\n\n\nGroup statistics')
	for (lupeg in 1:ngroups) { 
		dum <- subset(donnes, donnes[,1] == grpnames[lupeg])
		message('\n\nGroup ', paste(grpnames[lupeg]))
		normres <- umvn(dum[,2:ncol(dum)])

		grpnom <- paste(grpnames[lupeg])
		
		NORMALITYoutput[[grpnom]]$descriptives   <- normres$descriptives
		NORMALITYoutput[[grpnom]]$Shapiro_Wilk   <- normres$Shapiro_Wilk
		NORMALITYoutput[[grpnom]]$Mardia         <- normres$Mardia
		NORMALITYoutput[[grpnom]]$Henze_Zirkler  <- normres$Henze_Zirkler
		NORMALITYoutput[[grpnom]]$Royston        <- normres$Royston
		NORMALITYoutput[[grpnom]]$Doornik_Hansen <- normres$Doornik_Hansen

		if (verbose == TRUE) {
			message('\n\nDescriptive Statistics:\n')
			print(round(normres$descriptives,3))
			message('\n\nShapiro-Wilk tests of univariate normality:\n')
			print(normres$Shapiro_Wilk)		
			message('\n\nTests of multivariate normality:\n')
			print(normres$Mardia, row.names = FALSE); message('\n')
			print(normres$Henze_Zirkler, row.names = FALSE); message('\n')
			print(normres$Royston, row.names = FALSE); message('\n')
			print(normres$Doornik_Hansen, row.names = FALSE); message('\n')
		}
	}
}	

return(invisible(NORMALITYoutput))

}

