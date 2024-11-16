
NORMALITY <- function(data, groups=NULL, variables=NULL, verbose = TRUE) {

# descriptive statistics & tests of univariate & multivariate normality -- from the MVN package
# uses the normality.R function; runs it for all variables and for variables within groups (optional)

#  if (!is.matrix(data) & !is.data.frame(data))  data <- as.matrix(data)


  

      
# no groups
if (is.null(groups)) {

	if (is.null(variables)) { 
	  
	  donnes <- as.matrix(data)
	
	  if (is.null(names(data))) colnames(donnes) <- 'data'
	}
  
	if (!is.null(variables)) {
	  
	  donnes <- as.matrix(data[,variables])
	  
	  colnames(donnes) <- variables
  }

	if (anyNA(donnes)) {
		donnes <- na.omit(donnes)
		message('\n\nCases with missing values were found and removed from the data matrix.\n\n')
	}

  if (!is.numeric(donnes))
    message('\nThe data are not numeric, expect errors.\n\n')
  
	normres <- umvn(donnes)

	NORMALITYoutput <- list(  
	   descriptives = normres$descriptives,
	   univariate_tests = normres$univariate_tests,
	   multivariate_tests = normres$multivariate_tests
	)

	if (verbose) {		
		message('\n\nDescriptive Statistics:\n')
		# print(round(normres$descriptives,3), print.gap=3)		
		print(t(round(normres$descriptives,3)), print.gap=3)
		message('\n\nTests of univariate normality:\n')
		print(normres$univariate_tests, row.names=FALSE, print.gap=4)
		if (ncol(donnes) > 1) {
		  message('\nTests of multivariate normality:\n')
		  print(round_boc(normres$multivariate_tests, round_non_p = 3, round_p = 5), 
		        row.names=FALSE, print.gap=4)
		}
	}
}

  
# when there are groups
if (!is.null(groups)) {
  
  if (is.null(variables)) 
    message('\nA "groups" variable was provided but without any "variables" to be analyzed. Expect errors.\n')
    
	donnes <- as.data.frame(data[,c(groups,variables)])

	if (anyNA(donnes)) {
		donnes <- na.omit(donnes)
		message('\n\nCases with missing values were found and removed from the data matrix.\n\n')
	}

	grpnames <- as.vector(as.matrix(donnes[,1])) # group names, in the same order as in data
	grpnames <- unique(grpnames)
	grpnums  <- seq(1:length(grpnames))
	ngroups  <- length(grpnames)
	
	# whole sample / all groups combined
	message('\n\nWhole-sample statistics')	
	normres <- umvn(donnes[,variables])

	NORMALITYoutput <- list(  
	   descriptives = normres$descriptives,
	   univariate_tests = normres$univariate_tests,
	   multivariate_tests = normres$multivariate_tests
	)

	if (verbose) {		
	  message('\n\nDescriptive Statistics:\n')
	  # print(round(normres$descriptives,3), print.gap=3)		
	  print(t(round(normres$descriptives,3)), print.gap=3)
	  message('\n\nTests of univariate normality:\n')
	  print(normres$univariate_tests, row.names=FALSE, print.gap=4)
	  if (ncol(donnes) > 1) {
	    message('\nTests of multivariate normality:\n')
	    print(round_boc(normres$multivariate_tests, round_non_p = 3, round_p = 5), 
	          row.names=FALSE, print.gap=4)
	  }
	}

	# separate stats for each group
	message('\n\n\nGroup statistics')
	for (lupeg in 1:ngroups) { 
		dum <- subset(donnes, donnes[,1] == grpnames[lupeg])
		message('\n\nGroup ', paste(grpnames[lupeg]))
		
		normres <- umvn(dum[,2:ncol(dum)])

		grpnom <- paste(grpnames[lupeg])
		
		NORMALITYoutput[[grpnom]]$descriptives        <- normres$descriptives
		NORMALITYoutput[[grpnom]]$univariate_tests    <- normres$univariate_tests
		NORMALITYoutput[[grpnom]]$multivariate_tests  <- normres$multivariate_tests

		if (verbose) {		
		  message('\n\nDescriptive Statistics:\n')
		  # print(round(normres$descriptives,3), print.gap=3)		
		  print(t(round(normres$descriptives,3)), print.gap=3)
		  message('\n\nTests of univariate normality:\n')
		  print(normres$univariate_tests, row.names=FALSE, print.gap=4)
		  if (ncol(donnes) > 1) {
		    message('\nTests of multivariate normality:\n')
		    print(round_boc(normres$multivariate_tests, round_non_p = 3, round_p = 5), 
		          row.names=FALSE, print.gap=4)
		  }
		}
	}
}	

return(invisible(NORMALITYoutput))

}

