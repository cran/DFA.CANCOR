\name{LINEARITY}
\alias{LINEARITY}
\title{Linearity}
\description{Provides tests of the possible linear and quadratic
	         associations between two continuous variables.}
\usage{LINEARITY(data, variables, groups, idvs, dv, verbose)}
\arguments{
 \item{data}{
  \code{}A dataframe where the rows are cases & the columns are the variables.}

 \item{variables}{ 
  \code{}(optional) The names of the continuous variables in the dataframe for the analyses, 
         e.g., variables = c('varA', 'varB', 'varC').}
  	    
  \item{groups}{
  \code{}(optional) The name of the groups variable in the dataframe (if there is one), 
  		 \cr e.g., groups = 'Group'.}
  	    
  \item{idvs}{
  \code{}(optional) The names of the predictor variables, 
  		  e.g., variables = c('varA', 'varB', 'varC').}
  	    
  \item{dv}{
  \code{}(optional) The name of the dependent variable, if output for just 
  		  one dependent variable is desired.}
  	    
  \item{verbose}{
  \code{}(optional) Should detailed results be displayed in the console? 
  		\cr The options are: TRUE (default) or FALSE.}
}
\value{If "variables" is specified, the analyses will be run on the "variables" in "data".
	   If "groups" is specified, the analyses will be run for every value of "groups".
	   If verbose = TRUE, the linear and quadratic regression
	   coefficients and their statistical tests are displayed.
	   	   
	   The returned output is a list with the regression 
	   coefficients and their statistical tests.}
\references{
Tabachnik, B. G., & Fidell, L. S. (2019). \emph{Using multivariate statistics (7th ed.).} New York, NY: Pearson.
}
\author{Brian P. O'Connor }
\examples{
# data from Sherry (2006), using all variables
LINEARITY(data=data_DFA_Sherry, groups='Group',
          variables=c('Neuroticism','Extroversion','Openness',
                      'Agreeableness','Conscientiousness') )

\donttest{
# data from Sherry (2006), specifying independent variables and a dependent variable
LINEARITY(data=data_DFA_Sherry, groups='Group',
          idvs=c('Neuroticism','Extroversion','Openness','Agreeableness'),
          dv=c('Conscientiousness'), 
          verbose=TRUE )

# data that simulate those from De Leo & Wulfert (2013)
LINEARITY(data=data_CCA_De_Leo,
          variables=c('Tobacco_Use','Alcohol_Use','Illicit_Drug_Use',
                      'Gambling_Behavior', 'Unprotected_Sex','CIAS_Total', 
                      'Impulsivity','Social_Interaction_Anxiety','Depression',
                      'Social_Support','Intolerance_of_Deviance','Family_Morals',
                      'Family_Conflict','Grade_Point_Average'), 
          verbose=TRUE )
}

}
