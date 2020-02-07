\name{normality}
\alias{normality}
\title{Univariate and multivariate normality}
\description{P tests of univariate and multivariate normality using the MVN package.}
\usage{normality(data, groups, variables, verbose)}
\arguments{
  \item{data}{
  \code{}A dataframe where the rows are cases & the columns are the variables.}

  \item{groups}{
  \code{}(optional) The name of the groups variable in the dataframe, \cr e.g., groups = 'Group'.}
  	    
  \item{variables}{ 
  \code{}(optional) The names of the continuous variables in the dataframe for the analyses, 
         e.g., variables = c('varA', 'varB', 'varC').}
  	    
  \item{verbose}{
  \code{}Should detailed results be displayed in the console? \cr The options are: TRUE (default) or FALSE.}
}
\value{If "groups" is not specified, the analyses will be run on all of the variables in "data".
	   If "groups" is specified, the analyses will be run for every value of "groups".
	   If "variables" is specified, the analyses will be run on the "variables" in "data".
	   If verbose = TRUE, the displayed output includes descriptive statistics and
	   tests of univariate and multivariate normality.
	   	   
	   The returned output is a list with elements
  \item{descriptives}{descriptive statistics, including skewness and kurtosis}
  \item{Shapiro_Wilk}{the Shapiro_Wilk test of univariate normality}
  \item{Mardia}{the Mardia test of multivariate normality}
  \item{Henze_Zirkler}{the Henze-Zirkler test of multivariate normality}
  \item{Royston}{the Royston test of multivariate normality}
  \item{Doornik_Hansen}{the Doornik_Hansen test of multivariate normality}
}
\references{
	{Korkmaz, S., Goksuluk, D., Zararsiz, G. (2014). MVN: An R Package for Assessing 
	 Multivariate Normality. \emph{The R Journal, 6(2),} 151-162.}
	
	{Szekely,G. J., & Rizzo, M. L. (2017). The Energy of Data. \emph{The Annual 
	 Review of Statistics and Its Application 4,} 447-79.}
	 
	{Tabachnik, B. G., & Fidell, L. S. (2013). \emph{Using multivariate statistics 
	 (6th. ed., pp. 79-83).} New York, NY: Pearson.}

}
\author{Brian P. O'Connor }
\examples{
# data that simulate those from De Leo & Wulfert (2013)
normality(data = na.omit(data_CCA_De_Leo[c(
          'Unprotected_Sex','Tobacco_Use','Alcohol_Use','Illicit_Drug_Use',
          'Gambling_Behavior','CIAS_Total','Impulsivity','Social_Interaction_Anxiety',
          'Depression','Social_Support','Intolerance_of_Deviance','Family_Morals',
          'Family_Conflict','Grade_Point_Average')]))
 
\donttest{
# data from Field et al. (2012)
normality(data = data_DFA_Field, 
          groups = 'Group', 
          variables = c('Actions','Thoughts'))

# data from Tabachnik & Fidell (2013, p. 589)
normality(data = na.omit(data_CCA_Tabachnik[c('TS','TC','BS','BC')]))

# UCLA dataset
UCLA_CCA_data <- read.csv("https://stats.idre.ucla.edu/stat/data/mmreg.csv")
colnames(UCLA_CCA_data) <- c("LocusControl", "SelfConcept", "Motivation",
                             "read", "write", "math", "science", "female")
summary(UCLA_CCA_data)
normality(data = na.omit(UCLA_CCA_data[c("LocusControl","SelfConcept","Motivation",
                                         "read","write","math","science","female")]))
}

}
