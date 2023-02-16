\name{NORMALITY}
\alias{NORMALITY}
\title{Univariate and multivariate normality}
\description{Produces tests of univariate and multivariate normality using the MVN package.}
\usage{NORMALITY(data, groups, variables, verbose)}
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
	{Korkmaz, S., Goksuluk, D., Zararsiz, G. (2014). MVN: An R package for assessing 
	 multivariate normality. \emph{The R Journal, 6(2),} 151-162.}	
\cr\cr {Szekely,G. J., & Rizzo, M. L. (2017). The energy of data. \emph{Annual 
	 Review of Statistics and Its Application 4,} 447-79.}	 
\cr\cr {Tabachnik, B. G., & Fidell, L. S. (2019). \emph{Using multivariate statistics (7th ed.).} New York, NY: Pearson.}
}
\author{Brian P. O'Connor }
\examples{
# data that simulate those from De Leo & Wulfert (2013)
NORMALITY(data = na.omit(data_CANCOR$DeLeo_2013[c(
          'Unprotected_Sex','Tobacco_Use','Alcohol_Use','Illicit_Drug_Use',
          'Gambling_Behavior','CIAS_Total','Impulsivity','Social_Interaction_Anxiety',
          'Depression','Social_Support','Intolerance_of_Deviance','Family_Morals',
          'Family_Conflict','Grade_Point_Average')]))
 
\donttest{
# data from Field et al. (2012)
NORMALITY(data = data_DFA$Field_2012, 
          groups = 'Group', 
          variables = c('Actions','Thoughts'))

# data from Tabachnik & Fidell (2013, p. 589)
NORMALITY(data = na.omit(data_CANCOR$TabFid_2019_small[c('TS','TC','BS','BC')]))

# UCLA dataset
UCLA_CCA_data <- read.csv("https://stats.idre.ucla.edu/stat/data/mmreg.csv")
colnames(UCLA_CCA_data) <- c("LocusControl", "SelfConcept", "Motivation",
                             "read", "write", "math", "science", "female")
summary(UCLA_CCA_data)
NORMALITY(data = na.omit(UCLA_CCA_data[c("LocusControl","SelfConcept","Motivation",
                                         "read","write","math","science")]))
}

}
