\name{NORMALITY}
\alias{NORMALITY}
\title{Univariate and multivariate normality}
\description{Produces tests of univariate and multivariate normality using the MVN package.}
\usage{NORMALITY(data, groups, variables, verbose)}
\arguments{
  \item{data}{
  \code{}A dataframe or numeric matrix where the rows are cases & the columns are the variables.}

  \item{groups}{
  \code{}(optional) The name of the groups variable in the dataframe, \cr e.g., groups = 'Group'.}
  	    
  \item{variables}{ 
  \code{}(optional) The names of the continuous variables in the dataframe for the analyses, 
         e.g., variables = c('varA', 'varB', 'varC').}
  	    
  \item{verbose}{
  \code{}Should detailed results be displayed in the console? \cr The options are: TRUE (default) or FALSE.}
}
\details{If "groups" is not specified, the analyses will be run on all of the variables in "data".
     If "variables" is specified, the analyses will be run on the "variables" in "data".
	   If "groups" is specified, the analyses will be run for every value of "groups".
	   If verbose = TRUE, the displayed output includes descriptive statistics and
	   tests of univariate and multivariate normality.
}
\value{The returned output is a list with the following elements:
  \item{descriptives}{descriptive statistics, including skewness and kurtosis}
  \item{univariate_tests}{the univariate normality tests}
  \item{multivariate_tests}{the multivariate normality tests}
}
\references{
  {Doornik, J. A. & Hansen, H. (2008). An Omnibus test for univariate and 
   multivariate normality. \emph{Oxford Bulletin of Economics and Statistics 70,} 927-939.}
\cr\cr {Henze, N., & Wagner, T. (1997), A new approach to the BHEP tests for multivariate 
   normality. \emph{Journal of Multivariate Analysis, 62,} 1-23.}
\cr\cr {Johnson, R. A., & Wichern, D. W. (2007). \emph{Applied Multivariate 
   Statistical Analysis (3rd. ed.).} New Jersey, NJ: Prentice Hall.}
\cr\cr {Korkmaz, S., Goksuluk, D., Zararsiz, G. (2014). MVN: An R package for assessing 
	 multivariate normality. \emph{The R Journal, 6(2),} 151-162.}	
\cr\cr {Mardia, K. V. (1970), Measures of multivariate skewnees and kurtosis with 
   applications. \emph{Biometrika, 57(3),} 519-530.}
\cr\cr {Mardia, K. V. (1974), Applications of some measures of multivariate skewness 
   and kurtosis for testing normality and robustness studies. \emph{Sankhy A, 36,} 115-128.}
\cr\cr {Royston, J. P. (1992). Approximating the Shapiro-Wilk W-Test for non-normality. 
   \emph{Statistics and Computing, 2,} 117-119.}
\cr\cr {Shapiro, S., & Wilk, M. (1965). An analysis of variance test for normality. 
   \emph{Biometrika, 52,} 591611.}
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
