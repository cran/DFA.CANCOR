\name{CANCOR}
\alias{CANCOR}
\title{Canonical correlation analysis}
\description{Produces SPSS- and SAS-like output for canonical 
	correlation analysis. Portions of the code were adapted from James Steiger (www.statpower.net).}
\usage{CANCOR(data, set1, set2, plot, plotCV, plotcoefs, verbose)}
\arguments{
  \item{data}{
  \code{}A dataframe where the rows are cases & the columns are the variables.}

  \item{set1}{
  \code{}The names of the continuous variables for the first set,
   \cr e.g., set1 = c('varA', 'varB', 'varC').}
  	    
  \item{set2}{
  \code{}The names of the continuous variables for the second set,
   \cr e.g., set2 = c('varD', 'varE', 'varF').}

  \item{plot}{
  \code{}Should a plot of the coefficients be produced? \cr The options are: TRUE (default) or FALSE.}

  \item{plotCV}{
  \code{}The canonical variate number for the plot, e.g., plotCV = 1.}

  \item{plotcoefs}{
  \code{}The coefficient for the plots. \cr The options are 'structure' (default) or 'standardized'.}

  \item{verbose}{
  \code{}Should detailed results be displayed in the console? \cr The options are: TRUE (default) or FALSE.}
}
\value{If verbose = TRUE, the displayed output includes Pearson correlations, multivariate significance
	   tests, canonical function correlations and bivariate significance tests, raw canonical 
	   coefficients, structure coefficients, standardized coefficients, and a bar plot  
	   of the structure or standardized coefficients.
	   
 	   The returned output is a list with elements
  \item{cancorrels}{canonical correlations and their significance tests}
  \item{CoefRawSet1}{raw canonical coefficients for Set 1}
  \item{CoefRawSet2}{raw canonical coefficients for Set 2}
  \item{CoefStruct11}{structure coefficients for Set 1 variables with the Set 1 variates}
  \item{CoefStruct21}{structure coefficients for Set 2 variables with the Set 1 variates}
  \item{CoefStruct12}{structure coefficients for Set 1 variables with the Set 2 variates}
  \item{CoefStruct22}{structure coefficients for Set 2 variables with the Set 2 variates}
  \item{CoefStandSet1}{standardized coefficients for Set 1 variables}
  \item{CoefStandSet2}{standardized coefficients for Set 2 variables}
  \item{mv_Wilk}{Wilk's multivariate significance test}
  \item{mv_Pillai}{Pillai-Bartlett multivariate significance test}
  \item{mv_Hotelling}{Hotelling-Lawley multivariate significance test}
  \item{mv_Roy}{Roy's Largest Root multivariate significance test}
  \item{mv_BartlettV}{Bartlett's V multivariate significance test}
  \item{mv_Rao}{Rao's' multivariate significance test} 
  \item{CorrelSet1}{Pearson correlations for Set 1}
  \item{CorrelSet2}{Pearson correlations for Set 2}
  \item{CorrelSet1n2}{Pearson correlations between Set 1 & Set 2}
}
\references{
	   {Manly, B. F. J., & Alberto, J. A. (2017). \emph{Multivariate statistical methods: 
	   	A primer (4th Edition).} Chapman & Hall/CRC, Boca Raton, FL.}
\cr\cr {Rencher, A. C. (2002). \emph{Methods of Multivariate Analysis} (2nd ed.). New York, NY: John Wiley & Sons.}
\cr\cr {Sherry, A., & Henson, R. K. (2005). Conducting and interpreting canonical correlation analysis
        in personality research: A user-friendly primer. \emph{Journal of Personality Assessment, 84,} 37-48.}
\cr\cr {Steiger, J. (2019). \emph{Canonical correlation analysis.}
	    \cr{www.statpower.net/Content/312/Lecture\%20Slides/CanonicalCorrelation.pdf}}
\cr\cr {Tabachnik, B. G., & Fidell, L. S. (2019). \emph{Using multivariate statistics (7th ed.).} New York, NY: Pearson.}
}
\author{Brian P. O'Connor }
\examples{
# data that simulate those from De Leo & Wulfert (2013)
CANCOR(data = na.omit(data_CCA_De_Leo),
       set1 = c('Tobacco_Use','Alcohol_Use','Illicit_Drug_Use','Gambling_Behavior',
                'Unprotected_Sex','CIAS_Total'),
       set2 = c('Impulsivity','Social_Interaction_Anxiety','Depression',
                'Social_Support','Intolerance_of_Deviance','Family_Morals',
                'Family_Conflict','Grade_Point_Average'),
       plot = TRUE, plotCV = 1, plotcoefs='structure',
       verbose = TRUE)

\donttest{
# data from Tabachnik & Fidell (2013, p. 589)
CANCOR(data = data_CCA_Tabachnik,
       set1 = c('TS','TC'),
       set2 = c('BS','BC'),
       plot = TRUE, plotCV = 1, plotcoefs='structure',
       verbose = TRUE)

# UCLA dataset
UCLA_CCA_data <- read.csv("https://stats.idre.ucla.edu/stat/data/mmreg.csv")
colnames(UCLA_CCA_data) <- c("LocusControl", "SelfConcept", "Motivation",
                             "read", "write", "math", "science", "female")
summary(UCLA_CCA_data)
CANCOR(data = UCLA_CCA_data, 
       set1 = c("LocusControl","SelfConcept","Motivation"),
       set2 = c("read","write","math","science","female"),
       plot = TRUE, plotCV = 1, plotcoefs='standardized',
       verbose = TRUE)
}

}
