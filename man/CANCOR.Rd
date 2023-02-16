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
 \item{mv_Wilk}{The Wilk's lambda multivariate test}
  \item{mv_Pillai}{The Pillai-Bartlett multivariate test}
  \item{mv_Hotelling}{The Lawley-Hotelling multivariate test}
  \item{mv_Roy}{Roy's greatest characteristic root multivariate test}
  \item{mv_BartlettV}{Bartlett's V multivariate significance test}
  \item{mv_Rao}{Rao's' multivariate significance test} 
  \item{CoefRawSet1}{raw canonical coefficients for Set 1}
  \item{CoefRawSet2}{raw canonical coefficients for Set 2}
  \item{CoefStruct11}{structure coefficients for Set 1 variables with the Set 1 variates}
  \item{CoefStruct21}{structure coefficients for Set 2 variables with the Set 1 variates}
  \item{CoefStruct12}{structure coefficients for Set 1 variables with the Set 2 variates}
  \item{CoefStruct22}{structure coefficients for Set 2 variables with the Set 2 variates}
  \item{CoefStandSet1}{standardized coefficients for Set 1 variables}
  \item{CoefStandSet2}{standardized coefficients for Set 2 variables}
   \item{CorrelSet1}{Pearson correlations for Set 1}
  \item{CorrelSet2}{Pearson correlations for Set 2}
  \item{CorrelSet1n2}{Pearson correlations between Set 1 & Set 2}
  \item{set1_scores}{Canonical variate scores for Set 1}
  \item{set2_scores}{Canonical variate scores for Set 2}
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
CANCOR(data = data_CANCOR$DeLeo_2013,
       set1 = c('Tobacco_Use','Alcohol_Use','Illicit_Drug_Use','Gambling_Behavior',
                'Unprotected_Sex','CIAS_Total'),
       set2 = c('Impulsivity','Social_Interaction_Anxiety','Depression',
                'Social_Support','Intolerance_of_Deviance','Family_Morals',
                'Family_Conflict','Grade_Point_Average'),
       plot = TRUE, plotCV = 1, plotcoefs='structure',
       verbose = TRUE)

\donttest{
	
# data from Ho (2014, Chapter 17)
CANCOR(data = data_CANCOR$Ho_2014,
       set1 = c("willing_use","likely_use","intend_use","certain_use"),
       set2 = c("perceived_risk","perceived_severity","self_efficacy",
                "response_efficacy","maladaptive_coping","fear"),
       plot = 'yes', plotCV = 1)

	
# data from Rencher (2002, pp. 366, 369, 372)
CANCOR(data = data_CANCOR$Rencher_2002,
       set1 = c("y1","y2","y3"),
       set2 = c("x1","x2","x3","x1x2","x1x3","x2x3","x1sq","x2sq","x3sq"),
       plot = 'yes', plotCV = 1)


# data from Tabachnik & Fidell (2019, p. 451, 460)    small dataset
CANCOR(data = data_CANCOR$TabFid_2019_small,
       set1 = c('TS','TC'),
       set2 = c('BS','BC'),
       plot = TRUE, plotCV = 1, plotcoefs='structure',
       verbose = TRUE)


# data from Tabachnik & Fidell (2019, p. 463)     complete dataset
CANCOR(data = data_CANCOR$TabFid_2019_complete,
       set1 = c("esteem","control","attmar","attrole"),
       set2 = c("timedrs","attdrug","phyheal","menheal","druguse"),
       plot = TRUE, plotCV = 1, plotcoefs='structure',
       verbose = TRUE)


# UCLA dataset   https://stats.oarc.ucla.edu/r/dae/canonical-correlation-analysis/
CANCOR(data = data_CANCOR$UCLA, 
       set1 = c("Locus_Control","Self_Concept","Motivation"),
       set2 = c("Read","Write","Math","Science","Sex"),
       plot = TRUE, plotCV = 1, plotcoefs='standardized',
       verbose = TRUE)
}

}
