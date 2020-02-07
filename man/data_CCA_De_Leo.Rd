\name{data_CCA_De_Leo}
\alias{data_CCA_De_Leo}
\docType{data}
\title{data_CCA_De_Leo}
\description{A data frame with scores on 14 variables that have the same correlational 
	structure, and which produce the same canonical correlation analysis results, 
	as those reported in De Leo and Wulfert (2013).}
\usage{data(data_CCA_De_Leo)}

\source{De Leo, J. A., & Wulfert, E. (2013). Problematic internet use and other risky 
	behaviors in college students: An application of problem-behavior theory. \emph{Psychology 
	of Addictive Behaviors, 27(1),} 133-141. }

\examples{ 
\donttest{
head(data_CCA_De_Leo)
     
CANCOR(data = na.omit(data_CCA_De_Leo),
       set1 = c('Tobacco_Use','Alcohol_Use','Illicit_Drug_Use','Gambling_Behavior',
                'Unprotected_Sex','CIAS_Total'),
       set2 = c('Impulsivity','Social_Interaction_Anxiety','Depression',
                'Social_Support','Intolerance_of_Deviance','Family_Morals',
                'Family_Conflict','Grade_Point_Average'),
       plot = 'yes', plotCV = 1,
       verbose=TRUE)
}
}