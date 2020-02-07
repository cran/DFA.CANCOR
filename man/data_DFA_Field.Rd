\name{data_DFA_Field}
\alias{data_DFA_Field}
\docType{data}
\title{data_DFA_Field}
\description{A data frame with scores on 2 variables for 10 cases in each of 3 
	groups. Used by Field et al. (2012) in their chapter on MANOVA and 
	discriminant function analysis.}
\usage{data(data_DFA_Field)}

\source{
Field, A., Miles, J., & Field, Z. (2012). \emph{Discovering statistics using R.} Los Angeles, CA: Sage. }

\examples{
\donttest{
head(data_DFA_Field)
	 
DFA(data = data_DFA_Field, 
    groups = 'Group', 
    variables = c('Actions','Thoughts'),
    predictive = TRUE, priorprob = 'SIZES', 
    verbose = TRUE)
}
}