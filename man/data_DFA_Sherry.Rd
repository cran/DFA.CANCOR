\name{data_DFA_Sherry}
\alias{data_DFA_Sherry}
\docType{data}
\title{data_DFA_Sherry}
\description{A data frame with scores on 5 variables for 10 cases in each of 3 groups. 
	Used by Sherry (2006) in her article on discriminant function analysis.}
\usage{data(data_DFA_Sherry)}
\source{Sherry, A. (2006). Discriminant analysis in counseling research. 
	\emph{Counseling Psychologist, 34,} 661-683.}
\examples{
\donttest{
head(data_DFA_Sherry)
     
DFA(data = data_DFA_Sherry, 
    groups = 'Group',
    variables = c('Neuroticism','Extroversion','Openness',
                  'Agreeableness','Conscientiousness'),
    predictive = TRUE, priorprob = 'SIZES', 
    verbose=TRUE)
}
}
