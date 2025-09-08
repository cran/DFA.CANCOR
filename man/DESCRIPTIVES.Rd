\name{DESCRIPTIVES}
\alias{DESCRIPTIVES}
\title{Descriptive statistcs for numeric variables}
\description{Produces descriptive statistcs for numeric variables, possibly by a grouping variable.}
\usage{DESCRIPTIVES(data, groups, variables, CI_level = 95, verbose)}
\arguments{
  \item{data}{
  \code{}A dataframe or numeric matrix where the rows are cases & the columns are the variables.}

  \item{groups}{
  \code{}(optional) The name of the groups variable in the dataframe, if there is one, \cr e.g., groups = 'Group'.}
  	    
  \item{variables}{ 
  \code{}(optional) The names of the continuous variables in the dataframe for the analyses, 
         e.g., variables = c('varA', 'varB', 'varC').}
  	    
 \item{CI_level}{
  \code{}(optional) The confidence interval for the output, in whole numbers, e.g., CI_level = 95.
  The default is 95.}

  \item{verbose}{
  \code{}Should detailed results be displayed in the console? 
  \cr The options are: TRUE (default) or FALSE.}
}
\details{If "groups" is not specified, the analyses will be run on all of the variables in "data".
     If "variables" is specified, the analyses will be run on the "variables" in "data".
	   If "groups" is specified, the analyses will be run for every value of "groups".
}
\value{The returned output is a list with the descriptive statistics.}
\author{Brian P. O'Connor }
\examples{
# without a grouping variable
DESCRIPTIVES(data = data_DFA$Field_2012, variables = c('Actions','Thoughts'))

\donttest{
# with a grouping variable
DESCRIPTIVES(data = data_DFA$Field_2012, 
          groups = 'Group', 
          variables = c('Actions','Thoughts'))
}
}
