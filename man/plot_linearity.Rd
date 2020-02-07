\name{plot_linearity}
\alias{plot_linearity}
\title{Plot for linearity}
\description{Plots the linear, quadratic, and loess regression lines for the
	association between two continuous variables.}
\usage{plot_linearity(data, idv, dv, groups=NULL, groupNAME=NULL, verbose = TRUE)}
\arguments{
  \item{data}{
  \code{}A dataframe where the rows are cases & the columns are the variables.}

  \item{idv}{
  \code{}The name of the predictor variable.}
  	    
  \item{dv}{
  \code{}The name of the dependent variable.}
  	    
  \item{groups}{
  \code{}(optional) The name of the groups variable in the dataframe, 
  		 \cr e.g., groups = 'Group'.}
  	    
  \item{groupNAME}{
  \code{}(optional) The value (level, name, or number) from the groups 
  	variable that identifies the subset group whose data will be used 
  	for the analyses, \cr e.g., groupNAME = 1.}
  	    
  \item{verbose}{
  \code{}Should detailed results be displayed in the console? 
  		\cr The options are: TRUE (default) or FALSE.}
}
\value{If verbose = TRUE, the linear and quadratic regression
	   coefficients and their statistical tests are displayed.
	   	   
	   The returned output is a list with the regression 
	   coefficients and the plot data.}
\references{
 {Tabachnik, B. G., & Fidell, L. S. (2013). \emph{Using multivariate 
 		statistics (6th. ed., pp. 83-84).} New York, NY: Pearson.}
}
\author{Brian P. O'Connor }
\examples{
# data that simulate those from De Leo & Wulfert (2013)
plot_linearity(data=data_CCA_De_Leo, groups=NULL, 
               idv='Family_Conflict', dv='Grade_Point_Average', verbose=TRUE)

\donttest{
# data from Sherry (2006), ignoring the groups
plot_linearity(data=data_DFA_Sherry, groups=NULL, groupNAME=NULL,
               idv='Neuroticism', dv='Conscientiousness', verbose=TRUE)

# data from Sherry (2006), group 2 only
plot_linearity(data=data_DFA_Sherry, groups ='Group', groupNAME=2,
               idv='Neuroticism', dv='Conscientiousness', verbose=TRUE)
 }

}
