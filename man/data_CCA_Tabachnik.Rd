\name{data_CCA_Tabachnik}
\alias{data_CCA_Tabachnik}
\docType{data}
\title{data_CCA_Tabachnik}
\description{A data frame with scores on 4 variables for 8 cases. Used by 
	Tabachnik & Fidell (2019, p. 451) in their chapter on canonical correlation.}
\usage{data(data_CCA_Tabachnik)}
\source{Tabachnik, B. G., & Fidell, L. S. (2019). \emph{Using multivariate 
	statistics (7th ed.).} New York, NY: Pearson.}
\examples{
\donttest{
head(data_CCA_Tabachnik)
	
CANCOR(data = data_CCA_Tabachnik,
       set1 = c('TS','TC'),
       set2 = c('BS','BC'),
       plot = 'yes', plotCV = 1,
       verbose=TRUE)	
}
}