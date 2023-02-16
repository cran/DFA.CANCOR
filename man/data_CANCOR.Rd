\name{data_CANCOR}
\alias{data_CANCOR}
\docType{data}
\title{data_CANCOR}
\description{A list with example data that were used in various
	presentations of canonical correlation analysis
}
\usage{data(data_CANCOR)}
\details{
A list with the example data that were used in the following presentations 
of canonical correlation analysis: De Leo and Wulfert (2013), Ho (2014), 
Rencher (2002), Tabachnick and Fidell (2019), and by the UCLA statistics tutorial 
at https://stats.oarc.ucla.edu/r/dae/canonical-correlation-analysis/.
}

\references{
	{De Leo, J. A., & Wulfert, E. (2013). Problematic internet use and other risky 
	behaviors in college students: An application of problem-behavior theory. \emph{Psychology 
	of Addictive Behaviors, 27(1),} 133-141. }
\cr\cr {Ho, R. (2014). \emph{Handbook of univariate and multivariate data analysis with 
	IBM SPSS.} Boca Raton, FL: CRC Press.}
\cr\cr {Rencher, A. (2002). \emph{Methods of multivariate analysis} (2nd ed.). 
	New York, NY: John Wiley & Sons.}
\cr\cr {Tabachnick, B. G., & Fidell, L. S. (2019). Chapter 16: Multiway 
	frequency analysis. \emph{Using multivariate statistics.} New York, NY: Pearson.}
}

\examples{
names(data_CANCOR)
\donttest{ 
head(data_CANCOR$DeLeo_2013)

head(data_CANCOR$Ho_2014)

head(data_CANCOR$Rencher_2002)

head(data_CANCOR$TabFid_2019_small)

head(data_CANCOR$TabFid_2019_complete)
}
}
