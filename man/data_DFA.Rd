\name{data_DFA}
\alias{data_DFA}
\docType{data}
\title{data_DFA}
\description{A list with example data that were used in various
	presentations of discrimination function analysis
}
\usage{data(data_DFA)}
\details{
A list with the example data that were used in the following presentations 
of discrimination function analysis: Field (2012), Green and Salkind (2008), 
Ho (2014), Huberty and Olejnik (2006), Noursis (2012), Rencher (2002), Sherry (2006), and
Tabachnick and Fidell (2019).
}

\references{
	{Field, A., Miles, J., & Field, Z. (2012). Chapter 18 Categorical data. 
		\emph{Discovering statistics using R.} Los Angeles, CA: Sage.}
\cr\cr {Green, S. B., & Salkind, N. J. (2008). Lesson 35: Discriminant analysis
(pp. 300-311). In, \emph{Using SPSS for Windows and Macintosh: Analyzing and 
understanding data.} New York, NY: Pearson.}
\cr\cr {Ho, R. (2014). \emph{Handbook of univariate and multivariate data analysis with 
	IBM SPSS.} Boca Raton, FL: CRC Press.}
\cr\cr {Huberty, C. J., & Olejnik, S. (2019). \emph{Applied MANOVA and discriminant 
	analysis} (2nd. ed.). New York, NY: John Wiley & Sons.}
\cr\cr {Noursis, M. J. (2012). \emph{IBM SPSS Statistics 19 advanced statistical 
	procedures companion.} Upper Saddle River, NJ: Prentice Hall.}
\cr\cr {Rencher, A. (2002). \emph{Methods of multivariate analysis} (2nd ed.). 
	New York, NY: John Wiley & Sons.}
\cr\cr {Sherry, A. (2006). Discriminant analysis in counseling research.
	\emph{Counseling Psychologist, 34,} 661-683.}
\cr\cr {Tabachnick, B. G., & Fidell, L. S. (2019). Chapter 16: Multiway 
	frequency analysis. \emph{Using multivariate statistics.} New York, NY: Pearson.}
}

\examples{
names(data_DFA)
\donttest{ 
head(data_DFA$Field_2012)

head(data_DFA$Green_2008)

head(data_DFA$Ho_2014)

head(data_DFA$Huberty_2019_p45)

head(data_DFA$Huberty_2019_p285)

head(data_DFA$Norusis_2012)

head(data_DFA$Rencher_2002_football)

head(data_DFA$Rencher_2002_root)

head(data_DFA$Sherry_2006)

head(data_DFA$TabFid_2019_complete)

head(data_DFA$TabFid_2019_small)
}
}
