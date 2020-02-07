\name{homogeneity}
\alias{homogeneity}
\title{Homogeneity of variances and covariances}
\description{P tests of the homogeneity of variances and covariances.}
\usage{homogeneity(data, groups, variables, verbose)}
\arguments{
\item{data}{
\code{}A dataframe where the rows are cases & the columns are the variables.}

\item{groups}{
\code{}(optional) The name of the groups variable in the dataframe (if there is one) \cr e.g., groups = 'Group'.}

\item{variables}{
\code{}(optional) The names of the continuous variables in the dataframe for the analyses,
e.g., variables = c('varA', 'varB', 'varC').}

\item{verbose}{
\code{}Should detailed results be displayed in the console? \cr The options are: TRUE (default) or FALSE.}
}
\value{If "variables" is specified, the analyses will be run on the "variables" in "data".
If verbose = TRUE, the displayed output includes descriptive statistics and
tests of univariate and multivariate homogeneity.

Bartlett's test compares the variances of k samples. The data must be normally distributed.

The non-parametric Fligner-Killeen test also compares the variances of k samples and
it is robust when there are departures from normality.

Box's M test is a multivariate statistical test of the equality of multiple
variance-covariance matrices. The test is prone to errors when the sample sizes are small or when
the data do not meet model assumptions, especially the assumption of multivariate normality.
For large samples, Box's M test may be too strict, indicating heterogeneity when the covariance 
matrices are not very different.

The returned output is a list with elements
\item{covmatrix}{The variance-covariance matrix for each group}
\item{Bartlett}{Bartlett test of homogeneity of variances (parametric)}
\item{Figner_Killeen}{Figner-Killeen test of homogeneity of variances (non parametric)}
\item{PooledWithinCovarSPSS}{the pooled within groups covariance matrix from SPSS}
\item{PooledWithinCorrelSPSS}{the pooled within groups correlation matrix from SPSS}
\item{sscpWithin}{the within sums of squares and cross-products matrix}
\item{sscpBetween}{the between sums of squares and cross-products matrix}
\item{BoxLogdets}{the log determinants for Box's test}
\item{BoxMtest}{Box's' test of the equality of covariance matrices}
}
\references{
{Box, G. E. P. (1949). A general distribution theory for a class of likelihood
criteria. \emph{Biometrika, 36 (3-4),} 317-346.}

{Bartlett, M. S. (1937). Properties of sufficiency and statistical tests. \emph{Proceedings of the Royal Society
of London Series A 160,} 268-282.}

{Conover, W. J., Johnson, M. E., & Johnson, M. M. (1981). A comparative study of tests for homogeneity
of variances, with applications to the outer continental shelf bidding data. \emph{Technometrics, 23,} 351-361.}

{Warner, R. M. (2013). \emph{Applied statistics: From bivariate through multivariate techniques.} Thousand Oaks, CA: SAGE.}
}
\author{Brian P. O'Connor }
\examples{
# data from Field et al. (2012)
homogeneity(data = data_DFA_Field,
            groups = 'Group', variables = c('Actions','Thoughts'))

\donttest{
# data from Sherry (2006)
homogeneity(data = data_DFA_Sherry,
            groups = 'Group',
            variables = c('Neuroticism','Extroversion','Openness',
                          'Agreeableness','Conscientiousness'))
}

}