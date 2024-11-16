\name{OUTLIERS}
\alias{OUTLIERS}
\title{OUTLIERS}
\description{Provides tests and qqplots for multivariate outliers.}
\usage{OUTLIERS(data, variables, ID=NULL, iterate=TRUE,
            alpha_univ=.05, plot_univariates=TRUE,
            MCD=TRUE, MCD.quantile = .75, alpha=0.025, cutoff_type = 'adjusted',
            qqplot=TRUE, plot_iters=NULL, 
            verbose=TRUE)
}
\arguments{
 \item{data}{
  \code{}A dataframe where the rows are cases & the columns are the variables.}

 \item{variables}{
  \code{}The names of the continuous variables in the dataframe for the analyses, 
         e.g., variables = c('varA', 'varB', 'varC').}
  	    
 \item{ID}{
  \code{}(optional) The names of the case identification variable in data, if
  there is one. If ID is not specified, then the sequence of row numbers will 
  be used as the case IDs.}
  
 \item{iterate}{
  \code{}(optional) Should multiple iterations be conducted when searching for
         outliers? \cr The options are: TRUE (default) or FALSE.}

  \item{alpha_univ}{
  \code{}(optional) The p (alpha) level for univariate outliers.
       The default = .05.}

 \item{plot_univariates}{
  \code{}(optional) Should univariate plots be provided? 
  \cr The options are: TRUE (default) or FALSE.}

  \item{MCD}{
  \code{}(optional) Should the Minimum Covariance Determinant method be used
         to compute the means and covariances?
         \cr The options are: TRUE (default) or FALSE.}
  	    
  \item{MCD.quantile}{
  \code{}(optional) The MCD quantile, which is the the minimum 
         number of the data points regarded as good points (MASS package).
         The default = .75, as recommended by Leys et al. (2018).}
  	    
 \item{alpha}{
  \code{}(optional) alpha}
  	    
  \item{cutoff_type}{
  \code{}(optional) The kind of cutoff to be computed. The options are
         adjusted' (the default) or 'quan'.}
  	    
  \item{qqplot}{
  \code{}(optional) Should qqplots be provided? 
       \cr The options are: TRUE (default) or FALSE.}

  \item{plot_iters}{
  \code{}(optional) A vector with the iterations for the qqplot. For example,
  "plot_iters = c(1,2,6,7)" will produce a qqplot for each of iterations 1, 2, 6, and 7
  on the output figure. The default is "plot_iters = c(1,2,3,4)".}

  \item{verbose}{
  \code{}(optional) Should detailed results be displayed in console? TRUE (default) or FALSE}
}
\details{This function provides both statistical and graphical methods of
identifying multivariate outliers. Both methods are based on
Mahalanobis distances.

A Mahalanobis distance is an estimate of how far each case is from the
center of the joint distribution of the variables in multivariate
space. Cases that are distant from the swarm of most other cases may
be multivariate outliers.

Squared Mahalanobis distances have an approximate chi-squared
distribution (when there is multivariate normality). Statistically, a
multivariate outlier is said to exist when the squared Mahalanobis
distance for a case is greater than a specified cut-off value that is
derived from the chi-square distribution.

The computations for Mahalanobis distances are based on estimates of
the means and covariances for the dataset. However, the means and
covariances that are based on all of the data are affected by the
existence of multivariate outliers. This renders the simple,
whole-sample estimates of Mahalanobis distances, and thus the
identification of outliers, problematic.

Better estimates of the means and covariances are obtained using the
Minimum Covariance Determinant (MCD) method, which identifies the most
central subset of the data. Mahalanobis distances are considered more
"robust" when they are computed using the MCD means and covariances.
The default for the \strong{MCD argument} for this function is set to TRUE
for this reason. Setting it to FALSE will result in the procedure
using the whole-sample based means and covariances, which is not
recommended.

Once obtained, Mahalanobis distances (robust or not) are assessed for
statistical significance by comparing them with a specified quantile
from the chi-squared distribution. There are two options for
determining the specified quantile cutoff value. The simple,
traditional approach is to use the alpha quantile of the chi-squared
distribution with the degrees of freedom equal to the number of
variables. In the present function, the default alpha threshold is
0.025.

A modern, alternative method of determining cutoff values is to use the
adaptive reweighted estimator procedure (Filzmoser, Garrett, &
Reimann, 2005), which derives a cutoff value that is appropriate for
each specific dataset and sample size. These threshold values are
called "adjusted quantiles".

The \strong{cutoff_type argument} for this function can be set to "adjusted"
for an adjusted quantile, or to "quan" for the traditional alpha
quantile.

A "qqplot" of the squared Mahalanobis distances can be used to
graphically assess multivariate normality and the existence of
outliers. In this case, the (sorted) squared Mahalanobis distances are
plotted against the corresponding quantiles of the chi-square
distribution. When the the squared Mahalanobis distances fit the
hypothesized distribution, the points in the Q-Q plot will fall on a
straight, y = x line (chi-squared values are squared z scores).
Deviations from the straight line suggest violations of multivariate
normality and the possible existence of multivariate outliers.

The search for multivariate outliers can be conducted more than once
for a given dataset. If outliers are identified on the first step
(iteration), they can be removed from the dataset and another search
for outliers can be conducted on the remaining data. It is not
uncommon for multiple iterations to be required before no further
outliers are found. Bigger outliers can mask smaller but still
possibly important outliers. It is probably best to run the analyses
for multiple iterations. In the present function, multiple iterations
are conducted when the \strong{iterate argument} is set to TRUE.

The present function provides up to four possible qqplots in the
one-page output figure for a data analysis. By default, these plots
will be for the first four interations that produced outliers.  Use
the \strong{plot_iters argument} to produce plots from alternative iterations.
For example, "plot_iters = c(1,2,6,7)" will place the qqplots from
iterations 1, 2, 6, and 7 on the output figure.
}
\value{The returned output is a list with the outliers.
}
\references{
       {Filzmoser, P., Garrett, R. G., & Reimann, C. (2005). Multivariate outlier 
       detection in exploration geochemistry. \emph{Computers & Geosciences, 
       31,} 579-587.}
\cr\cr {Leys, C., Klein, O., Dominicy, Y., & Ley, C. (2018). Detecting 
        multivariate outliers: Use a robust variant of the Mahalanobis distance. 
        \emph{Journal of Experimental Social Psychology, 74,} 150-156.}
\cr\cr {Rodrigues, I. M., & Boente, G. (2011). Multivariate outliers. 
        \emph{International Encyclopedia of Statistical Science} (pp. 910-912). 
        Berlin:Springer-Verlag.}
\cr\cr {Rousseeuw, P. J., & Leroy, A. M. (1987). \emph{Robust Regression and Outlier 
        Detection}. New York, NY: John Wiley & Sons.}
}
\author{Brian P. O'Connor }
\examples{
OUTLIERS(data = iris, variables = c('Sepal.Length','Sepal.Width','Petal.Length'), 
         ID=NULL, iterate=TRUE,
         alpha_univ=.05, plot_univariates=TRUE,
         MCD=TRUE, MCD.quantile = .75, alpha=0.025, cutoff_type = 'adjusted',
         qqplot=TRUE, plot_iters=c(1,2,5,6), 
         verbose=TRUE)
\donttest{

}
}

