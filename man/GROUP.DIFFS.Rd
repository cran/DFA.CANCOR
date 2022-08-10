\name{GROUP.DIFFS}
\alias{GROUP.DIFFS}
\title{Group Mean Differences on a Continuous Outcome Variable}
\description{Produces a variety of statistics for all possible pairwise independent groups comparisons 
	of means on a continuous outcome variable.}
\usage{GROUP.DIFFS(data, GROUPS=NULL, DV=NULL, var.equal=FALSE, p.adjust.method="holm", 
	            Ncomps=NULL, verbose=TRUE)}
\arguments{
  \item{data}{
  \code{}A dataframe where the rows are cases & the columns are the variables. If GROUPS and DV are not 
  specified, then the GROUPS variable should be in the first column and the DV should be in the second
  column of data.}

  \item{GROUPS}{
  \code{}The name of the groups variable in the dataframe, e.g., groups = 'Group'.}
  	    
  \item{DV}{
  \code{}The name of the dependent (outcome) variable in the dataframe, e.g., DV = 'esteem'.}

  \item{var.equal}{ 
  \code{}(from stats::t.test) A logical variable indicating whether to treat the two variances as being equal.  
  If TRUE then the pooled variance is used to estimate the variance otherwise the Welch (or Satterthwaite) 
  approximation to the degrees of freedom is used.}
  	    
  \item{p.adjust.method}{
  \code{}The method to be used to adjust the p values for the number of comparisons. The options are "holm" 
  (the default), "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".}

  \item{Ncomps}{
  \code{}The number of pairwise comparisons for the adjusted p values. If unspecified, it will be the 
  number of all possible comparisons (i.e., the family-wise number of number of comparisons). Ncomps
  could alternatively be set to, e.g., the experiment-wise number of number of comparisons.}

  \item{verbose}{
  \code{}Should detailed results be displayed in console? \cr The options are: TRUE (default) or FALSE.}
}
\details{The function conducts all possible pairwise comparisons of the levels of the GROUPS variable on
	the continuous outcome variable. It supplements independent groups t-test results with effect size statistics
	and with the Bayes factor for each pairwise comparison.
	
	The d values are the Cohen d effect sizes, i.e., the mean difference expressed in standard deviation units.
	
	The g values are the Hedges g value corrections to the Cohen d effect sizes.
	
	The r values are the effect sizes for the group mean difference expressed in the metric of Pearson's r.
	
	The BESD values are the binomial effect size values for the group mean differences.
    The BESD casts the effect size in terms of the success rate for the implementation of a hypothetical procedure
	(e.g., the percentage of cases that were cured, or who died.) For example, an r = .32 is equivalent to increasing 
	the success rate from 34\% to 66\% (or, possibly, reducing an illness or death rate from 66\% to 34\%).
	
	The Bayes factor values are obtained from the ttest.tstat function in the BayesFactor package.

	For example, a Bayes_Factor_alt_vs_null = 3 indicates that the data are 3 times \emph{more} likely under the  
	alternative hypothesis than under the null hypothesis. A Bayes_Factor_alt_vs_null = .2 indicates that 
	the data are five times \emph{less} likely under the alternative hypothesis than under the null hypothesis (1 / .2).

	Conversely, a Bayes_Factor_null_vs_alt = 3 indicates that the data are 3 times \emph{more} likely under the  
	null hypothesis than under the alternative hypothesis. A Bayes_Factor_null_vs_alt = .2 indicates that 
	the data are five times \emph{less} likely under the null hypothesis than under the alternative hypothesis (1 / .2).
}
\value{If verbose = TRUE, the displayed output includes the means, standard deviations, and Ns for the groups, 
	the t-test results for each pairwise comparison, the mean difference and its 95\% confidence interval, four 
	indices of effect size for each pairwise comparison (r, d, g, and BESD), and the Bayes factor. The returned 
	output is a matrix with these values.
}
\references{
       {Funder, D. C., & Ozer, D. J. (2019). Evaluating effect size in psychological research: Sense and nonsense. 
       	\emph{Advances in Methods and Practices in Psychological Science, 2(2),} 156168.}
\cr\cr {Jarosz, A. F., & Wiley, J. (2014). What are the odds? A practical guide to computing and reporting 
	    Bayes factors. \emph{Journal of Problem Solving, 7,} 29.}
\cr\cr {Randolph, J. & Edmondson, R.S. (2005). Using the binomial effect size display (BESD) to present the 
	    magnitude of effect sizes to the evaluation audience. \emph{Practical Assessment Research & Evaluation, 10,} 14.}
\cr\cr {Rosenthal, R., Rosnow, R.L., & Rubin, D.R. (2000). \emph{Contrasts and effect sizes in behavioral research: 
	    A correlational approach.} Cambridge UK: Cambridge University Press.}
\cr\cr {Rosenthal, R., & Rubin, D. B. (1982). A simple general purpose display of magnitude and experimental effect. 
	   \emph{Journal of Educational Psychology, 74,} 166-169.}
\cr\cr {Rouder, J. N., Haaf, J. M., & Vandekerckhove, J. (2018). Bayesian inference for psychology, part IV:
       	parameter estimation and Bayes factors. \emph{Psychonomic Bulletin & Review, 25(1),} 102113.}
}
\author{Brian P. O'Connor}
\examples{
GROUP.DIFFS(data_DFA$Field_2012, var.equal=FALSE, p.adjust.method="fdr")

\donttest{
GROUP.DIFFS(data = data_DFA$Sherry_2006, var.equal=FALSE, p.adjust.method="bonferroni")    
}
}
	% The Bayes factor tests the null hypothesis that the mean difference of a normal population is 0. 
	% Specifically, the Bayes factor compares two hypotheses: that the standardized effect size is 0, or that the 
	% standardized effect size is not 0. The standardized effect size in this case, for a two sample test, 
	% is (mu2-mu1)/sigma.
