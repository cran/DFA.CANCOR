\name{DFA}
\alias{DFA}
\title{Discriminant function analysis}
\description{Produces SPSS- and SAS-like output for linear discriminant 
	function analysis. It uses functions from the MASS package.}
\usage{DFA(data, groups, variables, plot, predictive, priorprob, verbose)}
\arguments{
  \item{data}{
  \code{}A dataframe where the rows are cases & the columns are the variables.}

  \item{groups}{
  \code{}The name of the groups variable in the dataframe, \cr e.g., groups = 'Group'.}
  	    
  \item{variables}{ 
  \code{}The names of the continuous variables in the dataframe that will be used in the DFA, 
         e.g., variables = c('varA', 'varB', 'varC').}
  	    
  \item{plot}{ 
  \code{}Should a plot of the mean standardized discriminant function scores  
         \cr for the groups be produced? The options are: TRUE (default) or FALSE.}
  	    
  \item{predictive}{
  \code{}Should a predictive DFA be conducted? \cr The options are: TRUE (default) or FALSE.}
  
  \item{priorprob}{
  \code{}If predictive = TRUE, how should the prior probabilities of the group sizes be computed? 
  		 The options are:\cr 
  	    'EQUAL' for equal group sizes; or\cr 
  	    'SIZES' (default) for the group sizes to be based on the sizes of the groups in the dataframe.}
                  
  \item{verbose}{
  \code{}Should detailed results be displayed in console? \cr The options are: TRUE (default) or FALSE.}
  
%  \item{CV}{
%  \code{}Should cross-validation (leave-one-out cross-validation) be conducted for the predictive DFA? 
%         The options are: TRUE (default) or FALSE.}

}
\value{If verbose = TRUE, the displayed output includes descriptive statistics for the groups, 
	   tests of univariate and multivariate normality, 
	   the results of tests of the homogeneity of the group variance-covariance matrices, 
	   eigenvalues & canonical correlations, Wilks lambda & peel-down statistics,
	   raw and standardized discriminant function coefficients, 
	   structure coefficients, functions at group centroids, 
	   one-way ANOVA tests of group differences in scores on each discriminant function,
	   one-way ANOVA tests of group differences in scores on each original DV,
	   significance tests for group differences on the original DVs according to Bird et al. (2014),
	   a plot of the group means on the standardized discriminant functions,
	   and extensive output from predictive discriminant function analyses (if requested).
	   		   	   
	   The returned output is a list with elements
  \item{rawCoef}{canonical discriminant function coefficients}
  \item{structCoef}{structure coefficients}
  \item{standCoef}{standardized coefficients}
  \item{standCoefSPSS}{standardized coefficients from SPSS}
  \item{centroids}{unstandardized canonical discriminant functions evaluated at the group means}
  \item{centroidSDs}{group standard deviations on the unstandardized functions}
  \item{centroidsZ}{standardized canonical discriminant functions evaluated at the group means}
  \item{centroidSDsZ}{group standard deviations on the standardized functions}
  \item{DFAscores}{scores on the discriminant functions}
  \item{anovaDFoutput}{One-way ANOVAs using the scores on a discriminant function as the DV}
  \item{anovaDVoutput}{One-way ANOVAs on the original DVs}
  \item{MFWER1.sigtest}{Significance tests when controlling the MFWER by (only) carrying out multiple t tests}
  \item{MFWER2.sigtest}{Significance tests for the two-stage approach to controling the MFWER}
  \item{ldaoutputCV}{Classifications from leave-one-out cross-validations}
  \item{freqs_OR}{Cross-Tabulation of the Original and Predicted Group Memberships}
  \item{PropOrigCorrect}{Proportion of original grouped cases correctly classified}
  \item{chi_square_OR}{Chi-square test of independence}
  \item{PressQ_OR}{Press's Q significance test of classifiation accuracy for original vs. predicted group memberships}
  \item{rowfreqs_OR}{Row Frequencies}
  \item{colfreqs_OR}{Column Frequencies}
  \item{cellprops_OR}{Cell Proportions}
  \item{rowprops_OR}{Row-Based Proportions}
  \item{colprops_OR}{Column-Based Proportions}
  \item{kappas_cvo_OR}{Agreement (kappas) between the Predicted and Original Group Memberships}
  \item{freqs_CV}{Cross-Tabulation of the Cross-Validated and Predicted Group Memberships}
  \item{PropCrossValCorrect}{Proportion of cross-validated grouped cases correctly classified}
  \item{chi_square_CV}{Chi-square test of indepedence}
  \item{PressQ_CV}{Press's Q significance test of classifiation accuracy for cross-validated vs. predicted group memberships}
  \item{rowfreqs_CV}{Row Frequencies}
  \item{colfreqs_CV}{Column Frequencies}
  \item{cellprops_CV}{Cell Proportions}
  \item{rowprops_CV}{Row-Based Proportions}
  \item{colprops_CV}{Column-Based Proportions}
  \item{kappas_cvoCV}{Agreement (kappas) between the Cross-Validated and Original Group Membership}
  \item{kappas_CVP}{Agreement (kappas) between the Cross-Validated and Predicted Group Memberships}
}
\references{
       {Bird, K. D., & Hadzi-Pavlovic, D. (2013). Controlling the maximum familywise Type I error 
       	rate in analyses of multivariate experiments. \emph{Psychological Methods, 19(2),} p. 265-280.}
\cr\cr {Manly, B. F. J., & Alberto, J. A. (2017). \emph{Multivariate statistical methods: 
	   	A primer (4th Edition).} Chapman & Hall/CRC, Boca Raton, FL.}
\cr\cr {Sherry, A. (2006). Discriminant analysis in counseling research. \emph{Counseling Psychologist, 34,} 661-683.}
\cr\cr {Tabachnik, B. G., & Fidell, L. S. (2013). \emph{Using multivariate statistics (6th ed.).} New York, NY: Pearson.}
\cr\cr {Venables, W. N. & Ripley, B. D. (2002). \emph{Modern Applied Statistics with S (4th ed.).} Springer, New York.}
}
\author{Brian P. O'Connor }
\examples{
DFA(data = data_DFA_Field, 
    groups = 'Group', 
    variables = c('Actions','Thoughts'),
    predictive = TRUE, priorprob = 'SIZES',
    verbose = TRUE)

\donttest{
DFA(data = data_DFA_Sherry, 
    groups = 'Group',
    variables = c('Neuroticism','Extroversion','Openness', 
                  'Agreeableness','Conscientiousness'),
    predictive = TRUE, priorprob = 'SIZES', 
    verbose = TRUE)
}
}
