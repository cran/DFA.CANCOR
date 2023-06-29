\name{DFA}
\alias{DFA}
\title{Discriminant function analysis}
\description{Produces SPSS- and SAS-like output for linear discriminant function analysis.}
\usage{DFA(data, groups, variables, plot, predictive, priorprob, covmat_type, CV, verbose)}
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
                  
  \item{covmat_type}{
  \code{}The kind of covariance to be used for a predictive DFA. The options are: 
  \cr 'within' (for the pooled within-groups covariance matrix, which is the default) or 
  \cr 'separate' (for separate-groups covariance matrices).}
  
 \item{CV}{
 \code{}If predictive = TRUE, should cross-validation (leave-one-out cross-validation) analyses also be conducted? 
        The options are: TRUE (default) or FALSE.}

  \item{verbose}{
  \code{}Should detailed results be displayed in console? \cr The options are: TRUE (default) or FALSE.}
}
\details{The predictive DFA option using separate-groups covariance matrices (which is often called 'quadratic DFA') 
	is conducted following the procedures described by Rencher (2002). The covariance matrices 
	in this case are based on the scores on the continuous variables. In contrast, the 'separate-groups'
	option in SPSS involves use of the group scores on the discriminant functions (not the original continuous variables), 
	which can produce different classifications.
	
	When data has many cases (e.g., > 1000), the leave-one-out cross-validation analyses can
	be time-consuming to run. Set CV = FALSE to bypass the predictive DFA cross-validation analyses.
	
	See the documentation below for the GROUP.DIFFS function for information on the interpretation of the
	Bayes factors and effect sizes that are produced for the group comparisons.
}
\value{If verbose = TRUE, the displayed output includes descriptive statistics for the groups, 
	   tests of univariate and multivariate normality, 
	   the results of tests of the homogeneity of the group variance-covariance matrices, 
	   eigenvalues & canonical correlations, Wilks' lambda & peel-down statistics,
	   raw and standardized discriminant function coefficients, 
	   structure coefficients, functions at group centroids, 
	   one-way ANOVA tests of group differences in scores on each discriminant function,
	   one-way ANOVA tests of group differences in scores on each original DV,
	   significance tests for group differences on the original DVs according to Bird et al. (2014),
	   a plot of the group means on the standardized discriminant functions,
	   and extensive output from predictive discriminant function analyses (if requested).
	   		   	   
	   The returned output is a list with elements

  \item{evals}{eigenvalues and canonical correlations}
  \item{mv_Wilks}{The Wilks' lambda multivariate test}
  \item{mv_Pillai}{The Pillai-Bartlett multivariate test}
  \item{mv_Hotelling}{The Lawley-Hotelling multivariate test}
  \item{mv_Roy}{Roy's greatest characteristic root multivariate test}
  \item{coefs_raw}{canonical discriminant function coefficients}
  \item{coefs_structure}{structure coefficients}
  \item{coefs_standardized}{standardized coefficients}
  \item{coefs_standardizedSPSS}{standardized coefficients from SPSS}
  \item{centroids}{unstandardized canonical discriminant functions evaluated at the group means}
  \item{centroidSDs}{group standard deviations on the unstandardized functions}
  \item{centroidsZ}{standardized canonical discriminant functions evaluated at the group means}
  \item{centroidSDsZ}{group standard deviations on the standardized functions}
  \item{dfa_scores}{scores on the discriminant functions}
  \item{anovaDFoutput}{One-way ANOVAs using the scores on a discriminant function as the DV}
  \item{anovaDVoutput}{One-way ANOVAs on the original DVs}
  \item{MFWER1.sigtest}{Significance tests when controlling the MFWER by (only) carrying out multiple t tests}
  \item{MFWER2.sigtest}{Significance tests for the two-stage approach to controling the MFWER}  
  \item{classes_PRED}{The predicted group classifications}  
  \item{classes_CV}{The classifications from leave-one-out cross-validations, if requested}  
  \item{posteriors}{The posterior probabilities for the predicted group classifications}
  \item{grp_post_stats}{Group mean posterior classification probabilities}
  \item{classes_CV}{Classifications from leave-one-out cross-validations}
  \item{freqs_ORIG_PRED}{Cross-tabulation of the original and predicted group memberships}
  \item{chi_square_ORIG_PRED}{Chi-square test of independence}
  \item{PressQ_ORIG_PRED}{Press's Q significance test of classifiation accuracy for original vs. predicted group memberships}
  \item{kappas_ORIG_PRED}{Agreement (kappas) between the predicted and original group memberships}
  \item{PropOrigCorrect}{Proportion of original grouped cases correctly classified}  
  \item{freqs_ORIG_CV}{Cross-Tabulation of the cross-validated and predicted group memberships}
  \item{chi_square_ORIG_CV}{Chi-square test of indepedence}
  \item{PressQ_ORIG_CV}{Press's Q significance test of classifiation accuracy for cross-validated vs. predicted group memberships}
  \item{kappas_ORIG_CV}{Agreement (kappas) between the cross-validated and original group memberships}
  \item{PropCrossValCorrect}{Proportion of cross-validated grouped cases correctly classified}
}
\references{
       {Bird, K. D., & Hadzi-Pavlovic, D. (2013). Controlling the maximum familywise Type I error 
       	rate in analyses of multivariate experiments. \emph{Psychological Methods, 19(2),} p. 265-280.}
\cr\cr {Manly, B. F. J., & Alberto, J. A. (2017). \emph{Multivariate statistical methods: 
	   	A primer (4th Edition).} Chapman & Hall/CRC, Boca Raton, FL.}
\cr\cr {Rencher, A. C. (2002). \emph{Methods of Multivariate Analysis} (2nd ed.). New York, NY: John Wiley & Sons.}
\cr\cr {Sherry, A. (2006). Discriminant analysis in counseling research. \emph{Counseling Psychologist, 34,} 661-683.}
\cr\cr {Tabachnik, B. G., & Fidell, L. S. (2019). \emph{Using multivariate statistics (7th ed.).} New York, NY: Pearson.}
}
\author{Brian P. O'Connor}
\examples{
# data from Field et al. (2012, Chapter 16 MANOVA)
DFA_Field=DFA(data = data_DFA$Field_2012, 
    groups = 'Group', 
    variables = c('Actions','Thoughts'),
    predictive = TRUE, 
    priorprob = 'EQUAL',   
    covmat_type='within', # altho better to use 'separate' for these data
    verbose = TRUE)

\donttest{

# plots of posterior probabilities by group
# hoping to see correct separations between cases from different groups

# first, display the posterior probabilities
print(cbind(round(DFA_Field$posteriors[1:3],3), DFA_Field$posteriors[4]))

# group NT vs CBT
plot(DFA_Field$posteriors$posterior_NT, DFA_Field$posteriors$posterior_CBT, 
     pch = 16, col = c('red', 'blue', 'green')[DFA_Field$posteriors$Group],
     xlim=c(0,1), ylim=c(0,1),
     main = 'DFA Posterior Probabilities by Original Group Memberships',
     xlab='Posterior Probability of Being in Group NT',
     ylab='Posterior Probability of Being in Group CBT' )
legend(x=.8, y=.99, c('CBT','BT','NT'), cex=1.2, col=c('red', 'blue', 'green'), pch=16, bty='n')

# group NT vs BT
plot(DFA_Field$posteriors$posterior_NT, DFA_Field$posteriors$posterior_BT, 
     pch = 16, col = c('red', 'blue', 'green')[DFA_Field$posteriors$Group],
     xlim=c(0,1), ylim=c(0,1),
     main = 'DFA Posterior Probabilities by Group Membership',
     xlab='Posterior Probability of Being in Group NT',
     ylab='Posterior Probability of Being in Group BT' )
legend(x=.8, y=.99, c('CBT','BT','NT'), cex=1.2,col=c('red', 'blue', 'green'), pch=16, bty='n')

# group CBT vs BT
plot(DFA_Field$posteriors$posterior_CBT, DFA_Field$posteriors$posterior_BT, 
     pch = 16, col = c('red', 'blue', 'green')[DFA_Field$posteriors$Group],
     xlim=c(0,1), ylim=c(0,1),
     main = 'DFA Posterior Probabilities by Group Membership',
     xlab='Posterior Probability of Being in Group CBT',
     ylab='Posterior Probability of Being in Group BT' )
legend(x=.8, y=.99, c('CBT','BT','NT'), cex=1.2, col=c('red', 'blue', 'green'), pch=16, bty='n')


# data from Green & Salkind (2008, Lesson 35)  
DFA(data = data_DFA$Green_2008,
    groups = 'job_cat', 
    variables = c('friendly','gpa','job_hist','job_test'),
    plot=TRUE,
    predictive = TRUE,
    priorprob = 'SIZES',  
    covmat_type='within', 
    CV=TRUE, 
    verbose=TRUE) 


# data from Ho (2014, Chapter 15) 
# with group_1 as numeric 
DFA(data = data_DFA$Ho_2014,
    groups = 'group_1_num',
    variables = c("fast_ris", "disresp", "sen_seek", "danger"),
    plot=TRUE,
    predictive = TRUE,
    priorprob = 'SIZES', 
    covmat_type='within', 
    CV=TRUE, 
    verbose=TRUE) 


# data from Ho (2014, Chapter 15)   
# with group_1 as a factor 
DFA(data = data_DFA$Ho_2014,
    groups = 'group_1_fac',
    variables = c("fast_ris", "disresp", "sen_seek", "danger"),
    plot=TRUE,
    predictive = TRUE,
    priorprob = 'SIZES', 
    covmat_type='within', 
    CV=TRUE, 
    verbose=TRUE) 


# data from Huberty (2006, p 45)
DFA_Huberty=DFA(data = data_DFA$Huberty_2019_p45, 
    groups = 'treatmnt_S', 
    variables = c('Y1','Y2'),
    predictive = TRUE, 
    priorprob = 'SIZES', 
    covmat_type='separate', # altho better to used 'separate' for these data
    verbose = TRUE)


# data from Huberty (2006, p 285)
DFA_Huberty=DFA(data = data_DFA$Huberty_2019_p285, 
    groups = 'Grade', 
    variables = c('counsum','gainsum','learnsum','qelib','qefac','qestacq',
                  'qeamt','qewrite','qesci'),
    predictive = TRUE, 
    priorprob = 'EQUAL', 
    covmat_type='within', 
    verbose = TRUE)


# data from Norusis (2012, Chaper 15)  
DFA_Norusis=DFA(data = data_DFA$Norusis_2012, 
    groups = 'internet', 
    variables = c('age','gender','income','kids','suburban','work','yearsed'),
    predictive = TRUE, 
    priorprob = 'EQUAL', 
    covmat_type='within', 
    verbose = TRUE)


# data from Rencher (2002, p 170)  - rootstock
DFA(data = data_DFA$Rencher_2002_root, 
    groups = 'rootstock', 
    variables = c('girth4','ext4','girth15','weight15'),
    predictive = TRUE, 
    priorprob = 'SIZES',
	covmat_type='within',     
	verbose = TRUE)


# data from Rencher (2002, p 280)  - football
DFA(data = data_DFA$Rencher_2002_football, 
    groups = 'grp', 
    variables = c('WDIM','CIRCUM','FBEYE','EYEHD','EARHD','JAW'),
    predictive = TRUE, 
    priorprob = 'SIZES',
	covmat_type='separate',     
	verbose = TRUE)


# Sherry (2006)   -  with Group as numeric 
DFA_Sherry <- DFA(data = data_DFA$Sherry_2006, 
                  groups = 'Group_num',
                  variables = c('Neuroticism','Extroversion','Openness', 
                                'Agreeableness','Conscientiousness'),
                  predictive = TRUE, 
                  priorprob = 'SIZES', 
                  covmat_type='separate', 
                  verbose = TRUE)


# Sherry (2006)   -  with Group as a factor 
DFA_Sherry <- DFA(data = data_DFA$Sherry_2006, 
                  groups = 'Group_fac',
                  variables = c('Neuroticism','Extroversion','Openness', 
                                'Agreeableness','Conscientiousness'),
                  predictive = TRUE, 
                  priorprob = 'SIZES', 
                  covmat_type='separate', 
                  verbose = TRUE)

# plots of posterior probabilities by group
# hoping to see correct separations between cases from different groups

# first, display the posterior probabilities
print(cbind(round(DFA_Sherry$posteriors[1:3],3), DFA_Sherry$posteriors[4]))

# group 1 vs 2
plot(DFA_Sherry$posteriors$posterior_1, DFA_Sherry$posteriors$posterior_2, 
     pch = 16, cex = 1, col = c('red', 'blue', 'green')[DFA_Sherry$posteriors$Group],
     xlim=c(0,1), ylim=c(0,1),
     main = 'DFA Posterior Probabilities by Original Group Memberships',
     xlab='Posterior Probability of Being in Group 1',
     ylab='Posterior Probability of Being in Group 2' )
legend(x=.8, y=.99, c('1','2','3'), cex=1.2, col=c('red', 'blue', 'green'), pch=16, bty='n')

# group 1 vs 3
plot(DFA_Sherry$posteriors$posterior_1, DFA_Sherry$posteriors$posterior_3, 
     pch = 16, col = c('red', 'blue', 'green')[DFA_Sherry$posteriors$Group],
     xlim=c(0,1), ylim=c(0,1),
     main = 'DFA Posterior Probabilities by Group Membership',
     xlab='Posterior Probability of Being in Group 1',
     ylab='Posterior Probability of Being in Group 3' )
legend(x=.8, y=.99, c('1','2','3'), cex=1.2,col=c('red', 'blue', 'green'), pch=16, bty='n')

# group 2 vs 3
plot(DFA_Sherry$posteriors$posterior_2, DFA_Sherry$posteriors$posterior_3, 
     pch = 16, col = c('red', 'blue', 'green')[DFA_Sherry$posteriors$Group],
     xlim=c(0,1), ylim=c(0,1),
     main = 'DFA Posterior Probabilities by Group Membership',
     xlab='Posterior Probability of Being in Group 2',
     ylab='Posterior Probability of Being in Group 3' )
legend(x=.8, y=.99, c('1','2','3'), cex=1.2, col=c('red', 'blue', 'green'), pch=16, bty='n')
    

# Tabachnik & Fiddel (2019, p 307, 311)   - small - with group as numeric
DFA(data = data_DFA$TabFid_2019_small, 
    groups = 'group_num', 
    variables = c('perf','info','verbexp','age'),
    predictive = TRUE, 
    priorprob = 'SIZES', 
    covmat_type='within', 
    verbose = TRUE)  
    
    
# Tabachnik & Fiddel (2019, p 307, 311)   - small - with group as a factor
DFA(data = data_DFA$TabFid_2019_small, 
    groups = 'group_fac', 
    variables = c('perf','info','verbexp','age'),
    predictive = TRUE, 
    priorprob = 'SIZES', 
    covmat_type='within', 
    verbose = TRUE)  


# Tabachnik & Fiddel (2019, p 324)   - complete  - with WORKSTAT as numeric 
DFA(data = data_DFA$TabFid_2019_complete, 
    groups = 'WORKSTAT_num', 
    variables = c('CONTROL','ATTMAR','ATTROLE','ATTHOUSE'),
    plot=TRUE,
    predictive = TRUE,
    priorprob = 'SIZES',  
    covmat_type='within', 
    CV=TRUE, 
    verbose=TRUE) 

# Tabachnik & Fiddel (2019, p 324)   - complete  -  with WORKSTAT as a factor  
DFA(data = data_DFA$TabFid_2019_complete, 
    groups = 'WORKSTAT_fac', 
    variables = c('CONTROL','ATTMAR','ATTROLE','ATTHOUSE'),
    plot=TRUE,
    predictive = TRUE,
    priorprob = 'SIZES',  
    covmat_type='within', 
    CV=TRUE, 
    verbose=TRUE) 
}
}
