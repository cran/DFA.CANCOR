\name{GROUP.PROFILES}
\alias{GROUP.PROFILES}
\title{Group Profile Plots}
\description{Produces profile plots of group means for one or more continuous outcome variables.}
\usage{GROUP.PROFILES(data, groups, variables,  
        plot_type ='bar', bar_type = 'all',
        rescale='standardize',
        CI_level= 95, ylim=NULL,
        plot_save = FALSE, plot_save_type = 'png', plot_title = NULL,
        cols_user = NULL,
        verbose=TRUE)}
\arguments{
  \item{data}{
  \code{}A dataframe where the rows are cases and the columns are the variables.}

  \item{groups}{
  \code{}The name of the groups variable in data, e.g., groups = 'Group'.}
  	    
  \item{variables}{
  \code{}The name of the dependent (outcome) variable(s) in data, 
  e.g., variables = c('esteem','anxiety').}

  \item{plot_type}{
  \code{}The options are 'bar' for bar plot, or 'profile' for a lines profile plot.}

  \item{bar_type}{
  \code{}When plot_type = 'bar', the options for bar_type are 'all', for placing 
  the bar plots for all variables in one plot, or 'separate', for placing the bar 
  plots for the variables in separate plots.}

  \item{rescale}{
  \code{}(optional) Should the variables be rescaled into a common metric? The options
  are 'no' (the default), or 'standardize'.}

  \item{CI_level}{
  \code{}(optional) The confidence interval for the input, if provided (in whole numbers). 
  The default is 95.}

  \item{ylim}{
  \code{}(optional) Limits for the y-axis, e.g., ylim = c(0, 5). Not implemented when multiple
  bar plots are requested.}

  \item{plot_save}{
  \code{}Should a plot be saved to disk? TRUE or FALSE (the default).}

  \item{plot_save_type}{
  \code{}The output format if plot_save = TRUE. The options are 'bitmap', 'tiff', 
         'png' (the default), 'jpeg', and 'bmp'.}

  \item{plot_title}{
  \code{}(optional) A title for the plot. 
  \cr Example: title = 'Group Profiles'}

  \item{cols_user}{
  \code{}A vector of colours for the groups. 
  \cr If NULL, the default colours are selected, in order, from this vector: 
   cols_user <- c('blue', 'red', 'cyan2', 'darkviolet', 'chartreuse1', 'yellow',
              'burlywood3','darkseagreen1', 
              'mediumvioletred', 'darkgreen','bisque','cyan3', 'deeppink4').}

  \item{verbose}{
  \code{}(optional) Should detailed results be displayed in console? \cr The options 
  are: TRUE (default) or FALSE.}
}
\details{The continuous 'variables' can be rescaled (standardized) into the same 
  metric to facilitate interpretation when the variables that are in
  different metrics are placed on one plot.

  When plot_type = 'bar' and bar_type = 'separate', a maximum of four plots will be 
  produced, for the first four 'variables'.
}
\value{If verbose = TRUE, the displayed output includes the means, standard 
  deviations, Ns, and confidence intervals for the groups on the variables.
}
\author{Brian P. O'Connor}
\examples{
GROUP.PROFILES(data = data_DFA$Ho_2014,
               groups = 'group_1_fac',
               variables = c("fast_ris", "disresp", "sen_seek", "danger"),
               rescale= 'data', 
               plot_type ='bar', 
               bar_type = 'separate')
\donttest{
#first run DFA
DFA_output <- DFA(data = data_DFA$Field_2012, 
                  groups = 'Group', 
                  variables = c('Actions','Thoughts'),
                  predictive = TRUE, 
                  priorprob = 'EQUAL',   
                  covmat_type='separate',
                  verbose = TRUE)

# then produce a profile plot of the group centroids on the discriminant functions
GROUP.PROFILES(data = DFA_output$dfa_scores,
               groups = 'group',
               variables = c('Function.1','Function.2'),
               rescale= 'no',
               plot_type ='profile',
               bar_type = 'separate')
}
}
