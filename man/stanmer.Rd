\name{stanmer}
\alias{stanmer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarize a Stan model fit by glmer2stan}
\description{
  Displays mean and standard deviations of fixed effects samples and variance-covariance samples. Returned value is a structured list of mean and standard deviations of all parameters.
}
\usage{
stanmer( fit , digits=2 , probs=c(0.025,0.975) , fixed_prefix="beta_" , vary_prefix="vary_" )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{\code{stanfit} object, produced with \code{\link{glmer2stan}}}
  \item{probs}{quantiles to display in summary}
  \item{digits}{Number of digits to display in summaries}
  \item{fixed_prefix}{Identifying prefix for fixed effects, in the samples}
  \item{vary_prefix}{Identifying prefix for varying effects, in the samples}
}
\details{
  This command displays summary statistics in a form resembling the output for a regular \code{glmer} model fit. Averages and standard deviations are shown for all fixed effects, and variance components are displayed as standard deviations and correlations. 
  
  The variance component display or each grouping variable shows standard deviation on the diagonal and correlations in the lower triangle. 
  
  Each grouping variable also displays a measure of imbalance in sampling across groups within it. This measure is simply the gini coeficient of the distribution of the number of samples in each unique group. If every group (individual, household, county, etc.) has the same number of cases in the data, then the gini is zero, indicating no imbalance in sampling. As imbalance increases, the gini approaches 1.
  
  The returned value is a list of mean and standard deviations of all parameters in the Stan samples, organized as lists of arrays.
}
\references{}
\author{Richard McElreath}
\seealso{\code{\link{glmer2stan}}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

