\name{map2stan}
\alias{map2stan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Define Stan model using likelihood/prior formula notation}
\description{
  Builds a Stan model with a single likelihood and custom priors for any of the parameters, using formula notation.
}
\usage{
map2stan( flist , data , start , fit=TRUE , chains=1 , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{flist}{Model formula or list of formulas, using \code{\link{map}} notation with explicit functions. See examples.}
  \item{data}{Data frame or list}
  \item{start}{A list of named parameters with starting values}
  \item{fit}{If \code{TRUE}, samples from the model using \code{\link{stan}}}
  \item{chains}{Number of independent chains to use in sampling}
  \item{...}{additional parameters to pass to \code{\link{stan}}}
}
\details{
  This function builds a simple fixed effect model from a list of function-explicit formulas of the form \code{y ~ dfoo(theta)}, where \code{y} is an outcome variable or parameter, \code{dfoo} is a density function such as \code{dnorm}, and \code{theta} is a parameter of the density. For example, to estimate the mean and standard deviation of a Gaussian distribution from a set of samples \code{y}, the formula specifying the likelihood might be:
  
  \code{y ~ dnorm(mu,sigma)}
  
  See the examples for the correspoding \code{start} list. Priors can be added simply by adding more formulas to the list, but now with the 'outcomes' being parameters:
  
  \code{y ~ dnorm(mu,sigma)}\cr
  \code{mu ~ dnorm(0,10)}\cr
  \code{sigma ~ dcauchy(0,1)}
  
  This interface allows for simple penalized likelihood models to be specified. For example, by placing an informative prior on a regression coefficient, the resulting posterior will shrink the coefficient towards to mean of the prior.
  
  Identical priors for multiple parameters can be specified by using the \code{c} function on the left side of the formula. For example, the same prior can be assigned to parameters \code{beta1} and \code{beta2} with:
  
  \code{c(beta1,beta2) ~ dnorm(0,10)}
  
  \code{map2stan} interprets this as:
  
  \code{beta1 ~ dnorm(0,10)}\cr
  \code{beta2 ~ dnorm(0,10)}
}
\references{
}
\author{Richard McElreath}
\seealso{\code{\link{map}},\code{\link{mle2}},\code{\link{stan}}}
\examples{
data(cars)
flist <- list(
    dist ~ dnorm( mean=a+b*speed , sd=sigma ) ,
    c(a,b) ~ dnorm(0,10) , 
    sigma ~ dcauchy(0,1)
)

# see Stan code
stan_code <- map2stan( flist , start=list(a=40,b=0.1,sigma=20) , data=cars , fit=FALSE )
cat( stan_code[[1]] )

# sample from the Stan model
stan_samples <- map2stan( flist , start=list(a=40,b=0.1,sigma=20) , data=cars , fit=TRUE , chains=4 )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

