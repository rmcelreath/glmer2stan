\name{stanpredict}
\alias{stanpredict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute predictions from a model fit with glmer2stan}
\description{
  This function computes predictions and prediction intervals from a model fit by \code{glmer2stan}. Predictions use the entire posterior for inference. Uncertainty due to both posterior probability and sampling distribution are included.
}
\usage{
stanpredict( stanfit , data , vary_prefix="vary_" , fixed_prefix="beta_" , 
    probs=c(0.025,0.975) , nsims=1e4 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{\code{stanfit} object, produced with \code{\link{glmer2stan}}}
  \item{data}{data frame of cases to use for prediction}
  \item{fixed_prefix}{Identifying prefix for fixed effects, in the samples}
  \item{vary_prefix}{Identifying prefix for varying effects, in the samples}
  \item{probs}{Prediction quantiles to report (default is 95 percent interval)}
  \item{nsims}{Number of simulations, for generating sampling distributions}
}
\details{
  This function uses the posterior samples from a \code{stanfit} object, fit with \code{glmer2stan}, to simulate model predictions.
  
  The returned value is a list with a named entry for each outcome variable in the model. Under each outcome variable is a list with: (1) \code{mu}, expectation of mean; (2) \code{mu.ci}, quantiles for the mean; and (3) \code{obs.ci}, quantiles for the predicted sampling distribution. 
  
  When passing new data to \code{stanpredict}, you can use a list with different length vectors. Short vectors will be grown to match the length of the longest vector, using repetition. For cluster variables, using a value of zero will be interpreted as omitting varying effects from predictions, generating predictions for average cases in the population. Otherwise, cluster variable values should be integers.
}
\references{}
\author{Richard McElreath}
\seealso{\code{\link{glmer2stan}}}
\examples{
# simple linear regression
data(cars)
m <- glmer2stan( dist ~ speed , data=cars )
speed.seq <- seq( from=min(cars$speed) , to=max(cars$speed) , length.out=20 )
pred.dist <- stanpredict( m , data=list( speed=speed.seq ) )$dist
plot( dist ~ speed , cars )
lines( speed.seq , pred.dist$mu )
lines( speed.seq , pred.dist$mu.ci[1,] , lty=2 )
lines( speed.seq , pred.dist$mu.ci[2,] , lty=2 )

# binomial example
data(UCBadmit)
m <- glmer2stan( cbind(admit,reject) ~ (1|dept) + male , 
    data=UCBadmit , family="binomial" )
pred.admit <- stanpredict( m , data=UCBadmit )$admit
prop.admit <- UCBadmit$admit / UCBadmit$applications
plot( prop.admit , ylim=c(0,1) , pch=ifelse(UCBadmit$male==1,2,1) ,
    ylab="probability admit" , xlab="case" )
lines( 1:12 , pred.admit$mu )
lines( 1:12 , pred.admit$mu.ci[1,] , lty=2 )
lines( 1:12 , pred.admit$mu.ci[2,] , lty=2 )
lines( 1:12 , pred.admit$obs.ci[1,] , lty=3 )
lines( 1:12 , pred.admit$obs.ci[2,] , lty=3 )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

