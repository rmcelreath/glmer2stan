\name{varef}
\alias{varef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract varying effect estimates from a Stan model fit by glmer2stan}
\description{
  Returns a list of mean and standard deviations of varying effects.
}
\usage{
varef( object )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{stanfit} object produced by \code{glmer2stan}}
}
\details{
  The returned list is nested, with \code{expectation} and \code{stddev} at the top level. Within each, there are separate matrices for each varying effect grouping variable. Columns are factors (named) and rows are individual groups.
}
\references{}
\author{Richard McElreath}
\seealso{\code{\link{glmer2stan}}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

