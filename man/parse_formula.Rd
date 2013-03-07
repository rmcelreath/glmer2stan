\name{parse_formula}
\alias{parse_formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parse a mixed model formula into components}
\description{
  Takes a formula and data frame and parses them into lists of fixed and varying effects, with corresponding model matrices. Used internally by \code{\link{glmer2stan}}.
}
\usage{
parse_formula( formula , data )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A formula with varying effects}
  \item{data}{A data frame or list}
}
\details{
  The result of \code{parse_formula} is a list containing: (1) \code{y}, the outcome variable; (2) \code{yname}, the character name of the outcome; (3) \code{fixef}, a vector of fixed effect names; (4) \code{ranef}, a named list of varying (aka random) effect grouping variables, each containing a vector of factor names; (5) \code{dat}, a data frame containing all of the needed variables in the fixed and varying effects, including the grouping variables.
}
\references{}
\author{Richard McElreath}
\seealso{\code{\link{glmer2stan}}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

