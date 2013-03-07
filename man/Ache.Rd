\name{Ache}
\alias{Ache}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Human forager hunting returns data}
\description{
  Hunting returns of individual Ache men, 1981 to 2007. Each row is an individual trip by an individual hunter.
}
\usage{
data(Ache)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
}
\format{
    \enumerate{
        \item month : Month of trip
        \item day : Day of trip
        \item dayid : Unique day index
        \item year : Year of record
        \item date : Julian date of trip
        \item date.s : Julian date, standardized
        \item hunter.id : Identifier of individual man
        \item age : Man's age at time of record
        \item age.s : Age, standardized
        \item kg.meat : Kilograms of meat returned from hunt
        \item iszero : Indicator for zero value kg.meat
        \item nonzero : Positive kg.meat values
        \item hours : Duration in hours of hunting trip
        \item datatype : 1 if duration of trip known, 3 otherwise
    }
}
\value{
}
\references{Hill and Kintigh. 2009. Current Anthropology 50:369-377.}
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

