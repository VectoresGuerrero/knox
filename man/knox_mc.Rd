\name{knox_mc}
\alias{knox_mc}
\docType{package}
\title{
Original Knox Test 
}
\description{
Knox Test for space-time interaction with Monte Carlo simulation
}
\usage{
knox_mc(x,y,time,sigma,tau,perm)}
\arguments{
\item{x}{Longitude - should be projected to a planar system} 
\item{y}{Latitude - should be projected to a planar system}
\item{time}{time}
\item{sigma}{A defined threshold in space distance}
\item{tau}{A defined threshold in time distance}
\item{perm}{number of replications default=999}
}
\value{
	\item{knox}{Knox statistic}
	\item{p_value}{p-value calculated from simulation}
	\item{RR}{Relative Risk - calculated by observed value (Knox test) over mean of simulated values}
}
\references{
\strong{Knox, E. (1964)}. \emph{The detection of space-time interactions.} Journal of the Royal Statistical Society. Series C (13(1), 25-30.

\strong{Tango, T. (2010)}. \emph{Statistical methods for disease clustering. Springer.}
}
\examples{
\dontrun{
}
}
