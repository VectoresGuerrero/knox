\name{iknox_mc}
\alias{iknox_mc}
\docType{package}
\title{
Interval Knox Test 
}
\description{
Interval Knox Test for space-time interaction with Monte Carlo simulation
}
\usage{
knox_mc(x,y,time,sigma,tau,perm)}
\arguments{
\item{x}{Longtitude - should be projected to a planar system} 
\item{y}{Latitude - should be projected to a planar system}
\item{time}{time (day)}
\item{sigma}{A defined threshold in space distance (m,km...)}
\item{tau}{An interval between case i and j}
\item{perm}{number of replications default=999}
}
\value{
	\item{IKT}{Interval Knox statistic Test}
	\item{Z_score}{Z score}
	\item{p_value}{p-value calculated from simulation}
	\item{RR}{Relative Risk - calculated by observed value (Knox test) over mean of simulated values}
}
\references{\strong{Aldstadt, J. (2007)}. \emph{An incremental Knox test for the determination of the serial interval between successive cases of an infectious disease}. Stochastic Environmental Research and Risk Assessment, 21(5), 487-500. doi:10.1007/s00477-007-0132-3}
\examples{
\dontrun{
}
}
