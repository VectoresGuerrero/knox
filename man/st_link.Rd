\name{st_link}
\alias{st_link}
\docType{package}
\title{
Create space-time links from given space and time distance
}
\description{
Create space-time links from given space and time distance
}
\usage{
st_link(x,y,time,ds,dt)}
\arguments{
\item{x}{Longitude - should be projected to a planar system} 
\item{y}{Latitude - should be projected to a planar system}
\item{time}{time}
\item{ds}{space distance}
\item{dt}{time distance}
}
\value{
	A resulted data frame with
	\item{Xo,Yo}{X, Y of starting point}
	\item{Xd,Yd}{X, Y of destination point}
}
\examples{
\dontrun{
}
}
