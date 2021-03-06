% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{st_link}
\alias{st_link}
\title{space-time link}
\usage{
st_link(x, y, time, ds, dt)
}
\arguments{
\item{x}{Longitude - should be projected to a planar system.}

\item{y}{Latitude - should be projected to a plannar system.}

\item{time}{time}

\item{ds}{a distance in space.}

\item{dt}{a distance in time.}
}
\value{
A data frame with \item{Xo, Yo}{ X, Y of the starting point}
\item{Xd, Yd}{ X, Y of the ending point}
}
\description{
Create space-time links between two points from a given space and time distance.
}

