\name{writeSHP}
\alias{writeSHP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Export space-time link to shapefile
}
\description{
Writing space-time link from given distance and time to a shapefile inside working directory
}
\usage{
writeSHP(x,y,time,ds,dt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Longitude - should be projected into a planar system}
  \item{y}{Latitude - should be projected into a planar system}
  \item{time}{time}
  \item{ds}{space distance}
  \item{dt}{time distance}
  \item{w}{boolean value, if True: a shapefile will be written out, otherwise it will be plot in R. default = TRUE}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~st_link }
\keyword{ ~shape }% __ONLY ONE__ keyword per line
