#' @title Writting shapefile of space-time links
#' @description Writing shapefile from a data frame of space-time links
#' @author Thanh Le Viet, \email{lethanhx2k@@gmail.com}
#' @param data A data frame holds at least 3 columns: longtitude, latitude and 
#'   time.
#' @param x Longitude, should be projected to a planar system.
#' @param y Latitude, should be projected to a planar system.
#' @param time time column.
#' @param ds  a cut-off distance in space.
#' @param dt  a cut-off distance in time.
#' @param w default is TRUE, to write out a shapefile.
#' @param file.out default is "st_links", shapfile name to be written.
#' @param ... parameters of the plot function of a line, i.e col, width, etc.
#' @return A plot of space-time links in lines. By default, \code{w} is set 
#'   TRUE, a shapefile will be written out as well.
#' @seealso \code{\link{st_link}}   
#' @export
#' @importFrom sp SpatialLines SpatialLinesDataFrame
#' @importFrom rgdal writeOGR
#'
writeSHP <- function(data, x, y, time, ds, dt, w = TRUE, file.out = "st_links", ...){
  require('sp')
  require('rgdal')
  st_link <- with(data,st_link(x = x, y = y,time = time, ds = ds, dt = dt))
  l <- list()
  for (i in 1:nrow(st_link)) {
    l[[i]] <- with(st_link[i,],Lines(Line(cbind(c(Xo,Xd),c(Yo,Yd))),ID = i))
  }
  id <- data.frame(id = c(1:nrow(st_link)))
  shape <- SpatialLines(l)
  shape.line.df <- SpatialLinesDataFrame(shape,id,match.ID = TRUE)
  if (w) {
    writeOGR(shape.line.df,dsn = ".",file.out,driver = "ESRI Shapefile")
  }
plot(shape,...)
}