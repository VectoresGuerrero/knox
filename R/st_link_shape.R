writeSHP <- function(data,x,y,time,ds,dt,w=TRUE,...){
  require(sp)
  require(rgdal)
  st_link <- with(data,st_link(x=x,y=y,time=time,ds=ds,dt=dt))
  l <- list()
  for(i in 1:nrow(st_link)){
    l[[i]] <- with(st_link[i,],Lines(Line(cbind(c(Xo,Xd),c(Yo,Yd))),ID=i))
  }
  id <- data.frame(id=c(1:nrow(st_link)))
  shape <- SpatialLines(l)
  shape.line.df <- SpatialLinesDataFrame(shape,id,match.ID = TRUE)
  if (w){
    writeOGR(shape.line.df,dsn=".","cat_ba_st_links",driver = "ESRI Shapefile")  
  }
  plot(shape,...)  
}