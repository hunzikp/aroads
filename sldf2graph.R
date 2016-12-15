######################################################################################
# SLDF 2 GRAPH
#
# 2016-12-07
######################################################################################

require(sp)
require(rgdal)
require(rgeos)
require(igraph)

#' @title 
#' Fast sldf2graph
#' 
#' @decription
#' Generates an iGraph graph from SpatialLinesDataFrame.
#'
#' @details 
#' Generates an iGraph graph from a SpatialLinesDataFrame object or a Shapefile containing LINE geometries. 
#' Line attribute data is transformed into edge attribute data, and an extra 'length' edge attribute with 
#' Euclidean line lengths is added.
#' 
#' Only accepts LINE geometries, no MULTILINE geometries.
#' 
#' If the \code{x} parameter is a path, its format should be \code{/foo/bar.shp}
#'
#' @param x Either a \code{sp::SpatialLinesDataFrame} object or a valid path to a shapefile containing LINE geometries.
#' @param plot.result Optional boolean triggering graph plotting.
#'
#' @return A \code{igraph::graph} object.
#'
sldf2graph <- function(x, plot.result=FALSE) {
  
  if (class(x)=="character") {
    file.name <- basename(x)
    file.name <- gsub("[.][[:graph:]]+", "", file.name)
    dir.name <- dirname(x)
    sldf <- readOGR(dir.name, file.name)
  } else if (class(x) == "SpatialLinesDataFrame") {
    sldf <- x
  } else {
    stop("x is not string or SLDF.")
  }
  
  ## Make vertex DF
  Lines.ls <- lapply(sldf@lines, function(x) x@Lines)
  Lcount <- unlist(lapply(Lines.ls, length))
  if (any(Lcount > 1)) {
    stop("MULTILINES not allowed.")
  }
  coord.ls <- lapply(unlist(Lines.ls, recursive=FALSE), function(x) x@coords)
  endcoords.ls <- lapply(coord.ls, function(x) x[c(1,nrow(x)),])
  endcoord.mat <- unique(do.call("rbind", endcoords.ls))
  vertex.df <- data.frame(vid=1:nrow(endcoord.mat), x=endcoord.mat[,1], y=endcoord.mat[,2])
  
  ## Make edgelist DF
  endcoords.el.ls <- lapply(endcoords.ls, function(x) matrix(c(x[1,], x[2,]), ncol=4))
  el.mat <- do.call("rbind", endcoords.el.ls)
  el.df <- as.data.frame(el.mat)
  names(el.df) <- c('from.x', 'from.y', 'to.x', 'to.y')
  el.df$order <- 1:nrow(el.df)
  fromvertex.df <- vertex.df
  names(fromvertex.df) <- c('from.vid', 'from.x', 'from.y')
  el.df <- merge(el.df, fromvertex.df, all.x=TRUE, all.y=FALSE)
  tovertex <- vertex.df
  names(tovertex) <- c('to.vid', 'to.x', 'to.y')
  el.df <- merge(el.df, tovertex, all.x=TRUE, all.y=FALSE)
  el.df <- el.df[order(el.df$order),]
  el.df <- el.df[,c('from.vid', 'to.vid', 'order')]
  
  ## Add attribute data, euclidean length
  el.df <- cbind(el.df, sldf@data)
  if (!('length' %in% names(el.df))){
    suppressWarnings(el.df$length <- gLength(sldf, byid=TRUE))
  }
  
  ## Make graph
  graph <- graph_from_data_frame(el.df, FALSE, vertex.df)
  
  ## Optional plotting
  if (plot.result) {
    asp <- (max(endcoord.mat[,2])-min(endcoord.mat[,2]))/(max(endcoord.mat[,1])-min(endcoord.mat[,1]))
    plot(graph, vertex.size=0.25, vertex.label=NA, layout=as.matrix(vertex.df[,2:3]), asp=asp)
  }

  return(graph)
}