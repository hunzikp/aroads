#' Node SpatialLine Objects
#'
#' Nodes a list of SpatialLine objects, i.e. inserts line breaks at intersections.
#' 
#' Notes: 
#' \itemize{
#' \item Requires GEOS version 3.4.0 or newer.
#' \item Resets Line IDs.
#' }
#' 
#' @param x A list of SpatialLine objects. 
#'
#' @return A SpatialLines object.
#'
#' @export
perform_noding <- function(x) {
  require(sp)
  require(rgeos)
  
  # Replace feature IDs to enable rbind
  cntr <- 1
  x <- lapply(x, {
    function(sl) {
      sl <- spChFIDs(sl, as.character(cntr:(cntr+length(sl)-1)))
      cntr <<- cntr + length(sl)
      sl
    }
  })

  # Combine
  sl <- do.call(rbind, x)
  
  # Node
  nd.sl <- gNode(sl)
  
  # Unpack
  nd.lines <- nd.sl@lines
  nd.line.ls <- unlist(lapply(nd.lines, function(x) x@Lines))
  cntr <- 0
  out.lines <- lapply(nd.line.ls, {
    function(x) {
      cntr <<- cntr + 1
      Lines(slinelist = list(x), ID=cntr)
    }
  })
  out.sl <- SpatialLines(out.lines)
  return(out.sl)
} 

# ## Test
# N <- 20
# sl.ls <- vector('list', N)
# for (i in 1:N) {
#   coords <- matrix(runif(4), 2, 2)
#   sl.ls[[i]] <- SpatialLines(list(Lines(slinelist = list(Line(coords = coords)), ID = i)))
# }
# node.sl <- perform_noding(sl.ls)
# plot(node.sl, col=sample(colors(), length(node.sl)))
