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

# # ## Test
# N <- 20
# sl.ls <- vector('list', N)
# for (i in 1:N) {
#   coords <- matrix(runif(4), 2, 2)
#   sl.ls[[i]] <- SpatialLines(list(Lines(slinelist = list(Line(coords = coords)), ID = i)))
# }
# node.sl <- perform_noding(sl.ls)
# plot(node.sl, col=sample(colors(), length(node.sl)))


#' Node SpatialLineDataFrame Objects
#'
#' Nodes a list of SpatialLineDataFrame objects, i.e. inserts line breaks at intersections.
#' New noded lines inherit attribute data of parent lines.
#' 
#' Notes: 
#' \itemize{
#' \item Requires GEOS version 3.4.0 or newer.
#' \item Resets Line IDs.
#' \item All data slots in \code{x} must have the same columns.
#' }
#' 
#' @param x A list of SpatialLineDataFrame objects. 
#'
#' @return A SpatialLineDataFrame object.
#'
#' @export
perform_noding_sldf <- function(x, sample_n=5) {
  require(sp)
  require(rgeos)

  # Combine data
  data <- do.call(rbind, lapply(x, function(x) x@data))
  
  # Get list of spatial lines
  sl.ls <- lapply(x, {
    function(x) {
      SpatialLines(x@lines)
    }
  })
  
  # Replace feature IDs to enable rbind
  cntr <- 1
  sl.ls <- lapply(sl.ls, {
    function(sl) {
      sl <- spChFIDs(sl, as.character(cntr:(cntr+length(sl)-1)))
      cntr <<- cntr + length(sl)
      sl
    }
  })
  
  # Combine
  sl <- do.call(rbind, sl.ls)
  
  # Node
  nd.sl <- gNode(sl)
  
  # Unpack nodified lines
  nd.lines <- nd.sl@lines
  nd.line.ls <- unlist(lapply(nd.lines, function(x) x@Lines))
  cntr <- 0
  up.lines <- lapply(nd.line.ls, {
    function(x) {
      cntr <<- cntr + 1
      Lines(slinelist = list(x), ID=cntr)
    }
  })
  nd.sl <- SpatialLines(up.lines)
  
  # Get reference ID for each new nodified line
  in_ref_ls <- lapply(nd.sl@lines, {
    function(x) {
      this_sp <- spsample(x, n=sample_n, type="regular")
      dmat <- gDistance(this_sp, sl, byid=TRUE)
      ref_ids <- apply(dmat, 2, which.min)
      unique_ref_ids <- unique(ref_ids)
      modal_ref_id <- unique_ref_ids[which.max(tabulate(match(ref_ids, unique_ref_ids)))]
      return(modal_ref_id)
    }
  })
  
  # Get attribute data of nodified lines
  nd_data_ls <- lapply(in_ref_ls, {
    function(x) {
      data[x,,FALSE]
    }
  })
  nd_data <- do.call(rbind, nd_data_ls)
  
  # Make output sldf
  out.sldf <- SpatialLinesDataFrame(nd.sl, nd_data, FALSE)
  
  return(out.sldf)
} 



## Test
# make_random_sl <- function(N, j) {
#   sl.ls <- vector('list', N)
#   for (i in 1:N) {
#     coords <- matrix(runif(4), 2, 2)
#     sl.ls[[i]] <- SpatialLines(list(Lines(slinelist = list(Line(coords = coords)), ID = paste0(i, '_', j))))
#   }
#   sl <- do.call(rbind, sl.ls)
#   return(sl)
# }
# 
# M <- 2
# cntr <- 1
# sldf.ls <- vector('list', M)
# for (j in 1:M) {
#   N <- sample(1:100, 1)
#   sl <- make_random_sl(N, j)
#   df <- data.frame(id=cntr:(cntr+N-1))
#   rownames(df) <- getSLLinesIDSlots(sl)
#   sldf.ls[[j]] <- SpatialLinesDataFrame(sl, df, FALSE)
#   cntr <- cntr + N
# }
# 
# system.time(node.sldf <- perform_noding_sldf(sldf.ls))
# sldf <- do.call(rbind, sldf.ls)
# shortest_id <- node.sldf$id[gLength(node.sldf, T) == min(gLength(node.sldf, T))]
# plot(sldf[sldf$id==shortest_id,])
# plot(node.sldf[node.sldf$id==shortest_id,],col='red', add=TRUE)



