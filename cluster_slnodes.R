######################################################################################
# Unify node-clusters in SpatialLinesDataFrame
#
# 2016-12-07
######################################################################################

require(sp)
require(rgdal)
require(rgeos)
require(fpc)

#' @title 
#' cluster_slnodes
#' 
#' @decription
#' Clusters start- and endpoints of SpatialLines in SpatialLinesDataFrame.
#'
#' @details 
#' Clusters neighboring start and endponts in SpatialLinesDataFrame using the DBSCAN Algorithm ( https://en.wikipedia.org/wiki/DBSCAN )
#' After clustering, each start-/endpoint of a line which is within a cluster is replaced by the clusters core, 
#' which is the start/endpoint in the data clostest to the euclidean mean of the cluster. Lines which have their start and endpoints
#' in the same cluster are deleted.
#' 
#' Only accepts LINE geometries, no MULTILINE geometries.
#' 
#'
#' @param sldf A \code{sp::SpatialLinesDataFrame} object with only LINE geometries
#' @param eps 'Reachability' distance; cutoff determining the geographic neighborhood size; in decimal degrees
#' @param MinPts Minimum number of points in a cluster.
#'
#' @return A \code{sp::SpatialLinesDataFrame} object. Data from deleted lines is deleted.
#'
#'
#'


cluster_slnodes <- function(sldf, eps, MinPts){
  # Points as start and endpoints of Lines
  points <- SpatialPoints(unique(rbind(do.call(rbind,lapply(sldf@lines, function(x){x@Lines[[1]]@coords[1,]})),
                                       do.call(rbind,lapply(sldf@lines, function(x){x@Lines[[1]]@coords[nrow(x@Lines[[1]]@coords),]})))))
  
  # Neighborhoods with DBSCAN Algorithm
  # see: https://en.wikipedia.org/wiki/DBSCAN
  DBSCAN <- dbscan(points@coords, eps = eps, MinPts = MinPts)
  
  #Identify core of each cluster
  clusters <- unique(DBSCAN$cluster[DBSCAN$cluster != 0]) # 0 is cluster of 'noise' / solitary points
  clusters <- clusters[order(clusters)]
  cores <- lapply(clusters,
                  function(x){cluster_core(points = points, cluster = x, DBSCAN = DBSCAN)})
  
  print(paste("INFO: Found", length(cores), "clusters in", length(points), "points."))
  
  # Rewire lines
  
  # a. Start and endpoints of each line
  lines.p1 <- lapply(sldf@lines, function(x){x@Lines[[1]]@coords[1,]})
  lines.p2 <- lapply(sldf@lines, function(x){x@Lines[[1]]@coords[nrow(x@Lines[[1]]@coords),]})
  # b. cluster of these points
  cluster.p1 <- lapply(lines.p1, function(x){DBSCAN$cluster[which(points@coords[,1] == x[1] & points@coords[,2] == x[2])]})
  cluster.p2 <- lapply(lines.p2, function(x){DBSCAN$cluster[which(points@coords[,1] == x[1] & points@coords[,2] == x[2])]})
  # c. rewire lines with rewire_line function
  new.lines <- lapply(c(1:length(sldf)), 
                      function(x){rewire_line(line = sldf@lines[[x]], cluster.p1[[x]], cluster.p2[[x]], cores, points)})
  
  
  # Create new SpatialLinesDataFrame
  sldf.new <- SpatialLinesDataFrame(SpatialLines(new.lines[!sapply(new.lines,is.null)]),
                                    data = sldf@data[!sapply(new.lines,is.null),])
  
  #Return
  return(sldf.new)
}

cluster_core <- function(points, cluster, DBSCAN){
  # Computes the euclidean mean of each cluster and returns the point.id of the closest point.
  clust.points <- which(DBSCAN$cluster == cluster)
  core.dist <- gDistance(points[clust.points], 
                         SpatialPoints(matrix(apply(points[clust.points]@coords, 2, mean), nrow = 1)), 
                         byid = T)
  core <- clust.points[which(core.dist == min(core.dist))[1]]
  return(core)
}

rewire_line <- function(line, cluster.1, cluster.2, cores, points){
  # Takes a line and the clusters its start and endpoints is in. Rewires the line with the cores of these clusters.
  # If both points are in the same cluster, but not solitary (cluster == 0), the line is deleted.
  if(cluster.1 == cluster.2 & cluster.1 != 0){ 
    # Line within same cluster, delete road
    return(NULL)
  }
  if(cluster.1 != 0){
    line@Lines[[1]]@coords[1,] <- points[cores[[cluster.1]],]@coords
  }
  if(cluster.2 != 0){
    line@Lines[[1]]@coords[nrow(line@Lines[[1]]@coords),] <- points[cores[[cluster.2]],]@coords
  }
  
  return(line)
  
}

