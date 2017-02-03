

######## Spatial Lines Network to gridded network
require(sp)
require(raster)
require(igraph)


# FUNCTIONS #####################
explode_multilines <- function(multiline){
  lines <- lapply(c(1:length(multiline@lines[[1]]@Lines)), function(x){Lines(list(multiline@lines[[1]]@Lines[[x]]), 
                                                                             ID = paste0(multiline@lines[[1]]@ID,".",x))})
  data <- multiline@data[rep(1, length(lines)),]
  row.names(data) <- paste0(multiline@lines[[1]]@ID,".",c(1:length(lines)))
  lines.spldf <- SpatialLinesDataFrame(SpatialLines(lines), data)
  return(lines.spldf)
}

road_edges <- function(x, id){
  # Retrieves road edges from ordered(!) list of vertex ids
  if(length(x) < 2){
    return(NULL)
  } else if(length(x) == 2){
    res <- data.frame(matrix(x[order(x)],ncol = 2))
    colnames(res) <- c("Vert1","Vert2")
    res$id <- id
    return(res)
  } else {
    res <- data.frame(do.call(rbind,lapply(c(1:(length(x)-1)), function(y){x[y:(y+1)][order(x[y:(y+1)])]})))
    colnames(res) <- c("Vert1","Vert2")
    res$id <- id
    return(res)
  }
}

find.adj <- function(ncol, nrow){
  # Makes adjacency 2*X matrix listing all paris of direct neighbors. Same as raster::adjacent(.., directions = 4), but faster.
  adj.mat <- cbind(c(c(1:(ncol*nrow))[-seq(ncol, ncol*nrow, by = ncol)], c(1:(ncol*nrow))),
                   c(c(1:(ncol*nrow))[-seq(ncol, ncol*nrow, by = ncol)]+1, c(1:(ncol*nrow))+ncol))
  return(adj.mat[adj.mat[,2]>0 & adj.mat[,2]<= ncol*nrow,])
}

grid_spnetwork <- function(raster, spldf.ls = list(), fun.ls = list(), obstacles.ls = list()){
  # Produces gridded network from spatial lines dataframes
  # 
  # Args:
  #   raster: a raster layer of raster stack, of the extent of the output network. 
  #     Centroids are transformed into verteces, raster values into vertex attributes
  #   spldf: names list of spatial line data frames. Lines will be mapped onto the gridded network. 
  #     Variables become edge-attributes named with the prefix "name." of the spatial lines data frame.
  #   obstacles.ls: named list of SpatialLinesDataFrames / SpatialPolygonsDataFrame. Objects are intersected with 
  #     edges, producing dummy-variables for whether an edge intersects with a particular obstacle or not. Named as name of data set. 
  #   fun.list: named list of lists (one for each entry in spldf), with sub.names being variables and entries functions for aggregation
  
  
  #Create empty raster of extension of roads
  if(nlayers(raster) == 1){
    road.rs <- raster
  } else {
    road.rs <- raster[[1]]
  }
  road.rs[road.rs != 0] <- 0
  road.rs[is.na(road.rs)] <- 0
  
  # Clean spatial lines data if exhibits multilines
  for(s in c(1:length(spldf.ls))){
    multiline <- unlist(lapply(spldf.ls[[s]]@lines, function(x){length(x@Lines) > 1}))
    if(any(multiline)){
      new.lines <- lapply(which(multiline),function(x){explode_multilines(spldf.ls[[s]][x,])})
      new.lines <- do.call(rbind, new.lines)
      spldf.ls[[s]] <- rbind(spldf.ls[[s]][!multiline,], new.lines)
    } 
  }

  
  # Get adjacency matrix
  adj.mat <- find.adj(ncol(road.rs), nrow(road.rs))
  adj.df <- data.frame(adj.mat)
  colnames(adj.df) <- c("Vert1","Vert2")
  
  # Make graph
  spl.nw <- graph_from_edgelist(adj.mat, directed = F)
  
  # Incorporate spatial lines data - loop over spldf.ls
  if(length(spldf.ls)>0){
    for(s in c(1:length(spldf.ls))){
      
      # get this entry
      spldf <- spldf.ls[[s]]
      this.name <- names(spldf.ls)[s]
      this.fun.ls <- fun.ls[[s]]
      print(paste("Integrate", this.name))
      
      # Get line paths and edges
      spl.paths <- raster::extract(road.rs, spldf, cellnumbers = T, along = T)
      spl.edges <- lapply(c(1:length(spl.paths)), function(x){road_edges(spl.paths[[x]][,1], id = x)})
      spl.edges <- do.call(rbind, spl.edges)
      
      # Join with spatial line data
      spl.edges.df <- cbind(spl.edges, spldf@data[spl.edges$id,])
      
      # Aggregate edge data
      spl.uni <- unique(spl.edges.df[spl.edges.df$Vert1 != spl.edges.df$Vert2,c("Vert1","Vert2")])
      spl.uni$update <- 1
      colnames(spl.uni)[colnames(spl.uni) == "update"] <- this.name
      for(var in names(this.fun.ls)){
        var.agg <- try(aggregate.data.frame(spl.edges.df[,var], spl.edges.df[,c("Vert1","Vert2")], FUN = this.fun.ls[[var]]))
        if(class(var.agg) == "try-error"){
          warning(paste("Can't aggregate", var,". Skipping."))
        } else {
          colnames(var.agg) <- c("Vert1","Vert2", var)
          spl.uni <- join(spl.uni, var.agg, by = c("Vert1","Vert2"), type = "left")
        }
      }
      
      # Join with adjacency data.frame
      spl.uni.edge <- join(adj.df, spl.uni, type = "left", by = c("Vert1","Vert2"))
      spl.uni.edge[,this.name][is.na(spl.uni.edge[,this.name])] <- 0
      
      # Set edge attributes
      for(var in colnames(spl.uni)[!colnames(spl.uni) %in% c("Vert1","Vert2")]){
        edge_attr(spl.nw, paste0(this.name,".",var), E(spl.nw)) <- spl.uni.edge[,var]
      }
    }
  }

  # Save raster data to verteces
  print("Integrate vertex data")
  vertex.data <- data.frame(rasterToPoints(stack(raster,road.rs), spatial=F))
  #   stacking is necessary here to prevent rasterToPoints from omitting cells with missing data.
  vertex.data <- vertex.data[,-ncol(vertex.data)]
  for(var in colnames(vertex.data)){
    vertex_attr(spl.nw, var) <- vertex.data[,var]
  }
  
  # Integrate obstacles
  if(length(obstacles.ls) > 0){
    # Crop to extent of raster
    print("Integrate obstacles")
    obstacles.ls <- lapply(obstacles.ls, function(x){raster::crop(x, extent(raster))})
    
    # Make spatial lines from edges
    edge.vert <- get.edges(spl.nw, E(spl.nw))
    edges.sl <- lapply(c(1:nrow(edge.vert)), function(x){Lines(list(Line(vertex.data[edge.vert[x,],c("x","y")])), ID = as.character(x))})
    edges.sl <- SpatialLines(edges.sl)
    
    # Check for intersections and create dummy for it as edge attribute
    count = 0
    for(o in obstacles.ls){
      count = count + 1
      if(!is.null(o) & length(o)>0){
        obs.intersects <- gIntersects(edges.sl, o, byid = T)
        edge_attr(spl.nw, names(obstacles.ls)[count]) <- as.numeric(apply(obs.intersects, 2, any))
      } else {
        edge_attr(spl.nw, names(obstacles.ls)[count]) <- as.numeric(FALSE)
      }
    }
  }

  # Return data
  return(spl.nw)
  
}

