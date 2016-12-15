

#### ROAD NETOWKR EXPLORATORY ANALYSIS

# Carl Müller-Crepon w/ Philipp Hunziker, Matthew Simonson, ETHZ, Northestern


#############################################

#Load functions and packages
library(gtools)
library(rgeos)
library(geosphere)
library(rgdal)
library(cshapes)
library(maptools)
library(igraph)
source("sldf2graph.R")

# Load data

road.spdf <- readShapeLines("~/Data/Roads/Michelin/182_1941_referenced")
plot(road.spdf)


#Make network
roads.nw <- sldf2graph(road.spdf)


# Generate simplified network with nodes at intersections & endpoints only
simplify_network <- function(graph){
  #Find all vertices with only two edges
  edge.num <- unlist(lapply(V(graph), function(x){length(incident(graph, v = x, mode = c("all")))}))
  vertices.bridge <- which(edge.num == 2)
  
  #Sequentially (!) delete old edges and replace them with a bridging one
  for(v in vertices.bridge){
    these.neighbors <- neighbors(graph, v, mode =  "all")
    old.edges <- E(graph)[from(v)]
    #Add new edge
    graph <- add_edges(graph, these.neighbors,
                       attr = list(weigt = sum(edge_attr(graph, "weight", old.edges))))
    #Delede old edges
    graph <- delete.edges(graph, E(graph)[from(v)])
  }
  #Delete now unconnected vertices
  graph <- delete.vertices(graph, c(vertices.bridge, which(edge.num == 0)))
  #Return
  return(graph)
}

roads.simple.nw <- simplify_network(graph = roads.nw)


# Make simplified spatial roads data.frame
coords <- cbind(vertex_attr(roads.simple.nw, "x", index = V(roads.simple.nw)), 
                vertex_attr(roads.simple.nw, "y", index = V(roads.simple.nw)))
edges <- get.edges(roads.simple.nw, E(roads.simple.nw))
lines <- lapply(c(1:nrow(edges)), function(x){Lines(list(Line(coords[edges[x,],])), as.character(x))})
roads.simple.spdf <- SpatialLines(lines) 

plot(road.spdf, col = "red")
lines(roads.simple.spdf, col = "blue")

############## ANALYSIS ###################

# Set an analysis network

analysis.nw <- roads.simple.nw 
analysis.spdf <- roads.simple.spdf

#Create weights and correct 0-weights to small number..
analysis.nw <- set_edge_attr(analysis.nw, "weight", index = E(analysis.nw), 
                             unlist(lapply(edge_attr(analysis.nw, "length"),function(x){max(x, 0.0000001)})))

# Correlation age - edge betweenness centrality #######################
# 1. Calculate edge betweenness centrality

edge.bc <- edge_betweenness(analysis.nw, e = E(analysis.nw), directed = F) # , weights = rep(1, length(E(analysis.nw)))
plot(density(log(1+edge.bc)))

# 2. We can now correlate this with the age of an edge...

# 3. Plot in space
col.pal <- colorRampPalette(c("yellow","orange","red"),space = "Lab")(256)
col.vec <- col.pal[as.numeric(cut(log(1+edge.bc/1000),breaks = 256))]

plot(analysis.spdf, col = col.vec)

# Vertice betweenness ######################################
edge.num <- unlist(lapply(V(analysis.nw), function(x){length(incident(analysis.nw, v = x, mode = c("all")))}))
plot(density(edge.num))
vertice.bc <- betweenness(analysis.nw, v = V(analysis.nw), directed = F)

col.pal <- colorRampPalette(c("yellow","orange","red"),space = "Lab")(256)
col.vec <- col.pal[as.numeric(cut(vertice.bc,breaks = 256))]
plot(analysis.spdf)
points(cbind(vertex_attr(analysis.nw, "x"), 
             vertex_attr(analysis.nw, "y")), pch = 20, col = col.vec)

# Exploratory vs. densifying roads #########################

# This doesn't really square with the results of Strano et al. 
# In particular, I get negative difference, when I remove deadends, which should be "exploratory" routes in their sense. 
# I don't know, how the mean of the betweenness centrality can decrease, when a deadend edge (and the attibuted vertices) is removed 
# or bridged (to avoid 2-edge vertices). 
# Maybe I understand their formula wrong, and they actually only use the absolute difference - which wouldn't make too much sense I think...

# 1. (function) Mean edge betweenness (formulas 4 & 5; Strano et al. 2012)
mean_edge_bc <- function(graph){
  # Simplify: bridge and delete vertices with 2 (or 0) edges only
  graph <- simplify_network(graph)
  #Calculate mean edge betweenness centrality
  edge.bc <- edge_betweenness(graph, e = E(graph), directed = F) #, weights = rep(1, length(E(graph)))
  mean.edge.bc <- sum(edge.bc)/((length(V(graph))-1)*(length(V(graph))-2))
  return(mean.edge.bc)
}
# 2. (function) Difference in mean edge betweenness (formulas 6; Strano et al. 2012)
diff_edge_bc <- function(edge, graph, mean.edge.bc){
  # Calculates difference in the mean edge betweenness in the graph, if edge edge gets deleted. 
  # Takes full mean.edge.bc as argument to speed things up
  red.edge.bc <- mean_edge_bc(graph = delete.edges(graph,edge))
  diff.edge.bc <- (mean.edge.bc - red.edge.bc)/mean.edge.bc
  return(diff.edge.bc)
}

# Calculate for new edges. Use country level data only, if not you wait for ages. 
# Consider parLapply to speed up things for final analysis. 
edge.num <- unlist(lapply(V(analysis.nw), function(x){length(incident(analysis.nw, v = x, mode = c("all")))}))

mean.edge.bc <- mean_edge_bc(analysis.nw)

#Deadends
set.seed(250589)
new.edges.expl <- which(E(analysis.nw) %in% E(analysis.nw)[from(sample(which(edge.num == 1), 50))])
diff.edge.bc.expl <- unlist(lapply(new.edges.expl, function(x){diff_edge_bc(edge = x, graph = analysis.nw, mean.edge.bc)}))

#Densifying roads
set.seed(250589)
new.edges.dens <- which(E(analysis.nw) %in% E(analysis.nw)[from(sample(which(edge.num != 1), 50))] & 
                     !E(analysis.nw) %in% E(analysis.nw)[from(which(edge.num == 1))])
diff.edge.bc.dens <- unlist(lapply(new.edges.dens, function(x){diff_edge_bc(x, analysis.nw, mean.edge.bc)}))



plot(density(diff.edge.bc.dens), col = "blue")
lines(density(diff.edge.bc.expl), col = "red")

# Plot in space
col.pal <- colorRampPalette(c("yellow","orange","red"),space = "Lab")(256)
col.vec <- col.pal[as.numeric(cut(diff.edge.bc.dens,breaks = 256))]
plot(analysis.spdf)
lines(analysis.spdf[new.edges.dens], col = col.vec, lwd = 5)

# Number of dead-ends and T-Junctions
# Forget about T-Junctions, every crossroad is a T Junctions with us.... 

edge.num <- unlist(lapply(V(analysis.nw), function(x){length(incident(analysis.nw, v = x, mode = c("all")))}))

plot(density(edge.num))
deadend.ratio <- sum(edge.num == 1) / (sum(edge.num != 1) - sum(edge.num == 2))



