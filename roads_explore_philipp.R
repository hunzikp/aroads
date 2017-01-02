
######################################################################################
# EXPLORATORY ROAD NETWORK ANALYSIS
#
# 2016-12-32
# Carl Müller-Crepon w/ Philipp Hunziker, Matthew Simonson, ETHZ, Northestern
######################################################################################

######################################################################################
# INIT
######################################################################################

setwd("/home/hunzikp/Projects/roadvec/aroads")

library(gtools)
library(rgeos)
library(geosphere)
library(rgdal)
library(cshapes)
library(maptools)
library(igraph)
library(viridisLite)
library(ggplot2)
source("sldf2graph.R")


######################################################################################
# FUNCTIONS
######################################################################################

# Generate simplified network with nodes at intersections & endpoints only
simplify_network <- function(graph){
  
  # Find all vertices with only two edges
  edge.num <- unlist(lapply(V(graph), function(x){length(incident(graph, v = x, mode = c("all")))}))
  vertices.bridge <- which(edge.num == 2)
  
  # Sequentially (!) delete old edges and replace them with a bridging one
  for(v in vertices.bridge){
    these.neighbors <- neighbors(graph, v, mode =  "all")
    old.edges <- E(graph)[from(v)]
    # Add new edge
    graph <- add_edges(graph, these.neighbors,
                       attr = list(weight = sum(edge_attr(graph, "weight", old.edges))))
    # Delede old edges
    graph <- delete.edges(graph, E(graph)[from(v)])
  }
  
  # Delete now unconnected vertices
  graph <- delete.vertices(graph, c(vertices.bridge, which(edge.num == 0)))
  
  return(graph)
}


######################################################################################
# LOAD RN DATA
######################################################################################

## Load data
load("/home/hunzikp/Projects/roadvec/output/networks/country/roads_uni_full.RData")
unified.ls <- full.list

######################################################################################
# AGE VS EBC ANALYSIS
######################################################################################

for (cc in 1:length(unified.ls)) {
  
  ## Get this country's sldf and polygon
  unified.sldf <- unified.ls[[cc]]
  this.gwid <- as.numeric(sub("c", "", names(unified.ls)[cc]))
  cshp.all <- cshp(as.Date("2003-12-31"))
  this.cshp <- cshp.all[cshp.all$GWCODE == this.gwid,]
  
  ## Split unified into yearly SLDFs and calculate road age
  unified.col.names <- names(unified.sldf)
  unified.yearcol.names <- unified.col.names[grepl("map_", unified.col.names)]
  unified.yearcol.years <- as.numeric(sub("[[:alnum:]]+_[[:graph:]]+_", "", unified.yearcol.names))
  unified.yearcol.names <- unified.yearcol.names[order(unified.yearcol.years)]
  unified.yearcol.years <- unified.yearcol.years[order(unified.yearcol.years)]
  
  flag.mat <- as.matrix(unified.sldf@data[,unified.yearcol.names])
  firsttrue.vec <- apply(flag.mat, 1, function(x) which(x)[1])
  firstbuilt.vec <- sapply(firsttrue.vec, function(x) unified.yearcol.years[x])
  unified.sldf$firstyear <- firstbuilt.vec
  
  road.panel.ls <- vector("list", length(unified.yearcol.names))
  for (i in 1:length(unified.yearcol.names)) {
    this.year <- unified.yearcol.years[i]
    this.sldf <- unified.sldf[unified.sldf@data[,unified.yearcol.names[i]]==1,]
    this.sldf@data <- this.sldf@data[,!(names(this.sldf) %in% unified.yearcol.names)]
    this.sldf$age <- this.year - this.sldf$firstyear
    road.panel.ls[[i]] <- this.sldf
  }
  
  ## Calculate EBC for most recent network
  newest.sldf <- road.panel.ls[[length(road.panel.ls)]]
  newest.graph <- sldf2graph(newest.sldf)
  edge.bc <- edge_betweenness(newest.graph, e = E(newest.graph), directed = F, weights=E(newest.graph)$length)
  edge.bc.norm <- (edge.bc / ((length(V(newest.graph))-1)*(length(V(newest.graph))-2)))*100
  
  ## Plot normalized edge betweenness
  newest.sldf <- spChFIDs(newest.sldf, rownames(newest.sldf@data))
  shp.df    <- data.frame(id=rownames(newest.sldf@data),
                          values=edge.bc.norm,
                          newest.sldf@data, stringsAsFactors=F)
  data_fort   <- fortify(newest.sldf)
  data_merged <- join(data_fort, shp.df, by="id")
  
  ggplot() +  
    geom_polygon(data=this.cshp, aes(x=long, y=lat, group=group), fill="grey40", colour="grey90", alpha=1) +
    geom_path(data=data_merged, size=0.5, aes(x=long, y=lat, group=group, color=values)) +
    scale_color_gradientn(name="", colours = inferno(32)) +
    coord_equal(ratio=1) + 
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          legend.position="right")
  
  
  plot(newest.sldf$age, edge.bc.norm)
  
  
  
  
  col.pal <- inferno(n=256)
  col.vec <- col.pal[as.numeric(cut(log(1+edge.bc/1000), breaks=256))]
  plot(newest.sldf, col = col.vec)
  
  
  
  plot(density(log(1+edge.bc)))
  
  
}











######################################################################################
# ANALYSIS
######################################################################################

## Set an analysis network
analysis.nw <- roads.nw 
analysis.spdf <- roads.sldf

# Correlation age - edge betweenness centrality #######################

## 1. Calculate edge betweenness centrality
edge.bc <- edge_betweenness(analysis.nw, e = E(analysis.nw), directed = F) # , weights = rep(1, length(E(analysis.nw)))
plot(density(log(1+edge.bc)))

## 2. We can now correlate this with the age of an edge...

# Simple correlation
cor(log(edge.bc+1), analysis.spdf$age)

# Mean age for each cumulative percentage of edge betweenness
efun <- ecdf(edge.bc)
e.vec <- efun(edge.bc)
plot(e.vec, analysis.spdf$age)
lines(lowess(e.vec, analysis.spdf$age), col="red")

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



