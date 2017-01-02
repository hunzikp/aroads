
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
library(gridExtra)
source("sldf2graph.R")


######################################################################################
# FUNCTIONS
######################################################################################

# Generate simplified network with nodes at intersections & endpoints only
simplify_network <- function(graph) {
  
  # Find all vertices with only two edges
  edge.num <- unlist(lapply(V(graph), function(x){length(incident(graph, v = x, mode = c("all")))}))
  vertices.bridge <- which(edge.num == 2)
  
  E(graph)$id <- 1:length(E(graph))
  reference.ls <- vector("list", length(vertices.bridge))
  
  out.graph <- graph
  
  # Sequentially (!) delete old edges and replace them with a bridging one
  new.id <- max(E(graph)$id) + 1
  cntr <- 1
  for(v in vertices.bridge){
    these.neighbors <- neighbors(out.graph, v, mode =  "all")
    old.edges <- E(out.graph)[from(v)]
    old.ids <- old.edges$id
    
    # Add new edge
    out.graph <- add_edges(out.graph, these.neighbors,
                       attr = list(weight = sum(edge_attr(out.graph, "weight", old.edges)), id=new.id))
    
    # Delede old edges
    out.graph <- delete.edges(out.graph, E(out.graph)[from(v)])
    
    # Update reference list
    reference.ls[[cntr]] <- cbind(new.id, old.ids)
    
    # Increase id
    new.id <- new.id + 1
    cntr <- cntr + 1
  }
  
  # Delete now unconnected vertices
  out.graph <- delete.vertices(out.graph, c(vertices.bridge, which(edge.num == 0)))
  
  # Make reference table
  reference.mat <- do.call("rbind", reference.ls)
  reference.df <- data.frame(original.id=E(graph)$id, new.id=NA)
  replaced.edges <- E(graph)$id[E(graph)$id %in% reference.mat[,2]]
  reference.df$new.id[!(reference.df$original.id %in% replaced.edges)] <- reference.df$original.id[!(reference.df$original.id %in% replaced.edges)]
  for (i in 1:length(replaced.edges)) {
    original.id <- replaced.edges[i]
    parent.id <- reference.mat[reference.mat[,2]==original.id,1]
    hasparent <- TRUE
    while (hasparent) {
      candidate.parent.id <- reference.mat[reference.mat[,2]==parent.id,1]
      hasparent <- length(candidate.parent.id) > 0
      if (hasparent) {
        parent.id <- candidate.parent.id
      }
    }
    reference.df$new.id[reference.df$original.id==original.id] <- parent.id
  }
  
  return(list(out.graph, reference.df))
}

realvertex_edge_betweenness <- function(graph, normalize=FALSE) {
  
  ## Simplify and get reference table
  simple.graph.ls <- simplify_network(graph)
  simple.graph <- simple.graph.ls[[1]]
  reference.df <- simple.graph.ls[[2]]
  
  ## Calculate edge betweenness on simplified graph
  ebc.simple <- edge_betweenness(simple.graph)
  if (normalize) {
    ebc.simple <- ebc.simple / ((length(V(simple.graph))-1)*(length(V(simple.graph))-2))
  }
  
  ## Merge simple ebc back to original edges
  simple.ebc.df <- data.frame(new.id=E(simple.graph)$id, ebc=ebc.simple)
  reference.df <- merge(reference.df, simple.ebc.df, by="new.id", all.x=TRUE, all.y=FALSE)
  
  ## Order and out
  reference.df <- reference.df[order(reference.df$original.id),]
  ebc <- reference.df$ebc
  return(ebc)
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

ebc.plot.ls <- vector("list", length(unified.ls))
age.plot.ls <- vector("list", length(unified.ls))
for (cc in 1:length(unified.ls)) {
  
  ## Get this country's sldf and polygon
  unified.sldf <- unified.ls[[cc]]
  this.gwid <- as.numeric(sub("c", "", names(unified.ls)[cc]))
  cshp.all <- cshp(as.Date("2003-12-31"))
  this.cshp <- cshp.all[cshp.all$GWCODE == this.gwid,]
  this.name <- as.character(this.cshp$CNTRY_NAME)
  
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
  E(newest.graph)$weight <- E(newest.graph)$length
  ebc <- realvertex_edge_betweenness(newest.graph, normalize=TRUE) 
  ebc <- ebc*100
  
  ## Plot normalized edge betweenness
  newest.sldf <- spChFIDs(newest.sldf, rownames(newest.sldf@data))
  shp.df    <- data.frame(id=rownames(newest.sldf@data),
                          values=ebc,
                          newest.sldf@data, stringsAsFactors=F)
  data_fort   <- fortify(newest.sldf)
  data_merged <- join(data_fort, shp.df, by="id")
  
  ebc.plot.ls[[cc]] <- ggplot() +  
    geom_polygon(data=this.cshp, aes(x=long, y=lat, group=group), fill="grey40", colour="grey90", alpha=1) +
    geom_path(data=data_merged, size=0.5, aes(x=long, y=lat, group=group, color=values)) +
    scale_color_gradientn(name="", colours = inferno(32)) +
    coord_equal(ratio=1) + 
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          legend.position="right") +
    labs(title=this.name)
  
  ## Prepare data for age vs betweenness plots
  pl.df <- data.frame(year=2003-newest.sldf$age, ebc=ebc, gwid=this.gwid, name=this.name )
  age.plot.ls[[cc]] <- pl.df
}

## Edge vs betweenness plots
age.plot.df <- do.call("rbind", age.plot.ls)
age.plot.df <- age.plot.df[order(age.plot.df$gwid),]
cors <- ddply(age.plot.df, c("gwid", "name"), summarise, cor = round(cor(year, ebc), 2))

ggplot(data = age.plot.df, aes(year, ebc)) + 
  geom_point() + 
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~name) + 
  labs(x="Road Age", y="Edge Betweenness") + 
  geom_text(data=cors, aes(label=paste("r=", cor, sep="")), x=1990, y=20)


## EBC maps
grid.arrange(grobs=ebc.plot.ls, nrow=3)






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



