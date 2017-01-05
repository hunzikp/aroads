
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
source("cluster_slnodes.R")


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

voronoipolygons = function(layer) {
  require(deldir)
  crds = layer@coords
  z = deldir(crds[,1], crds[,2])
  w = tile.list(z)
  polys = vector(mode='list', length=length(w))
  require(sp)
  for (i in seq(along=polys)) {
    pcrds = cbind(w[[i]]$x, w[[i]]$y)
    pcrds = rbind(pcrds, pcrds[1,])
    polys[[i]] = Polygons(list(Polygon(pcrds)), ID=as.character(i))
  }
  SP = SpatialPolygons(polys)
  voronoi = SpatialPolygonsDataFrame(SP, data=data.frame(x=crds[,1], 
                                                         y=crds[,2], row.names=sapply(slot(SP, 'polygons'), 
                                                                                      function(x) slot(x, 'ID'))))
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

pdf("/home/hunzikp/Projects/roadvec/output/plots/age_ebc.pdf", height=12, width=8)
ggplot(data = age.plot.df, aes(year, ebc)) + 
  geom_point() + 
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~name, nrow=4) + 
  labs(x="Road Age", y="Edge Betweenness") + 
  geom_text(data=cors, aes(label=paste("r=", cor, sep="")), x=1990, y=20)
dev.off()

## EBC maps
pdf("/home/hunzikp/Projects/roadvec/output/plots/ebc.pdf", height=12, width=8)
grid.arrange(grobs=ebc.plot.ls, nrow=4)
dev.off()


######################################################################################
# COMMUNITIES
######################################################################################

process.countries <- c(1:8)
community.sp.ls <- vector("list", length(process.countries))
for (j in 1:length(process.countries)) {
  
  cc <- process.countries[j]
  
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
    this.sldf <- unified.sldf[unified.sldf@data[,unified.yearcol.names[i]],]
    this.sldf@data <- this.sldf@data[,!(names(this.sldf) %in% unified.yearcol.names)]
    this.sldf$age <- this.year - this.sldf$firstyear
    road.panel.ls[[i]] <- this.sldf
  }
  
  ## Get most recent network, cluster nodes, and simplify
  newest.sldf <- road.panel.ls[[length(unified.yearcol.years)]]
  newest.sldf <- spChFIDs(newest.sldf, rownames(newest.sldf@data))
  clustered.sldf <- cluster_slnodes(newest.sldf , eps=0.05, MinPts=2)
  clustered.graph <- sldf2graph(clustered.sldf)
  E(clustered.graph)$weight <- E(clustered.graph)$length
  simple.graph.ls <- simplify_network(clustered.graph)
  simple.graph <- simple.graph.ls[[1]]
  
  ## Remove small connected clusters from simple graph
  cls <- clusters(simple.graph)
  if (any(cls$csize==2)) {
    single.edge.cls <- which(cls$csize==2)
    for (idx in single.edge.cls) {
      this.vertex.ids <- which(cls$membership==idx)
      this.verteces <- V(simple.graph)[this.vertex.ids[1]]
      this.edge <- E(simple.graph)[from(this.verteces[1])]
      simple.graph <- delete.edges(simple.graph, this.edge)
    }
    edge.num <- unlist(lapply(V(simple.graph), function(x){length(incident(simple.graph, v = x, mode = c("all")))}))
    simple.graph <- delete.vertices(simple.graph, which(edge.num == 0))
  }
  
  ## Make simple sldf from simple graph
  coords <- cbind(vertex_attr(simple.graph, "x", index = V(simple.graph)), 
                  vertex_attr(simple.graph, "y", index = V(simple.graph)))
  edges <- get.edges(simple.graph, E(simple.graph))
  lines <- lapply(c(1:nrow(edges)), function(x){Lines(list(Line(coords[edges[x,],])), as.character(x))})
  simple.sldf <- SpatialLinesDataFrame(SpatialLines(lines), data.frame(id=E(simple.graph)$id), FALSE)
  
  ## Make community dendrogram
  community <- cluster_edge_betweenness(simple.graph, weights=E(simple.graph)$weight, directed=FALSE)
  
  N <- 6
  split.voronoi.ls <- vector("list", N-1)
  split.network.ls <- split.voronoi.ls
  for (k in 2:N) {

    # Calculate communities for given k
    community$membership <- cut_at(community, no=k)
    
    # Crossing edges
    crossing.bool <- crossing(community, simple.graph)
    crossing.ids <- E(simple.graph)$id[crossing.bool]
    
    # Create community SLDF (with crossing edges removed)
    this.sldf <- simple.sldf
    this.sldf <- this.sldf[!(this.sldf$id %in% crossing.ids),]
    
    # Make areas for communities using voronoi tesselations
    nodes.df <- data.frame(x=V(simple.graph)$x, y=V(simple.graph)$y, community=community$membership)
    nds.sp <- SpatialPoints(nodes.df[,c(1:2)])
    vrn.spdf <- voronoipolygons(nds.sp)
    vrn.Polygons.ls <- vector("list", k)
    for (i in 1:k) {
      this.vrn.spdf <- vrn.spdf[nodes.df$community==i,]
      this.vrn.spdf <- gUnionCascaded(this.vrn.spdf)
      this.vrn.sp <- gIntersection(this.vrn.spdf, this.cshp)
      vrn.Polygons.ls[[i]] <- Polygons(this.vrn.sp@polygons[[1]]@Polygons, as.character(i))
    }
    vrn.sp <- SpatialPolygons(vrn.Polygons.ls)
    
    split.voronoi.ls[[k-1]] <- vrn.sp
    split.network.ls[[k-1]] <- this.sldf
  }
  
  community.sp.ls[[j]] <- list(split.voronoi.ls, split.network.ls)
  names(community.sp.ls)[j] <- paste0("c", this.gwid)
}


#### Plot two communities for all countries
tc.plot.ls <- vector("list", length(community.sp.ls))
for (cc in 1:length(community.sp.ls)) {
  
  k <- 2
  
  vrn.sp.ls <- community.sp.ls[[cc]][[1]]
  community.sldf.ls <- community.sp.ls[[cc]][[2]]
  vrn.sp <- vrn.sp.ls[[1]]
  community.sldf <- community.sldf.ls[[1]]
  
  vrn.df <- data.frame(id=1:length(vrn.sp))
  vrn.spdf <- SpatialPolygonsDataFrame(vrn.sp, vrn.df, FALSE)
  vrn.points = fortify(vrn.spdf, region="id")
  roads.fort <- fortify(community.sldf)
  
  tc.plot.ls[[cc]] <- ggplot() + 
    geom_polygon(data=vrn.points, aes_string("long", "lat", group="group", fill="id")) + 
    geom_path(data=roads.fort, aes(long, lat, group=group)) +
    coord_equal(ratio=1) + 
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          legend.position="none") +
    scale_fill_manual(values=viridis(k)[1:k], drop=F)
}

grid.arrange(grobs=tc.plot.ls)


#### Plot two communities for nigeria
vrn.sp.ls <- community.sp.ls$c475[[1]]
community.sldf.ls <- community.sp.ls$c475[[2]]
vrn.sp <- vrn.sp.ls[[1]]
community.sldf <- community.sldf.ls[[1]]

vrn.df <- data.frame(id=1:length(vrn.sp))
vrn.spdf <- SpatialPolygonsDataFrame(vrn.sp, vrn.df, FALSE)
vrn.points = fortify(vrn.spdf, region="id")
roads.fort <- fortify(community.sldf)

pdf("/home/hunzikp/Projects/roadvec/output/plots/nigeria_2.pdf")
ggplot() + 
geom_polygon(data=vrn.points, aes_string("long", "lat", group="group", fill="id")) + 
geom_path(data=roads.fort, aes(long, lat, group=group)) +
coord_equal(ratio=1) + 
theme(axis.title=element_blank(),
      axis.text=element_blank(),
      axis.ticks=element_blank(),
      legend.position="none") +
scale_fill_manual(values=viridis(k)[1:k], drop=F)
dev.off()

#### Plot Nigeria with six communities, next to six largest ethnic groups
# communties
k <- 6
vrn.sp.ls <- community.sp.ls$c475[[1]]
community.sldf.ls <- community.sp.ls$c475[[2]]
vrn.sp <- vrn.sp.ls[[k-1]]
community.sldf <- community.sldf.ls[[k-1]]

vrn.df <- data.frame(id=1:length(vrn.sp))
vrn.spdf <- SpatialPolygonsDataFrame(vrn.sp, vrn.df, FALSE)
vrn.points = fortify(vrn.spdf, region="id")
roads.fort <- fortify(community.sldf)

community.gplot <- ggplot() + 
  geom_polygon(data=vrn.points, aes_string("long", "lat", group="group", fill="id")) + 
  geom_path(data=roads.fort, aes(long, lat, group=group)) +
  coord_equal(ratio=1) + 
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position="none") +
  scale_fill_manual(values=viridis(k)[1:k], drop=F)

# ethnic groups
this.greg <- greg.spdf[greg.spdf$COW==475,]
grpnums <- unique(this.greg$GROUP1)
grpnames <- as.character(unique(this.greg$G1SHORTNAM))
greg.Polygons.ls <- vector("list", length(grpnums))
for (i in 1:length(grpnums)) {
  grpnum <- grpnums[i]
  this.group <- this.greg[this.greg$GROUP1==grpnum,]
  this.group <- gUnionCascaded(this.group)
  greg.Polygons.ls[[i]] <- Polygons(this.group@polygons[[1]]@Polygons, as.character(i))
}
this.greg.sp <- SpatialPolygons(greg.Polygons.ls)
this.greg.spdf <- SpatialPolygonsDataFrame(this.greg.sp, data.frame(id=1:length(grpnums), name=as.character(grpnames), stringsAsFactors = FALSE), FALSE)
this.greg.spdf$size <- gArea(this.greg.spdf, byid=TRUE)
largest.ids <- (this.greg.spdf$id[order(-this.greg.spdf$size)])[1:k]
largest.greg.spdf <- this.greg.spdf[this.greg.spdf$id %in% largest.ids,]

largest.points <- fortify(largest.greg.spdf, region="id")
largest.merged <- join(largest.points, largest.greg.spdf@data, by="id")
ethnic.gplot <- ggplot() + 
  geom_polygon(data=largest.merged, aes(long, lat, group=group, fill=as.factor(name))) + 
  coord_equal(ratio=1) + 
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position="right") +
  scale_fill_manual("", values=viridis(k)[c(3,1,6,5,2,4)], drop=F)

pdf("/home/hunzikp/Projects/roadvec/output/plots/nigeria_6.pdf", height=4, width=10)
grid.arrange(community.gplot, ethnic.gplot, ncol=2, widths=c(5.6,8))
dev.off()

