
######################################################################################
# EXPLORATORY ROAD NETWORK ANALYSIS
#
# 2016-12-32
# Carl Müller-Crepon w/ Philipp Hunziker, Matthew Simonson, ETHZ, Northestern
######################################################################################

######################################################################################
# INIT
######################################################################################
.libPaths(c(.libPaths(), "/icr/home/carlvs/R/x86_64-redhat-linux-gnu-library/3.2" ))
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

fig.path <- "figures/"

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

## Read Spatial DataFrame from PostgreSQL
dbReadSpatial <- function(con, schemaname="public", tablename, geomcol="the_geom", idcol=NULL, where="") {
  
  # con:          A PostgreSQL connection (from RPostgreSQL)
  # spatial.df:   A Spatial Data Frame object
  # schemaname:   Target schema name
  # tablename:    Target table name
  # geomcol:      Name of the geometry column in the target table (target table may not have more than one geometry column!)
  # idcol:        Name of the column with unique IDs to be used in the ID slot of the spatial objects (not relevant for point data)
  # Added possibility for subsetting table via where clause (Carl)
  
  ## Build query and fetch the target table
  # Get column names
  q.res <- dbSendQuery(con, statement=paste("SELECT column_name FROM information_schema.columns WHERE table_name ='", tablename, "' AND table_schema ='", schemaname, "';", sep=""))
  schema.table = paste(schemaname, ".", tablename, sep="")
  q.df <- fetch(q.res, -1)
  # Some safe programming
  if (!(geomcol %in% q.df[,1])) {stop(paste("No", geomcol, "column in specified table."))}
  if (!is.null(idcol)) {
    if (!(idcol %in% q.df[,1])) {stop(paste("Specified idname '", idcol, "' not found.", sep=""))}
  }
  # Get table
  query <- paste("SELECT", paste(q.df[,1][q.df[,1] != geomcol], collapse=", "), paste(", ST_ASTEXT(", geomcol, ") AS the_geom FROM", sep=""), schema.table, "  ", where, ";")
  t.res <- dbSendQuery(con, statement=query)
  t.df <- fetch(t.res, -1)
  
  ## Get geometry ID column number
  if (!is.null(idcol)) {
    idcolnum <- which(names(t.df) == idcol)
  } else {
    t.df$id.new <- 1:nrow(t.df)
    idcolnum <- which(names(t.df) == "id.new")
  }
  
  ## Get geometry column number
  geomcolnum <- which(names(t.df) == "the_geom")
  
  ## Build spatial data frame using OGR
  write.df <- t.df[,geomcolnum,drop=FALSE]
  names(write.df) <- "WKT"
  filename <- paste("vector_", as.character(format(Sys.time(), "%H_%M_%S")), sep="")
  filename.csv <- paste(filename, ".csv", sep="")
  write.csv(write.df, paste(gsub("[\\]", "/", tempdir()), "/", filename.csv, sep=""), row.names=TRUE)
  down.spdf <- readOGR(dsn=paste(gsub("[\\]", "/", tempdir()), "/", filename.csv, sep=""), layer=filename, verbose=FALSE)
  rv <- file.remove(paste(gsub("[\\]", "/", tempdir()), "/", filename.csv, sep=""))
  data.df <- data.frame(t.df[,-geomcolnum])
  names(data.df) <- names(t.df)[-geomcolnum]  
  
  # For Spatial Points Data Frame  
  if (grepl("POINT", t.df[1,geomcolnum])) {
    spatial.df <-  SpatialPointsDataFrame(down.spdf@coords, data.df, match.ID=FALSE)
  }
  # For Spatial Polygons/Lines Data Frame    
  if (grepl("POLYGON", t.df[1,geomcolnum]) | grepl("LINE", t.df[1,geomcolnum])) {
    spatial.df <- down.spdf
    spatial.df@data <- data.df
    spatial.df <- spChFIDs(spatial.df, paste(t.df[,idcolnum]))
  }
  return(spatial.df)
}

realvertex_vertex_betweenness <- function(graph, normalize=FALSE) {
  
  ## Simplify and get reference table
  simple.graph.ls <- simplify_network(graph)
  simple.graph <- simple.graph.ls[[1]]
  reference.df <- simple.graph.ls[[2]]
  
  ## Calculate edge betweenness on simplified graph
  vbc.simple <- betweenness(simple.graph)
  if (normalize) {
    vbc.simple <- vbc.simple / ((length(V(simple.graph))-1)*(length(V(simple.graph))-2))
  }
  
  ## Merge simple ebc back to original edges
  vertex_attr(simple.graph, "x")
  simple.vbc.df <- data.frame(x = vertex_attr(simple.graph, "x"),
                              y = vertex_attr(simple.graph, "y"),
                              vbc=vbc.simple)
  coordinates(simple.vbc.df) <- c("x","y")
  return(simple.vbc.df)
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
library(RPostgreSQL)
# con_growup <- dbConnect(dbDriver("PostgreSQL"), dbname="dbname", 
#                         port=5432, host="host", 
#                         user="user", password="pw")
## Load data
load("~/Data/Roads/Bluebooks/preliminary/roads_uni_full.RData")
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

  ## Calculate VBC for most recent network
  vbc <- realvertex_vertex_betweenness(newest.graph, normalize=TRUE) 
  vbc$vbc <- vbc$vbc*100
  
  ## Plot normalized edge betweenness
  
  #Fortify edges
  newest.sldf <- spChFIDs(newest.sldf, rownames(newest.sldf@data))
  shp.df    <- data.frame(id=rownames(newest.sldf@data),
                          values=ebc,
                          newest.sldf@data, stringsAsFactors=F)
  data_fort   <- fortify(newest.sldf)
  data_merged <- join(data_fort, shp.df, by="id")
  
  #Fortify vertices
  vbc.df    <- data.frame(id=rownames(vbc@data),
                          long = coordinates(vbc)[,1],
                          lat = coordinates(vbc)[,2],
                          values=vbc$vbc, stringsAsFactors=F)

  
  ebc.plot.ls[[cc]] <- ggplot() +  
    geom_polygon(data=this.cshp, aes(x=long, y=lat, group=group), fill="grey80", colour="grey90", alpha=1) +
    geom_path(data=data_merged, size=0.5, aes(x=long, y=lat, group=group, color=values)) +
    geom_point(data=vbc.df, size=2, aes(x=long, y=lat, color=values)) +
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
labels = c(`500` = "Uganda",`452` = "Ghana",`475` = "Nigeria",`501` = "Kenya",`553` = "Malawi",`551` = "Zambia",`451` = "Sierra Leone",`510` = "Tanzania")

age.plot.df <- do.call("rbind", age.plot.ls)
age.plot.df <- age.plot.df[order(age.plot.df$gwid),]
cors <- ddply(age.plot.df, c("gwid", "name"), summarise, cor = round(cor(year, ebc), 2))
age.plot.df$name <- as.character(age.plot.df$name)
pdf(file=paste0(fig.path,"age_ebc.pdf"), width = 12, height=5)
ggplot(data = age.plot.df, aes(year, ebc, group = name)) + 
  geom_point() + 
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~name, nrow=2) + 
  labs(x="Road Age", y="Edge Betweenness") + 
  geom_text(data=cors, aes(label=paste("r=", cor, sep="")), x=1990, y=20)
dev.off()

## EBC maps
ebc.plot.ord <- lapply(labels[order(labels)], function(x){ebc.plot.ls[[which(labels == x)]]})
pdf(file=paste0(fig.path,"ebc.pdf"), width = 8, height=12)
grid.arrange(grobs=ebc.plot.ord, nrow=4)
dev.off()


###########################################
# New roads and their edge betweenness and deadends
###########################################

new.ebc.all <- NULL
new.deadend.all <- NULL
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
    this.sldf <- unified.sldf[unified.sldf$firstyear <= this.year,]
    this.sldf@data <- this.sldf@data[,!(names(this.sldf) %in% unified.yearcol.names)]
    this.sldf$age <- this.year - this.sldf$firstyear
    road.panel.ls[[i]] <- list(this.sldf, sldf2graph(this.sldf))
  }
  
  # Edge betweenness of new roads
  newebc.panel.mat <- NULL
  for (i in 1:length(unified.yearcol.names)) {
    # EBC
    this.sldf <- road.panel.ls[[i]] 
    this.ebc <-  edge_betweenness(this.sldf[[2]], e = E(this.sldf[[2]]), 
                                  directed = F, weights=E(this.sldf[[2]])$length)
    newebc.panel.mat <- rbind(newebc.panel.mat,
                              cbind(year =  max(this.sldf[[1]]@data$firstyear),
                                    ebc = this.ebc[this.sldf[[1]]@data$firstyear == max(this.sldf[[1]]@data$firstyear)]/
                                      ((length(V(this.sldf[[2]]))-1)*(length(V(this.sldf[[2]]))-2))*100))
    
    #Deadends
    new.vert <- get.edges(this.sldf[[2]], E(this.sldf[[2]])[this.sldf[[1]]@data$firstyear == max(this.sldf[[1]]@data$firstyear)])
    edge.num <- unlist(lapply(new.vert, function(x){length(incident(this.sldf[[2]], v = x, mode = c("all")))}))
    new.deadend.all <- rbind(new.deadend.all,
                             cbind(year =  max(this.sldf[[1]]@data$firstyear),
                                   deadend.ratio = sum(edge.num == 1) / (sum(edge.num != 2)),
                                   gwid = this.gwid))
  }
  newebc.panel.mat <- data.frame(newebc.panel.mat)
  color.scale <- unique(inferno(256)[cut(newebc.panel.mat$year, 256)])
  names(color.scale) <- unique(newebc.panel.mat$year)
  newebc.panel.mat$year <- as.character(newebc.panel.mat$year)
  
  # Density plot per year -- one plot
  ggplot(newebc.panel.mat[newebc.panel.mat$year != "2003",], aes(ebc, color = year, fill = year)) +
    geom_density(alpha = 0.3) + 
    scale_color_manual(values=color.scale) + 
    scale_fill_manual(values=color.scale) + 
    xlab("Edge betweenness centrality") + ylab("Density")
  
  #  Density plot per year -- multi-plot
  ggplot(newebc.panel.mat, aes(ebc, color = year, fill = year)) +
    geom_density(alpha = 0.3) + 
    scale_color_manual(values=color.scale) + 
    scale_fill_manual(values=color.scale) + facet_wrap(~year, nrow = 1) + 
    xlab("Edge betweenness centrality") + ylab("Density")
  
  # Mean over years
  newebc.panel.mat$year <- as.numeric(newebc.panel.mat$year)
  newebc.panel.agg <- aggregate.data.frame(newebc.panel.mat$ebc, list(newebc.panel.mat$year), FUN = mean)
  colnames(newebc.panel.agg) <- c("year","ebc")
  ggplot(newebc.panel.agg, aes(x = year, y = ebc)) +
    geom_point() +  xlab("Year") + ylab("Edge betweenness centrality")
  
  #Save mean ebc of new roads for later cross-colony plotting
  new.ebc.all <- rbind(new.ebc.all,
                       cbind(gwid = this.gwid, newebc.panel.agg))
  
  
  
  #Exploratory vs densifying roads 
  
  ## Plot normalized edge betweenness of new edges
  #   plot.ls <- list()
  #   for( i in road.panel.ls){
  #     newest.sldf <- spChFIDs(i[[1]], rownames(i[[1]]@data))
  #     newest.graph <- sldf2graph(newest.sldf)
  #     edge.bc <- edge_betweenness(newest.graph, e = E(newest.graph), directed = F, weights=E(newest.graph)$length)
  #     edge.bc.norm <- (edge.bc / ((length(V(newest.graph))-1)*(length(V(newest.graph))-2)))*100
  #     shp.df    <- data.frame(id=rownames(newest.sldf@data),
  #                             values=edge.bc.norm,
  #                             newest.sldf@data, stringsAsFactors=F)
  #     shp.df$values[newest.sldf@data$firstyear != max(newest.sldf@data$firstyear, na.rm = T)] <- NA
  #     data_fort   <- fortify(newest.sldf)
  #     data_merged <- join(data_fort, shp.df, by="id")
  #     plot.ls <- c(plot.ls, list(data_merged))
  #   }
  # 
  #   
  #   ggplot() +  
  #     geom_polygon(data=this.cshp, aes(x=long, y=lat, group=group), fill="grey40", colour="grey90", alpha=1) +
  #     geom_path(data=plot.ls[[7]], size=0.5, aes(x=long, y=lat, group=group, color=values)) +
  #     scale_color_gradientn(limits = c(0,25), name="", colours = inferno(32)) +
  #     coord_equal(ratio=1) + 
  #     theme(axis.title=element_blank(),
  #           axis.text=element_blank(),
  #           axis.ticks=element_blank(),
  #           legend.position="right")
  
  
}


# Cross colony plot ebc
labels = c(`500` = "Uganda",`452` = "Ghana",`475` = "Nigeria",`501` = "Kenya",`553` = "Malawi",`551` = "Zambia",`451` = "Sierra Leone",`510` = "Tanzania")
new.ebc.all$cowname <- ""
for(l in 1:length(labels)){
  new.ebc.all$cowname[new.ebc.all$gwid == as.numeric(names(labels)[l])] <- labels[l]
}
pdf(file=paste0(fig.path,"newroads_ebc.pdf"), width = 12, height=5)
ggplot(new.ebc.all, aes(x = year, y = ebc, group = cowname)) +
  geom_point()  + facet_wrap(~cowname, nrow = 2) +  xlab("Year") + ylab("Edge betweenness centrality") + 
  geom_smooth(method = "lm",formula = y ~x + I(x^2), se = F) ## formula = y~x+I(x^2)
dev.off()
### Interpretable as decreasing marginal utility of new roads.


# Cross colony plot dead ends
labels = c(`500` = "Uganda",`452` = "Ghana",`475` = "Nigeria",`501` = "Kenya",`553` = "Malawi",`551` = "Zambia",`451` = "Sierra Leone",`510` = "Tanzania")
new.deadend.all <- data.frame(new.deadend.all)
new.deadend.all$cowname <- ""
for(l in 1:length(labels)){
  new.deadend.all$cowname[new.deadend.all$gwid == as.numeric(names(labels)[l])] <- labels[l]
}
pdf(file=paste0(fig.path,"newroads_deadends.pdf"), width = 12, height=5)
ggplot(new.deadend.all, aes(x = year, y = deadend.ratio, group = cowname)) +
  geom_point()  + facet_wrap(~cowname, nrow = 2) +  xlab("Year") + ylab("Fraction of deadends") + 
  geom_smooth(method = "lm",formula = y ~x , se = F) ## formula = y~x+I(x^2)
dev.off()

#################################
## Maps over time ##########
#################################
lakes.spdf <- dbReadSpatial(con_growup, schemaname="carlvs", tablename = "lakes", geomcol="the_geom", idcol=NULL, where = "")

length.year.mat <- NULL
for (cc in 1:length(unified.ls)) {
  ## Get this country's sldf and polygon
  unified.sldf <- unified.ls[[cc]]
  this.gwid <- as.numeric(sub("c", "", names(unified.ls)[cc]))
  cshp.all <- cshp(as.Date("2003-12-31"))
  this.cshp <- cshp.all[cshp.all$GWCODE == this.gwid,]
  unified.sldf <- raster::crop(unified.sldf, this.cshp)
  this.lakes <- raster::crop(lakes.spdf, this.cshp)
  ## Road Length
  length.km <- lapply(unified.sldf@lines, 
                      function(x){sum(unlist(lapply(c(2:nrow(x@Lines[[1]]@coords)), 
                                                    function(y){distHaversine(x@Lines[[1]]@coords[y-1,],
                                                                              x@Lines[[1]]@coords[y,])})))})
  unified.sldf$length.km <- unlist(length.km)/1000
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
    this.sldf <- unified.sldf[unified.sldf$firstyear <= this.year,]
    this.sldf@data <- this.sldf@data[,!(names(this.sldf) %in% unified.yearcol.names)]
    this.sldf$age <- this.year - this.sldf$firstyear
    road.panel.ls[[i]] <- list(this.sldf)
  }
  
  # Plot pure
  #   for(i in road.panel.ls){
  #     sldf <- fortify(i[[1]])
  #     p <- ggplot() +  
  #       geom_polygon(data=this.cshp, aes(x=long, y=lat, group=group), fill="grey80", colour="grey90", alpha=1) +
  #       geom_path(data=sldf, size=0.5, aes(x=long, y=lat, group = group)) + coord_equal(ratio=1) + 
  #       theme(axis.title=element_blank(),
  #             axis.text=element_blank(),
  #             axis.ticks=element_blank(),
  #             legend.position="right")
  #     print(p)
  #   }
  # Plot with new roads in red
  col.new.road = "red"
  p.list <- list()
  count = 0
  for(i in road.panel.ls){
    count = count + 1
    #Road plot
    shp.df    <- data.frame(id=i[[1]]@data$hlid,i[[1]]@data, stringsAsFactors=F)
    sldf <- join(fortify(i[[1]]), shp.df, by = "id")
    any.new <- any(sldf$firstyear == unified.yearcol.years[count])
    if(any.new){
      p <- ggplot() +  
        geom_polygon(data=this.cshp, aes(x=long, y=lat, group=group), fill="grey80", colour="grey90", alpha=1) +
        geom_path(data=sldf, size=0.5, aes(x=long, y=lat, group = group)) + coord_equal(ratio=1) + 
        geom_path(data=sldf[sldf$firstyear == unified.yearcol.years[count],], size=0.5, aes(x=long, y=lat, group = group, color = col.new.road)) + 
        ggtitle(unified.yearcol.years[count]) + 
        theme(axis.title=element_blank(),
              axis.text=element_blank(),
              axis.ticks=element_blank(),
              legend.position="none")
    } else {
      p <- ggplot() +  
        geom_polygon(data=this.cshp, aes(x=long, y=lat, group=group), fill="grey80", colour="grey90", alpha=1) +
        geom_path(data=sldf, size=0.5, aes(x=long, y=lat, group = group)) + coord_equal(ratio=1) + 
        ggtitle(unified.yearcol.years[count]) + 
        theme(axis.title=element_blank(),
              axis.text=element_blank(),
              axis.ticks=element_blank(),
              legend.position="none")
    }

    p.list <- c(p.list,list(p))
    
    ## Calculate Road density
    
    area <- try((areaPolygon(this.cshp) - sum(areaPolygon(this.lakes)))/1e6)
    if(class(area) == "try-error"){area <- (areaPolygon(this.cshp))/1e6}
    road.length <- sum(i[[1]]@data$length.km)
    length.year.mat <- rbind(length.year.mat,
                             c(area = area, road.length = road.length, year = max(i[[1]]@data$firstyear),gwid = this.gwid))
    
  }
  png(file=paste0(fig.path,this.gwid,"_mappanel.png"), width = 10, height=5, units = "in",res = 600)
  grid.arrange(grobs = p.list, nrow = 2)
  dev.off()
}

#Plot road length/density over time
length.year.mat <- data.frame(length.year.mat)
labels = c(`500` = "Uganda",`452` = "Ghana",`475` = "Nigeria",`501` = "Kenya",`553` = "Malawi",`551` = "Zambia",`451` = "Sierra Leone",`510` = "Tanzania")
length.year.mat$cowname <- ""
for(l in 1:length(labels)){
  length.year.mat$cowname[length.year.mat$gwid == as.numeric(names(labels)[l])] <- labels[l]
}
pdf(file=paste0(fig.path,"road_density.pdf"), width = 12, height=5)
ggplot(length.year.mat, aes(x = year, y = road.length/area, group = cowname)) +
  geom_point()  + facet_wrap(~cowname, nrow = 2) +  xlab("Year") + ylab("Road density (km/km^2)") + 
  geom_smooth(method = "loess", se = F) ## formula = y~x+I(x^2)
dev.off()
##################################
# Plot some (apparent) anomalies
##################################

# Rivers
rivers.spdf <- dbReadSpatial(con_growup, schemaname="carlvs", tablename = "rivers", geomcol="the_geom", idcol=NULL, where = "")
lakes.spdf <- dbReadSpatial(con_growup, schemaname="carlvs", tablename = "lakes", geomcol="the_geom", idcol=NULL, where = "")

# Uganda & Nigeria

# Uganda
cc <- "c500"
year <- "map_bb_1945"
## Get this country's sldf and polygon
unified.sldf <- unified.ls[[cc]]
unified.sldf <- unified.sldf[unified.sldf@data[,year],]
this.gwid <- as.numeric(sub("c", "", cc))
cshp.all <- cshp(as.Date("2003-12-31"))
this.cshp <- cshp.all[cshp.all$GWCODE == this.gwid,]

# Crop roads
unified.sldf <- raster::crop(unified.sldf, this.cshp)

## Crop rivers
this.rivers <- raster::crop(rivers.spdf, this.cshp)
shp.df    <- data.frame(id = rownames(this.rivers@data), this.rivers@data, stringsAsFactors = F)
this.rivers <- join(fortify(this.rivers), shp.df, by = "id")

## Crop Lakes
this.lakes <- raster::crop(lakes.spdf, this.cshp)

## Plot
p <- ggplot() +  
  geom_polygon(data=this.cshp, aes(x=long, y=lat, group=group), fill="grey80", colour="grey90", alpha=1) +
  geom_polygon(data=this.lakes, aes(x=long, y=lat, group=group, fill="Water", colour="Water"), alpha=1) +
  geom_path(data=this.rivers,  aes(x=long, y=lat, group = group, 
                                   size = strokeweig, colour="Water"))  + scale_size(range = c(0, 2)) +
  geom_path(data=unified.sldf, size=0.5, aes(x=long, y=lat, group = group, colour = "Road")) + coord_equal(ratio=1) + 
  scale_colour_manual("", 
                      breaks = c("Road", "Water"),
                      values = c("black", "dodgerblue3")) + 
  scale_fill_manual("", 
                    breaks = c("Water"),
                    values = c("dodgerblue3")) +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.justification=c(1,0), legend.position=c(1,0),
        legend.background = element_rect(fill="transparent"),
        legend.text=element_text(size=12)) + 
  guides(colour=guide_legend(override.aes=list(fill="grey90")), fill = F, size = F)
png(file=paste0(fig.path,"500_water.png"), width = 7, height=7, units = "in",res = 600)
print(p)
dev.off()

# Nigeria
cc <- "c475"
year <- "map_182_1941"
## Get this country's sldf and polygon
unified.sldf <- unified.ls[[cc]]
unified.sldf <- unified.sldf[unified.sldf@data[,year],]
this.gwid <- as.numeric(sub("c", "", cc))
cshp.all <- cshp(as.Date("2003-12-31"))
this.cshp <- cshp.all[cshp.all$GWCODE == this.gwid,]

# Crop roads
unified.sldf <- raster::crop(unified.sldf, this.cshp)

## Crop rivers
this.rivers <- raster::crop(rivers.spdf, this.cshp)
shp.df    <- data.frame(id = rownames(this.rivers@data), this.rivers@data, stringsAsFactors = F)
this.rivers <- join(fortify(this.rivers), shp.df, by = "id")

## Crop Lakes
this.lakes <- raster::crop(lakes.spdf, this.cshp)

## Plot
p <- ggplot() +  
  geom_polygon(data=this.cshp, aes(x=long, y=lat, group=group), fill="grey80", colour="grey90", alpha=1) +
  geom_polygon(data=this.lakes, aes(x=long, y=lat, group=group, fill="Water", colour="Water"), alpha=1) +
  geom_path(data=this.rivers,  aes(x=long, y=lat, group = group, 
                                   size = strokeweig, colour="Water"))  + scale_size(range = c(0, 2)) +
  geom_path(data=unified.sldf, size=0.5, aes(x=long, y=lat, group = group, colour = "Road")) + coord_equal(ratio=1) + 
  scale_colour_manual("", 
                      breaks = c("Road", "Water"),
                      values = c("black", "dodgerblue3")) + 
  scale_fill_manual("", 
                    breaks = c("Water"),
                    values = c("dodgerblue3")) +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.justification=c(1,0), legend.position=c(1,0),
        legend.background = element_rect(fill="transparent"),
        legend.text=element_text(size=12)) + 
  guides(colour=guide_legend(override.aes=list(fill="grey90")), fill = F, size = F)

png(file=paste0(fig.path,"475_water.png"), width = 9, height=7, units = "in",res = 600)
print(p)
dev.off()

# Something is wrong here - Rivers and roads do not allign... Check projections




# Sierra Leone and Railway
cc <- "c451"
year <- "map_bb_1917"
## Get this country's sldf and polygon
unified.sldf <- unified.ls[[cc]]
unified.sldf <- unified.sldf[unified.sldf@data[,year],]
this.gwid <- as.numeric(sub("c", "", cc))
cshp.all <- cshp(as.Date("2003-12-31"))
this.cshp <- cshp.all[cshp.all$GWCODE == this.gwid,]

# Railroad
rail_ne.spdf <- readShapeLines("~/Data/geodata/Railroads/Replication_Files_Jedwab_Moradi/GIS_Files_for_Africa/Railroads/Built/Railroads.shp")
rail_ne.spdf <- rail_ne.spdf[rail_ne.spdf$year_built < 1918, ]
rail_ne.spdf <- raster::crop(rail_ne.spdf, this.cshp)

## Plot
p <- ggplot() +  
  geom_polygon(data=this.cshp, aes(x=long, y=lat, group=group), fill="grey80", colour="grey90", alpha=1) +
  geom_path(data=rail_ne.spdf, aes(x=long, y=lat, group=group, color = "Railroad"), alpha=1) +
  geom_path(data=unified.sldf, size=0.5, aes(x=long, y=lat, group = group, color = "Road")) + coord_equal(ratio=1) + 
  scale_colour_manual("", 
                      breaks = c("Railroad", "Road"),
                      values = c("firebrick3", "black")) +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.justification=c(1,0), legend.position=c(1,0),
        legend.background = element_rect(fill="transparent"),
        legend.text=element_text(size=12))
png(file=paste0(fig.path,"451_railway.png"), width = 7, height=7, units = "in",res = 600)
print(p)
dev.off()


#################################
# ROADS TIMELINE / TABLE WITH YEARS per COUNTRY
######################################
labels = c(`500` = "Uganda",`452` = "Ghana",`475` = "Nigeria",`501` = "Kenya",
           `553` = "Malawi",`551` = "Zambia",`451` = "Sierra Leone",`510` = "Tanzania")
labels <- labels[order(labels, decreasing = F)]
# Make plot of country years for which road lists are converted
year.start <- 1890
year.end <- 2010
pdf(file=paste0(fig.path,"roads_timeline.pdf"), width = 12, height=5)
plot(NA, xlim=c(year.start,year.end),
     bty='n',yaxt="n", xaxt="n", ylab="", xlab="Year", col="black", ylim = c(0,length(labels)))

for(c in seq(1, length(labels), by = 2)){
  rect(year.start, c-0.5, year.end, c+0.5, 
       col = "lightgrey", border = NA)
}

count = length(labels) + 1
for(cc in 1:length(labels)){
  count = count - 1
  unified.sldf <- unified.ls[[paste0("c",names(labels)[cc])]]
  unified.col.names <- names(unified.sldf)
  unified.yearcol.names <- unified.col.names[grepl("map_", unified.col.names)]
  unified.yearcol.years <- as.numeric(sub("[[:alnum:]]+_[[:graph:]]+_", "", unified.yearcol.names))
  # Michelin
  points(unified.yearcol.years[!grepl("_bb_",unified.yearcol.names)],rep(count, length(unified.yearcol.years[!grepl("_bb_",unified.yearcol.names)])), 
         pch=3, col="black", cex = 1.75)
  # Bluebook
  points(unified.yearcol.years[grepl("_bb_",unified.yearcol.names)],rep(count, length(unified.yearcol.years[grepl("_bb_",unified.yearcol.names)])), 
         pch=4, col="black", cex = 1.75)
}
axis(1, at=c(year.start,seq(1900, 2000, by= 25), year.end))
labels <- labels[order(labels, decreasing = T)]
axis(2, at=c(1:length(labels)), labels=labels, las=1, tick=F, lwd=0, pos=1890)
dev.off()

