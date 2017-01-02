######################################################################################
# LINE UNIFICATION HELPER FUNCTIONS
#
# 2016-12-07
######################################################################################

require(sp)
require(data.table)
require(plyr)

###############################################################
# CREATE VERTEX TABLES
###############################################################

make_vertex_tables <- function(lines.sldf, vid.prefix) {
  
  vertex.coord.ls <- lapply(lines.sldf@lines, function(x) x@Lines[[1]]@coords)
  vertex.coord.mat <- unique(do.call('rbind', vertex.coord.ls))
  rownames(vertex.coord.mat) <- NULL
  vertex.df <- data.frame(vid=paste0(vid.prefix, 1:nrow(vertex.coord.mat)))
  vertex.spdf <- SpatialPointsDataFrame(vertex.coord.mat, vertex.df)
  
  return(vertex.spdf)
}


###############################################################
# CREATE LINE DATA FRAMES
###############################################################

make_line_df <- function(lines.sldf, vertex.spdf, lidname="lid") {
  
  vertex.df <- vertex.spdf@data
  vertex.coord.mat <- vertex.spdf@coords
  
  ## Make line.df
  line.coord.ls <- vector('list', nrow(lines.sldf))
  for (i in 1:nrow(lines.sldf)) {
    df <- as.data.frame(lines.sldf@lines[[i]]@Lines[[1]]@coords)
    names(df) <- c('x', 'y')
    df$lid <- lines.sldf@data[i,lidname]
    line.coord.ls[[i]] <- df
  }
  #line.df <- do.call('rbind', line.coord.ls)
  line.df <- rbind.fill(line.coord.ls)
  names(line.df)[3] <- lidname
  
  ## Add vertex IDs to line.df
  line.df$order <- 1:nrow(line.df)
  line.dt <- data.table(line.df)
  vertex.dt <- data.table(vid=vertex.df$vid, x=vertex.coord.mat[,1], y=vertex.coord.mat[,2])
  line.dt <- merge(line.dt, vertex.dt, by=c('x', 'y'), all.x=TRUE, all.y=FALSE)
  line.df <- as.data.frame(line.dt)
  line.df <- line.df[order(line.df$order),]
  
  return(line.df)
}


###############################################################
# REARRANGE LINES IN LINE.DF TO HYPERLINES
# Hyperline: Line that ends only at intersection or end point
# Hyperlines identified with new hlid identifier
###############################################################

make_hyperlines <- function(line.df) {
  ## Identify intersections and end points (= break points)
  line.df$first <- line.df$lid != c(-99,line.df$lid[-nrow(line.df)])
  line.df$last <- line.df$lid != c(line.df$lid[-1], -99)
  break.df <- line.df[line.df$first | line.df$last,]
  break.dt <- data.table(break.df)
  break.dt$unit <- 1
  break.dt <- break.dt[,list(cnt=sum(unit)), by='vid']
  break.df <- data.frame(break.dt)
  break.df <- break.df[break.df$cnt != 2,]
  break.df$end <- break.df$cnt==1
  break.df$intersection <- break.df$cnt > 2
  break.df <- break.df[,c('vid', 'end', 'intersection')]
  
  ## Add break info to line.df
  line.df <- merge(line.df, break.df, by="vid", all.x=TRUE, all.y=FALSE)
  line.df <- line.df[order(line.df$order),]
  line.df$end[is.na(line.df$end)] <- FALSE
  line.df$intersection[is.na(line.df$intersection)] <- FALSE
  line.df$brk <- ifelse(line.df$end | line.df$intersection, TRUE, FALSE)
  
  ## Iteratively build new lines
  new.line.ls <- list()
  fl.df <- line.df[line.df$first | line.df$last,]
  fl.df$hlid <- NA
  hlid <- 1
  fl.df$used <- FALSE
  start.vid <- (fl.df$vid[fl.df$brk])[1]
  start.order <- (fl.df$order[fl.df$brk])[1]
  this.lid <- (fl.df$lid[fl.df$brk])[1]
  while(any(!fl.df$used)) {
    
    fl.df$hlid[fl.df$lid==this.lid] <- hlid
    fl.df$used[fl.df$lid==this.lid] <- TRUE
    
    end.row <- fl.df[fl.df$lid == this.lid & fl.df$order != start.order,]
    end.vid <- end.row$vid
    
    is.reverse <- end.row$first
    this.line.df <- line.df[line.df$lid==this.lid,]
    this.line.df$hlid <- hlid
    if (is.reverse) {
      this.line.df <- this.line.df[nrow(this.line.df):1,]
    }
    new.line.ls[[length(new.line.ls)+1]] <- this.line.df
    
    if (end.row$brk) {
      start.row <- (fl.df[fl.df$used == FALSE & fl.df$brk,,drop=FALSE])
      if (nrow(start.row) == 0) {
        break
      } else {
        start.row <- start.row[1,]
      }
      this.lid <- start.row$lid
      start.vid <- start.row$vid
      start.order <- start.row$order
      hlid <- hlid + 1
    } else {
      start.row <- fl.df[fl.df$lid != this.lid & fl.df$vid == end.vid,,drop=FALSE]
      start.row <- start.row[1,]
      start.vid <- start.row$vid
      start.order <- start.row$order
      this.lid <- start.row$lid
    }
  }
  #new.line.df <- do.call('rbind', new.line.ls)
  new.line.df <- rbind.fill(new.line.ls)
  
  ## Remove duplicate verteces
  drop <- rep(FALSE, nrow(new.line.df))
  for (i in 2:(nrow(new.line.df))) {
    last.hlid <- new.line.df$hlid[i-1]
    this.hlid <- new.line.df$hlid[i]
    last.vid <- new.line.df$vid[i-1]
    this.vid <- new.line.df$vid[i]
    if (this.vid == last.vid & this.hlid == last.hlid) {
      drop[i] <- TRUE
    }
  }
  new.line.df <- new.line.df[!drop,]
  new.line.df$order <- 1:nrow(new.line.df)
  
  return(new.line.df)
}

###############################################################
# MAKE SLDF FROM LINE DATA FRAME
###############################################################

make_sldf_from_line_df <- function(line.df, lidname) {
  lid.vec <- unique(line.df[,lidname])
  Lines.ls <- vector('list', length(lid.vec))
  for (i in 1:length(lid.vec)) {
    this.Lines <- Lines(list(Line(line.df[line.df[,lidname]==lid.vec[i],c('x', 'y')])), i)
    Lines.ls[[i]] <- this.Lines
  }
  new.lines.sl <- SpatialLines(Lines.ls)
  new.lines.df <- data.frame(lid.vec)
  names(new.lines.df) <- lidname
  new.lines.sldf <- SpatialLinesDataFrame(new.lines.sl, data=new.lines.df, FALSE)
  return(new.lines.sldf)
}

###############################################################
# MAKE HYPERLINE SLDF
###############################################################

make_hyperline_sldf <- function(line.sldf) {
  
  line.sldf$lid <- 1:nrow(line.sldf)
  line.vertex.spdf <- make_vertex_tables(lines.sldf = line.sldf, vid.prefix = "v")
  line.df <- make_line_df(lines.sldf = line.sldf, vertex.spdf = line.vertex.spdf, lidname = "lid")
  hline.df <- make_hyperlines(line.df = line.df)
  hline.sldf <- make_sldf_from_line_df(line.df = hline.df, lidname = "hlid")
  
  return(hline.sldf)
}



