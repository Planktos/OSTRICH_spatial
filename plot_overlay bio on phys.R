
#PURPOSE: Plot physical data for each transect

# modified from "plot_physical_data_2014.R"
# late modified: 3 May 2016 Kelly Robinson
#------

#libraries
library("plyr")
library("stringr")
library("ggplot2")
library("reshape2")
library("grid")
library("gridExtra")
library("scales")
library("akima")
library("dplyr")
library("oce")
library("pastecs")

# setup R to keep decimal seconds in the times
options("digits.secs"=3)

#Call functions:
#-----------------
# Compute the straight line distance (km) from the starting point of a lat,lon trajectory
  dist.from.start <- function(lat, lon) {
    library("oce")
    geodDist(lat1=lat, lon1=lon, lat2=na.omit(lat)[1], lon2=na.omit(lon)[1])/1.852 # use if you want to convert from km to nautical miles
  }

# Spectral colour map from ColorBrewer
  spectral <- function(n=6) {
    library("RColorBrewer")
    rev(brewer.pal(name="Spectral", n=n))
  }
  
  scale_fill_spectral <- function(...) {
    scale_fill_gradientn(colours=spectral(...))
  }
  scale_colour_spectral <- function(...) {
    scale_colour_gradientn(colours=spectral(...))
  }

# Interpolate a slice of data for which the x-axis is a distance in nautical miles
interp.dist <- function(x, y, z, anisotropy=1000, x.step=500, y.step=2.5, smooth=FALSE, theta=0.2, ...) {
    #
    # Interpolate data over a distance coordinate
    #
    # x   vector of distance *IN KILOMETERS*
    # y   vector of depth in m
    # z   vector of measured variable
    # anisotropy  anisotropy ratio between x and y
    # x/y.step    interpolation grid steps in m
    # smooth      boolean, wether to smooth the first interpolation using fields::image.smooth
    # x/y.step.smooth   interpolation grid step for the smoothing
    # grid.smooth intepolation grid for the smoothing, overrides x/y.step.smooth
    # theta       bandwidth for the kernel smoother in fields::image.smooth
  
    library("akima")
    library("reshape2")
  
    # correct x-axis for anisotropy between horizontal and vertical
    x <- x <- x*1852/anisotropy #if x unit is nautical miles
  
    # interpolate
    i <- interp(x=x, y=y, z=z, xo=seq(0, max(x), by=x.step/anisotropy), yo=seq(0, max(y), by=y.step), ...)
  
    # smooth
    if ( smooth ) {
      library("fields")
      i <- image.smooth(i, grid=list(x=i$x, y=i$y), theta=theta)
    }
  
    # extract a data.frame
    out <- melt(i$z, varnames=c("x","y"))
    out$x <- i$x[out$x] * anisotropy/1852
    out$y <- i$y[out$y]
  
    return(out)
  }

## Detect up and down casts in a depth yo
  detect.casts <- function(depth, order=200) {
    # smoothing the depth profile using a moving average and find the turning points
    
    # smooth depths
    library("pastecs")
    depth_avg <- decaverage(-depth, times=3, weights=c(seq(1, order), order+1, seq(order, 1, -1)))
    # plot(depth_avg)
    depth_avg <- as.numeric(pastecs::extract(depth_avg, component="filtered"))
    
    # detect turning points
    TP <- suppressWarnings(turnpoints(depth_avg))
    
    # set cast numbers (different for up and down casts)
    cast <- cumsum(TP$peaks | TP$pits) + 1
    
    # detect which are up and which are down casts:
    # if the first turning point is a peak, then the first cast (and all odd casts) are upcasts
    if ( TP$firstispeak ) {
      # these are the types for
      #              even  & odd   cast numbers
      castTypes <- c("down", "up")
    } else {
      castTypes <- c("up", "down")
    }
    down.up <- castTypes[cast %% 2 + 1]
    
    return(data.frame(cast, down.up))
  }

#rowShift
rowShift <- function(x, shiftLen = 1L) {
  r <- (1L + shiftLen):(length(x) + shiftLen)
  r[r<1] <- NA
  return(x[r])
}

#get legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
  #-------------------

p.vars <- c("temp", "salinity", "fluoro", "sw.density", "oxygen", "irradiance")

b.vars <- c("fish", "copepod calanoid", "copepod oithona", "appendicularian")

# select bin size
bin.size <- 1
density.metric <- paste0("density.per.",bin.size,".m3")

# list merged biophys data
files <- list.files(path = "data/3w/biophys merge", pattern = paste0("bin_",bin.size,"m3.Robj"))


plot.densities <- function(x){

binned_biophy$group <- as.character(binned_biophy$group)

data <- binned_biophy[,c("depth", "lat", "lon", "temp", "salinity", "pressure", "fluoro", "oxygen", "irradiance", "dateTime", "group","density.per.1.m3")]

# physical data ----

  # calculate seawater density
  data$sw.density <- swRho(salinity = data$salinity, temperature = data$temp, pressure = data$pressure, eos = "unesco")
  data$sw.density <- data$sw.density-1000
  
  # interpolate physical data
  data$distanceFromStart <- dist.from.start(data$lat, data$lon)
  
  # detect yos
  if(tow == "und"){
      
      casts <- detect.casts(data$depth, order=50)
      data <- cbind(data, casts)
      
      dm <- melt(data, id.vars=c("depth", "down.up", "distanceFromStart"), measure.vars=p.vars)
      
      di <- ddply(dm, ~variable, function(x) {
        x <- na.omit(x[which(x$down.up=="up"),])
        xi <- interp.dist(x=x$distanceFromStart, y=x$depth, z=x$value, duplicate="mean", x.step=300, y.step=1, anisotropy=1000)
      })
  
  } else {
      
      dm <- melt(data, id.vars=c("depth", "distanceFromStart"), measure.vars=p.vars)
      
        di <- ddply(dm, ~variable, function(x) {
        xi <- interp.dist(x=x$distanceFromStart, y=x$depth, z=x$value, duplicate="mean", x.step=300, y.step=1, anisotropy=1000)
      })
  }
  
  names(di)[names(di)=="x"] <- "distance"
  names(di)[names(di)=="y"] <- "depth"
  
  di$distance <- di$distance*1.852 #convert from nautical miles to kilometers

# create physical data section plots ----
  
  # set maximum & minimum depths for tow types
  if(tow == "und"){
    max.depth <- -100
    min.depth <- 0
  } else {
    max.depth <- floor(max(data$depth)*-1)
    min.depth <- ceiling(min(data$depth)*-1)
  }
  
  # physical data plots
  t <- ggplot(di[di$variable == "temp",], aes(x=distance, y=-depth)) +
    geom_tile(aes(fill=value), na.rm=T) +
    stat_contour(aes(z=value), colour="black", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
    scale_fill_gradientn(colours = spectral(), limits = c(12, 29), na.value = NA, name = "°C") +
    scale_color_gradient(high = spectral()) +
    scale_x_continuous(expand=c(0.02,0.02), "distance (km)", breaks = seq(1,max(di$distance),1)) + 
    scale_y_continuous(expand=c(0,0), "depth (m)", limits = c(max.depth,min.depth)) +
    theme(strip.text.y = element_text(size = 10)) +
    theme(axis.text = element_text(size = 10)) + theme(axis.title = element_text(size = 12))
  
   
  s <- ggplot(di[di$variable == "salinity",], aes(x=distance, y=-depth)) +
    geom_tile(aes(fill=value), na.rm=T) +
    stat_contour(aes(z=value), colour="black", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
    scale_fill_gradientn(colours = spectral(), limits = c(34,37), na.value=NA, name = "sal(ppt)") +
    scale_color_gradient(high = spectral()) +
    scale_x_continuous("distance (km)", expand=c(0.02,0.02), breaks = seq(1,max(di$distance),1)) +
    scale_y_continuous(expand=c(0,0), "depth (m)", limits = c(max.depth,min.depth)) +
    theme(strip.text.y = element_text(size = 10)) +
    theme(axis.text = element_text(size = 10)) + theme(axis.title = element_text(size = 12))
  
  
  d <- ggplot(di[di$variable == "sw.density",], aes(x=distance, y=-depth)) +
    geom_tile(aes(fill=value), na.rm=T) +
    stat_contour(aes(z=value), colour="black", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
    scale_fill_gradientn(colours = spectral(), na.value=NA, limits = c(22,28), name = "d(kg/m3") +
    scale_color_gradient(high = spectral()) +
    scale_x_continuous("distance (km)", expand=c(0.02,0.02),breaks = seq(1,max(di$distance),1)) +
    scale_y_continuous(expand=c(0,0), "depth (m)", limits = c(max.depth,min.depth))  +
    theme(strip.text.y = element_text(size = 10)) +
    theme(axis.text = element_text(size = 10)) + theme(axis.title = element_text(size = 12))
  
  f <- ggplot(di[di$variable == "fluoro",], aes(x=distance, y=-depth)) +
    geom_tile(aes(fill=value), na.rm=T) +
    stat_contour(aes(z=value), colour="black", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
    scale_fill_gradientn(colours = spectral(), na.value=NA,limits = c(0.09, 0.7), name = "fluoro(v))", oob = squish) +
    scale_color_gradient(high = spectral()) +
    scale_x_continuous("distance (km)", expand=c(0.02,0.02), breaks = seq(1,max(di$distance),1))+
    scale_y_continuous(expand=c(0,0), "depth (m)", limits = c(max.depth,min.depth)) +
    theme(strip.text.y = element_text(size = 10)) +
    theme(axis.text = element_text(size = 10)) + theme(axis.title = element_text(size = 12))
  
    
  o <- ggplot(di[di$variable == "oxygen",], aes(x=distance, y=-depth)) +
    geom_tile(aes(fill=value), na.rm=T) +
    stat_contour(aes(z=value), colour="black", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
    scale_fill_gradientn(colours = spectral(), na.value=NA, limits = c(0,5.0), name = "dO2(ml/l)") +
    scale_color_gradient(high = spectral()) +
    scale_x_continuous("distance (km)", expand=c(0.02,0.02),breaks = seq(1,max(di$distance),1)) +
    scale_y_continuous(expand=c(0,0), "depth (m)", limits = c(max.depth,min.depth)) +
    theme(strip.text.y = element_text(size = 10)) +
    theme(axis.text = element_text(size = 10)) + theme(axis.title = element_text(size = 12))
    

# biological data -----
  
  b <- dcast(data, dateTime + lat + lon + depth + distanceFromStart ~ group, value.var = density.metric)
  
  db <- melt(b, id.vars=c("depth","distanceFromStart"), measure.vars=b.vars)
  
  names(db)[names(db)=="value"] <- "density"
  names(db)[names(db)=="variable"] <- "group"
  names(db)[names(db)=="distanceFromStart"] <- "distance"
  
  db$distance <- db$distance*1.852 #convert from nautical miles to kilometers
  
  
# overlay plots -----
  
 ddply(.data = db, .variables = "group", function(x){
   
  taxa <- ggplot(data = x[x$density > 0,]) + geom_point(aes(x = distance, y = -depth, size = density), alpha=0.5) +
          scale_size_area(x$density, max_size=8,guide = guide_legend(title = density.metric, direction ="horizontal", title.position = "left", nrow = 1)) + labs(title=(paste0(transect," ",x$group))) + theme(legend.position="bottom")
   
  # get legend
  density.guide <- g_legend(taxa)
   
# temperature 
  
    taxa.t <- t + geom_point(data = x[x$density > 0,], aes(x = distance, y = -depth, size = density), alpha=0.1) +
              scale_size_area(density, max_size=6, guide = "none") +  ggtitle(label = (paste0(transect,": ",x$group," Indiv/",bin.size,"m3"))) + theme(plot.title = element_text(lineheight=1.0, face="bold"))
            
# salinity
    taxa.s <- s + geom_point(data = x[x$density > 0,], aes(x = distance, y = -depth, size = density), alpha=0.1) +
              scale_size_area(density, max_size=6, guide = "none")

# density
    taxa.d <- d + geom_point(data = x[x$density > 0,], aes(x = distance, y = -depth, size = density), alpha=0.1) +
                scale_size_area(density, max_size=6, guide = "none")
#fluro
    taxa.f <- f + geom_point(data = x[x$density > 0,], aes(x = distance, y = -depth, size = density), alpha=0.1) +
      scale_size_area(density, max_size=6, guide = "none")
     
# OXYGEN
    taxa.o <- o + geom_point(data = x[x$density > 0,], aes(x = distance, y = -depth, size = density), alpha=0.1) +
      scale_size_area(density, max_size=6, guide = "none")
   
g <- grid.arrange(taxa.t, taxa.s, taxa.d, taxa.f, taxa.o, density.guide, ncol=1)

#print to output directory
    
    png(file = paste0(output.dir,"/",transect,"_",x$group,".png"), width = 9, height = 14, units = "in", res = 300)
    plot(g) #change for each taxon
    dev.off()
     
  }, .progress = "text")
  
}

#loop through all binned data files
  
for (i in 1:length(files)){
    
    load(paste0("data/3w/biophys merge/",files[i]))
  
    # identify region, tow & transect
    temp <- str_split_fixed(string = files[i], pattern = "_", n = 6)
    region <- str_sub(temp[,3])
    tow <- str_sub(temp[,4])
    transect <- str_c(region,tow,sep = "_")
    rm(temp)

    # prepare storage
    output.dir <- paste("plots/",region,"/distributions/",tow,"_",bin.size,"m3", sep="")
    dir.create(output.dir)
    
    plot.densities(files[i])
    }
  
