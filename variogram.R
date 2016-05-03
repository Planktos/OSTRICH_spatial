
options(digit.secs = 3)

library(geoR)
library(plyr)
library(dplyr)

# list merged biophys data
files <- list.files(path = "data/3w/biophys merge", pattern = paste0("bin_",bin.size,"m3.Robj"))

# select bin size
bin.size <- 0.25
# select lag size
lags <- 10
density.metric <- paste0("indiv.per.",bin.size,"m3")

#output directory
output.dir = paste0("plots/",region,"/variograms_",lags," lags/")
dir.create(output.dir)

plot.variogram <- function(x){

  d <- binned_biophy
  
  data <- d[,c(2:length(d))]
  
#variogram
  v <- data[,c("group", density.metric, "lon", "lat")]
  
  ddply(.data = v, .variables = "group", function(x){
    
    x.geo <- as.geodata(x, coords.col = 3:4, data.col = 2)
    
    group.var <- variog(x.geo, estimator.type = "classical", max.dist = 300)
    
    group.var.summary <- cbind(c(1:length(group.var$v)), group.var$v, group.var$n)
    colnames(group.var.summary) <- c("lag", "semi-variance", "# of pairs")
    
    title <- paste0(transect,"_",(unique(x$group)))
    
    #print to directory
    file = paste0(output.dir,title,".png")
    
    png(file=file, width = 11, height = 8.5, units = "in", res = 150)
    par(mfrow = c(2, 2))
    plot(group.var, main = title) #change for each taxon
    dev.off()
  
  }, .progress = "text")

}

for (i in 1:length(files)){
  
  load(paste0("data/3w/biophys merge/",files[i]))
  
  # identify region, tow & transect
  temp <- str_split_fixed(string = files[i], pattern = "_", n = 6)
  region <- str_sub(temp[,3])
  tow <- str_sub(temp[,4])
  transect <- str_c(region,tow,sep = "_")
  rm(temp)
  
  # prepare storage
  output.dir <- paste("plots/",region,"/variograms_",lags,"/",transect,"_",bin.size,"m3/", sep="")
  dir.create(output.dir)
  
  plot.variogram(files[i])
}

