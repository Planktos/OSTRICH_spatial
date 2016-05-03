#
# Calculate taxa densities (#/m3) for OSTRICH REU
#
## The purpose of this script is to estimate the taxa densities (individuals per m3) for OSTRICH 2014.

# 1. read in individual transect data (output from 'merging_phy_bio_20151125.R")
# 2. estimate the average number of seconds it took to image 1-m3 of water in 1E and 1W
# 3. create time bins
# 4. calculate the total number of each taxon in a given time bin (each representing 1-m3), so individuals/m3

# Created by: Kelly Robinson
# Created on: 5 August 2015
# Date last modified: 28 April 2016
#
#---------------------------------------------------------------------------

library("plyr")
library("stringr")
library("reshape2")
library("pastecs")

# setup to keep decimal seconds in the times
options("digits.secs"=3)

# load merged biophys data
files <- list.files(path = "data/3w/biophys merge", pattern = "_nobin.Robj")

# set bin size
bin.size = 1 # OPTIONS: '1' for 1-m3; '0.5' for 0.5-m3

bin.data <- function(x){

  x <- bio.f
  d <- x

# bin taxon counts into 0.5-m3 bins using time it takes to image 1-m3 (4.6 sec)
  #seconds it took to image 0.5-m3 of water ('see 'lib_process.R" scripts for calculations of 'm3.sec' field)
  
  binSize <-  ((1/mean(d$m3.sec, na.rm = TRUE))*bin.size)
  
  #calculate the number of bins
  d.maxT <- max(d$dateTime, na.rm=T)
  d.minT <- min(d$dateTime, na.rm=T)
  d.bins=seq(d.minT, d.maxT, by=binSize)
  d$timeBin <- cut(d$dateTime, breaks=d.bins, labels=1:(length(d.bins)-1))

  # convert from factor to numeric
  d$timeBin <- as.numeric(d$timeBin)
  
  # binning
    # find total counts of each taxa in a given time bin, which represents the number of seconds it took to image 1 m3 of water.
    d.count.bin <- aggregate(count~timeBin+group, data = d, FUN = sum)
    
    # bin the dateTimes for each time bin & taxon combination
    d.timeBin <- aggregate(dateTime~timeBin+group, data = d, FUN = mean)
    
    # bin the physical data for each time bin & taxon combination
    d.phyBin <- aggregate(cbind(depth,temp,salinity,pressure,fluoro,oxygen,irradiance,lat,lon)~timeBin+group, data = d, FUN = mean)
  
    # cHECK for time off-set in d.timeBin against minimum and maximum times in 'bio'
#       summary(d.timeBin)
#       summary(bio.f)
      
        # OPTIONAL: fix the time offset if present (e.g. if a shift from 11:26:00 to 8:26:00 happened)
          d.timeBin$dateTime <- d.timeBin$dateTime + 3*3600
          #check again to make sure the dateTime was corrected
          # summary(d.timeBin)
    
    #merge all three together (each should have the exact same number of rows)
    d.taxa <- merge(d.count.bin, d.timeBin, by = c("timeBin", "group"))
    d.taxa <- merge(d.taxa, d.phyBin, by = c("timeBin", "group"))
  
#rename 'count' field to density.m3 because now it represents the total number of organisms in 1 m3.
names(d.taxa)[names(d.taxa)=="count"] <- (str_c("indiv.per.", bin.size, "m3", sep = ""))
  
#order by time bin
binned_biophy <- d.taxa[order(d.taxa$timeBin),]

return(binned_biophy)

}

for (i in 1:length(files)){
  
  load(paste0("data/3w/biophys merge/",files[i]))
  
  # identify transect
  temp <- str_split_fixed(string = files[i], pattern = "_", n = 5)
  region <- str_sub(temp[,3])
  tow <- str_sub(temp[,4])
  transect <- str_c(region,tow,sep = "_")
  rm(temp)
  
  binned_biophy <- bin.data(files[i])
  
  binned_biophy$transect <- transect
  
  #save the data
  save(binned_biophy, file=paste0("data/3w/biophys merge/all_joined_",transect,"_bin_",bin.size,"m3.Robj"))
  }




