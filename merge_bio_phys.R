
library("plyr")
library("dplyr")
library("stringr")
library("reshape2")

# setup R to keep decimal seconds in the times
options("digits.secs"=3)

# Sparse ConvNets v2016 test predictions
#d <- read.csv(file = "data/3w/plankton_predictions_3w_shallow.csv", header = T, stringsAsFactors = F)
#d <- read.csv(file = "data/3w/plankton_predictions_3w_deep.csv", header = T, stringsAsFactors = F)
#d <- read.csv(file = "data/3w/plankton_predictions_3w_mid.csv", header = T, stringsAsFactors = F)
d <- read.csv(file = "data/3w/predictions/plankton_predictions_3w_und.csv", header = T, stringsAsFactors = F)

# identify transect
transect <- "3W_und"
transect.ID <- "OST14-3W-Und"


# Read processed physical data

  load("ost14_phy_t.R")
  
  # select transect
  phy <- subset(phy_t, transect.id == transect.ID)
  
  #grab fields of interest
  phy <- phy[,c("dateTime", "depth","temp","salinity","pressure","fluoro","oxygen","irradiance","m3.sec","lat","lon")]

  #get transect start and end times
  time.start <- as.character(min(phy$dateTime))
  time.end <- as.character(max(phy$dateTime))

# Read in classification predictions

    p <- as.data.frame(d)
    
    names(p)[names(p)=="image"] <- "path"
    temp <- str_split_fixed(string = p$path, pattern = "/", n = 4)
    p$image <- str_sub(temp[,4])
    rm(temp)
    p$path <- NULL
    
    #re-order column names
    p <- p %>% select(image, everything())
    
    #get class names from test predictions
    classes <- names(p)
    measured.variables <- classes[2:length(classes)]
    
    pm <- melt(data = p, measure.vars = measured.variables, id.vars = "image")
    
    #return the class with highest probability for each segment
    class.pred <- ddply(.data = pm, .variables = "image", function(x){
      d <- subset(x[x$value == max(x$value),])
      return(d)
    }, .progress="text")
    
    names(class.pred)[names(class.pred)=="variable"] <- "class"
    names(class.pred)[names(class.pred)=="value"] <- "probability"
    
    class.pred$probability <- as.character(class.pred$probability)
    class.pred$class <- as.character(class.pred$class)
    
    # add count for each image
    class.pred$count <- 1
    
    save(class.pred, file = paste0("data/3w/",transect,"_preds.Robj"))

    
# Get time-stamp from segment file names
    
    class.pred$dateTime <- as.POSIXct(str_sub(class.pred$image, 1, 14), format="%Y%m%d%H%M%OS", tz="America/New_York") + as.numeric(str_sub(class.pred$image, 16, 18))/1000
    
    # LARGE Camera: 430 frames per stack means 24.46 seconds between each stack (averaged over ~600 stacks)
    # if 430 frames then 0.05688609 sec per frame
    class.pred$dateTime <- class.pred$dateTime + as.numeric(str_sub(class.pred$image, 20, 23)) * 0.05688609
    
    
# assign classes to taxonomic groups
    train.groups <- read.csv("ostrich_training_set_v3_class_groups.csv", header = T)
    train.groups <- as.data.frame(train.groups)
    train.groups$class <- as.character(train.groups$class)
    train.groups$group <- as.character(train.groups$group)
    
    bio <- merge(x = class.pred, y = train.groups, by.x = "class", by.y = "class", all.x = T)
    
    # total of all same class organisms in a single dateTime instance
    bio <- dcast(bio, dateTime~group, sum, value.var="count")
    
    
# Add TRUE Zeros to bio frame:
    
    # 1) Get unique frame times so all true zeros are accounted for
      frames <- read.csv(file = "data/3w/3w_frames.csv", stringsAsFactors = F) # changes with each transect
      
      f <- as.data.frame(frames)
      colnames(f) <- c("frame")
      
      f$dT <- str_sub(f$frame, 1, 18)
      f$yy <- str_sub(f$dT, 1, 4)
      f$mm <- str_sub(f$dT, 5, 6)
      f$dd <- str_sub(f$dT, 7, 8)
      f$hh <- str_sub(f$dT, 9, 10)
      f$min <- str_sub(f$dT, 11, 12)
      f$sec <- str_sub(f$dT, 13, 18)
      f$date <- str_c(f$yy, f$mm, f$dd,  sep="-")
      f$time <- str_c(f$hh, f$min, f$sec, sep=":")
      f$date.time <- str_c(f$date, f$time, sep=" ")
    
    # 2) extract frames from transect
    
      #start time of undulation transect
      start <- subset(f, date.time > time.start)
      
      #end time of undulation transect
      transect.frames <- subset(start, date.time < time.end)
    
    # 3) get date-time of individual frames
      t.frames <- transect.frames[, c("frame")]
      t.frames <- as.data.frame(t.frames)
      colnames(t.frames) <- c("frame")
      t.frames$frame <- as.character(t.frames$frame)
      
    
      t.frames$stack <- as.POSIXct(str_sub(t.frames$frame, 1, 14), format="%Y%m%d%H%M%OS", tz="America/New_York") + as.numeric(str_sub(t.frames$frame, 16, 18))/1000
      t.frames$dateTime <- t.frames$stack + as.numeric(str_sub(t.frames$frame, 20, 23)) * 0.05688609 #get seconds
    
      f.unique <- as.data.frame(unique(t.frames$dateTime)) # get unique frame times
      colnames(f.unique) <- "dateTime"
    
      
    # 4) Merge the entire set of unique dateTimes with bio 'dateTimes'
    bio.f <- merge(x = f.unique, y = bio, by = "dateTime", all.x = TRUE)
    
    # 4) Set 'NA' counts associated with taxonomic group at dateTime to zero
    bio.f[,2:length(bio.f)][is.na(bio.f[,2:length(bio.f)])] <- 0
  

# Melt in prep for joining with physical data
    bio.f <- melt(bio.f, id.vars = "dateTime")
    
    names(bio.f)[names(bio.f)=="variable"] <- "group"
    names(bio.f)[names(bio.f)=="value"] <- "count"    


# Interpolate physical parameters for each biological sample
    
    # unique biological sample times
    times <- sort(unique(bio.f$dateTime))
    
    # prepare storage
    interp <- data.frame(dateTime=times)
    
    # interpolate all numeric columns linearly
    cols <- names(phy)[which(llply(phy, class)=="numeric")]
    for (col in cols) {
      interp[,col] <- approx(phy$dateTime, phy[,col], times)$y
    }
    
    # get transect, cast etc. using nearest neighbour interpolation
    closestIndex <- round(approx(phy$dateTime, 1:nrow(phy), times)$y)
    interp <- cbind(interp, phy[closestIndex,])
    
    # original JL code (August 2015)
    # put all that back in bio.f
    bio.f <- join(x = bio.f, y = interp, by = "dateTime")
    
    # just grab the data for which there exists physical data
    bio.f <- bio.f[complete.cases(bio.f),]
    
    #save bio data frame as an R object
    save(bio.f, file=paste0("data/3w/all_joined_",transect,"_nobin.Robj")) #name of file will change with each transect
    