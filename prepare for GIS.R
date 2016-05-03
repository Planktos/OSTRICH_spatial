
library("plyr")
library("dplyr")
library("stringr")
library("reshape2")
library("pastecs")

# load transect data
load("data/3w/biophys merge/all_joined_3w_shallow_bin.Robj")

d <- binned_biophy

d <- d[,c(2:length(d))]

transect <- "3W-shallow"

# Marked Points
d$mark <- ifelse(d$group == "fish", 1, 0)

d$depth <- d$depth*-1

d$dateTime <- round(d$dateTime, units = c("secs"))

write.csv(d, file = paste0("data/3w/GIS input/",transect,"_fish-marked.csv"), row.names = F)

