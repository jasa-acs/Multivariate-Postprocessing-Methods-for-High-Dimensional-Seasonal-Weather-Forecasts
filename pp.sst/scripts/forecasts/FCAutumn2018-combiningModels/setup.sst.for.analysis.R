#############################################################################################

###### master script part 1 - setting up and creating or loading of wide data table #########

#############################################################################################

# This script sets up for a full post-processing analysis 
# for a lon/lat window to be specified below

rm(list = ls())

# get SST Data of GCFS1:

setwd("~/PostClimDataNoBackup/SFE/GCFS1")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

load(file = "./DT_SST.RData")

# get the NorCPM forecast

setwd("~/NR/SFE")

dt = load_combined_wide()

#only keep what we also have as NorCPM forecast

dt[,YM:= year * 12 + month]
DT[,YM:= year * 12 + month]

ym = unique(dt[,YM])
DT = DT[YM %in% ym,]

# get SSTs, there is something weird happening in 03/2018 so we leave it out

SST_1 = dt[month %in% c(1:3,7:12) & YM < 24219,SST_bar]
SST_2 = dt[month %in% c(1:3,7:12) & YM > 24219,SST_bar]

#check:
length(SST_1)+ length(SST_2) - DT[,.N] # -64800, the missing year

#fill in:
DT[YM < 24219, SST_bar := SST_1]
DT[YM > 24219, SST_bar := SST_2]

#get grid_ids

gids = rep(dt[year == 1985 & month == 1, grid_id],length(unique(DT[,YM])))
DT[,grid_id:= gids]


name_abbr = "Aut_2018_small/GCFS1" # for northern atlantic ocean

lat_box = c(50,85)
lon_box = c(-20,40)


ens_size = 15 # size of forecast ensemble

training_years = 1985:2013
validation_years = 2014:2017 # all previous years are used for training 
months = 7:12

DT = DT[month %in% months,]
DT = DT[Lon >= lon_box[1] & Lon <= lon_box[2] & Lat >= lat_box[1] & Lat <= lat_box[2],]


# create directories

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")
dir.create(save_dir, showWarnings = FALSE)

plot_dir = paste0("./figures/", name_abbr,"/")
dir.create(plot_dir, showWarnings = FALSE)

rm(dt,gids,SST_1,SST_2,ym)

save.image(file = paste0(save_dir,"setup.RData"))

