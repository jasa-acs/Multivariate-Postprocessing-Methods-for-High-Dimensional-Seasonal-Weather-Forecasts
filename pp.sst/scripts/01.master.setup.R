

#############################################################################################

###### master script part 1 - setting up and creating or loading of wide data table #########

#############################################################################################

# This script sets up a run of the post-processing analysis 

# installing the R-package:

rm(list = ls())

setwd("~/pkg/paper/PostprocessingSST/")
options(max.print = 1e3)

#install.packages('devtools')
library(devtools)

devtools::install_github('ClaudioHeinrich/pp.sst')


library(pp.sst)

#start timer:
time_s1 = proc.time()


# choose an abbreviation for this run and give it a description, see the README file for more details. 

name_abbr = "pp.sst/Full"

description = 'Working on the full data set.'

#### specify your directories: ###

# Directory for derived datasets, this should change when you change name_abbr
save_dir = file.path('~','SST','Derived', name_abbr)
dir.create(save_dir,recursive = TRUE,showWarnings = FALSE)

# Directory for plots, this should change when you change name_abbr
plot_dir = file.path('~','SST','Figures', name_abbr)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)


# choose the area of the globe to consider: Reduce this to a smaller window when testing scripts.

lat_box = c(0, 65)
lon_box = c(-90, 40)


 # a couple of parameters:

ens_size = 9 # size of forecast ensemble

training_years = 1985:2000
validation_years = 2001:2016 

months = 1:2

mc_cores = 6 # number of cores for parallelization

### subset data set ###

# note that some example data DT is included in the package, check the documentation by typing ?DT. 
# If you are interested in getting access to the full dataset considered in the paper please contact the authors.

DT = DT[Lon >= lon_box[1] & Lon <= lon_box[2] & Lat >= lat_box[1] & Lat <= lat_box[2]][month %in% months][year %in% c(training_years,validation_years)]

# tidy up DT:

setkey(x = DT,year,month,Lon,Lat)
DT[, YM := 12*year + month]
DT = DT[order(year,month,Lon,Lat)]

setcolorder(DT,c("year","month",'Lon','Lat','YM','grid_id','SST_bar','Ens_bar','Ens_sd'))

#### time, update script counter, save ####

time_s1 = proc.time() - time_s1

script_counter = 1

save.image(file = paste0(save_dir,"setup.RData"))

