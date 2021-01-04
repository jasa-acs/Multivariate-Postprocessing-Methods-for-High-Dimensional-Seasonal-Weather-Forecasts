

##### setting up ######

rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "Full/lv" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))



# get locally adaptive estimate 

DT= compute_clim(DT)
DT = bias_lr_bm_la(DT,validation_years = validation_years)


save.image(paste0(save_dir,'setup.RData'))