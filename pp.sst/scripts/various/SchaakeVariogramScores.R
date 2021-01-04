##### setting up ######

rm(list = ls())



setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "NAO/lv" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")


load(paste0(save_dir,'vs.RData'))

###### Schaake ######

pp = 0.5
validation_years = 2001:2016
months = 1:12
training_years = 1985:2000
mc_cores = 6

Schaake_dir = "~/PostClimDataNoBackup/SFE/Derived/NAO/Schaake/"

load(paste0(Schaake_dir,"fc.RData"))

vs_Schaake = var_sc_par(dt_fc = Schaake_fc,  p = pp,
                        years = validation_years, ms = months, n = length(training_years), 
                        save_dir = NULL,
                        mc_cores = mc_cores)

###########

# combine, boxplot and save

vs_dt = data.table( PCA_mc = vs_PCA_mc[,mean(vs)],
                    PCA_ac = vs_PCA_ac[,mean(vs)],
                    GS = vs_GS[,mean(vs)], 
                    ECC = vs_ECC[,mean(vs)],
                    Schaake = vs_Schaake[,mean(vs)])


save(vs_dt,vs_PCA_ac,vs_PCA_mc,vs_GS,vs_ECC,vs_Schaake, file = '~/PostClimDataNoBackup/SFE/Derived/NAO/vs_backup.RData') # save here in order to not get an overwrite when a current job finishes.


