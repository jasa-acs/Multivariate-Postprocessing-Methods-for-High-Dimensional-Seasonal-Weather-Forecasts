
##########################################################################################

###################  master script part 5 - multivariate modelling  ######################

##########################################################################################

# This script generates and saves forecasts for all methods.

#### setting up ######

rm(list = ls())

options(max.print = 1e3)

library(pp.sst)

name_abbr = "pp.sst/Full" 

save_dir = file.path('~','SST','Derived', name_abbr)

load(file = paste0(save_dir,"setup.RData"))

#make reproducible:
set.seed(101010)

time_s5 = proc.time()

# specifications for the desired forecasts:

fc_years = validation_years
fc_months = months
fc_ens_size = 50 

mod_vec = c('PCA_mc','PCA_ac','GS','ECC','Schaake')


###################################################
###################### PCA  #######################
###################################################

PCA_fc_mc = forecast_PCA_mult_corr(DT, 
                                   Y = fc_years,
                                   M = fc_months,
                                   n = fc_ens_size,
                                   nPCs = nPCs,
                                   cov_dir = PCA_dir)

save(PCA_fc_mc,file = paste0(PCA_dir,"fc_mc.RData"))

# When fc_ens_size is chosen large, the forecast datatables get large, so let's make room...
rm(PCA_fc_mc)
gc()

PCA_fc_ac = forecast_PCA_add_corr(DT, 
                                  Y = fc_years,
                                  M = fc_months,
                                  n = fc_ens_size,
                                  nPCs = nPCs,
                                  cov_dir = PCA_dir)

save(PCA_fc_ac,file = paste0(PCA_dir,"fc_ac.RData"))

rm(PCA_fc_ac)
gc()

###################################################
################## geostationary ##################
###################################################

GS_fc = forecast_GS(DT,
                    Y = validation_years,
                    M = months,
                    n = fc_ens_size,
                    var_dir = GS_dir,
                    mc_cores = mc_cores)

save(GS_fc,file = paste0(GS_dir,"fc.RData"))

rm(GS_fc)
gc()

########################################
################ ECC  ##################
########################################

ECC_fc = forecast_ECC(DT,
                      Y = validation_years,
                      M = months,
                      ens_size = ens_size)

save(ECC_fc,file = paste0(ECC_dir,"fc.RData"))
rm(ECC_fc)


####################################################
################ Schaake shuffle  ##################
####################################################

Schaake_fc = forecast_Schaake(DT,
                              obs_ens_years = training_years,
                              Y = validation_years,
                              M = months)

save(Schaake_fc,file = paste0(Schaake_dir,"fc.RData"))
rm(Schaake_fc)

#### time, update script counter, save ####

time_s5 = proc.time() - time_s5

script_counter = 5

save.image(file = paste0(save_dir,"setup.RData"))


