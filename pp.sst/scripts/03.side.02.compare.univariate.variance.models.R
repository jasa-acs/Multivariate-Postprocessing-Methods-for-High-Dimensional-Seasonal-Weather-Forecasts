
########################################################################################

#############  side script 3.2 - comparing univariate variance methods  ################

########################################################################################

# This script compares the univariate modelling of variance by moving averages to the 
# most commonly used methods in NGR: this is c^2 + d^2 SD, where SD is the ensemble spread.



##### setting up ######

rm(list = ls())

options(max.print = 1e3)

library(pp.sst)

name_abbr = "Full" 

save_dir = file.path('~','SST','Derived', name_abbr)

load(file = paste0(save_dir,"setup.RData"))


time_s32 = proc.time()

###### Estimate SD in NGR fashion: #####

DT = var_est_NGR_bm(DT, months = months, validation_years = validation_years)

DT = var_est_NGR_bl(DT, months = months, validation_years = validation_years,mc.cores = mc_cores)

DT = var_est_NGR_bb(DT, months = months, validation_years = validation_years,mc.cores = mc_cores)

# get CRPS for all methods:

CRPS_bm = DT[year %in% validation_years,mean(crps_na_rm(SST_bar,SST_hat,SD_hat_lr_bm),na.rm = TRUE),.(year,month)]
CRPS_bl = DT[year %in% validation_years,mean(crps_na_rm(SST_bar,SST_hat,SD_hat_lr_bl),na.rm = TRUE),.(year,month)]
CRPS_bb = DT[year %in% validation_years,mean(crps_na_rm(SST_bar,SST_hat,SD_hat_lr_bb),na.rm = TRUE),.(year,month)]
CRPS_sma = DT[year %in% validation_years,mean(crps_na_rm(SST_bar,SST_hat,SD_hat_sma),na.rm = TRUE),.(year,month)]
CRPS_ema = DT[year %in% validation_years,mean(crps_na_rm(SST_bar,SST_hat,SD_hat_ema),na.rm = TRUE),.(year,month)]

CRPS_comparison = data.table(CRPS_bm,CRPS_bl[,3],CRPS_bb[,3],CRPS_sma[,3],CRPS_ema[,3])
setnames(CRPS_comparison, c('year','month','lr_bm','lr_bl','lr_bb','sma','ema'))

save(CRPS_comparison,file = paste0(save_dir,'univ_scores_comp_NGR.RData'))


#### save stuff ####

time_s32 = proc.time() - time_s32

save.image(file = paste0(save_dir,"setup.RData"))


