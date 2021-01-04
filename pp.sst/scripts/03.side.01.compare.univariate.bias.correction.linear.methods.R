###############################################################################################

#############  side script 3.1 - comparing univariate bias correction methods  ################

###############################################################################################

# This script compares the univariate modelling of bias and variance by moving averages to the 
# most commonly used methods in Postprocessing: Three variants of NGR, as well as locally adaptive NGR. 



##### setting up ######

rm(list = ls())

options(max.print = 1e3)

library(pp.sst)

name_abbr = "Full" 

save_dir = file.path('~','SST','Derived', name_abbr)

load(file = paste0(save_dir,"setup.RData"))

###### compare to linear regression models ######


# get estimates, grouped by month, location and both
DT = bias_lr_bm(DT,months = months,validation_years = validation_years)
DT = bias_lr_bl(DT,months = months,validation_years = validation_years)
DT = bias_lr_bb(DT,months = months,validation_years = validation_years)

# get locally adaptive estimates 

DT= compute_clim(DT)

DT = bias_lr_bm_la(DT,months = months,validation_years = validation_years)

#### get MSEs ####

mod_vec = c('lr_bm','lr_bl','lr_bb','lr_bm_la','sma','ema')

MSE_linear_models = DT[year %in% validation_years,mean((T_hat_lr_bm-SST_bar)^2, na.rm = TRUE),.(year,month)]
setnames(MSE_linear_models, old = 3, new = 'MSE_lr_bm')

for(mean_mod in mod_vec[2:length(mod_vec)])
{
  var_name = paste0('T_hat_',mean_mod)
  if(mean_mod %in% c('sma','ema')) var_name = paste0('SST_hat_',mean_mod)
  
  MSE_mod0 = DT[year %in% validation_years,mean((get(var_name)-SST_bar)^2, na.rm = TRUE),.(year,month)]
  MSE_linear_models = data.table(MSE_linear_models,MSE_mod0[,3])
  setnames(MSE_linear_models, old = ncol(MSE_linear_models), new = paste0('MSE_',mean_mod))
}



MSE_compact = MSE_linear_models[,lapply(.SD,mean),.SDcols = 3:ncol(MSE_linear_models),by = month]

#### save stuff ####

time_s31 = proc.time() - time_s31

save(MSE_linear_models,MSE_compact,file = paste0(save_dir,'bc_scores_comp_NGR.RData'))

save.image(file = paste0(save_dir,"setup.RData"))
