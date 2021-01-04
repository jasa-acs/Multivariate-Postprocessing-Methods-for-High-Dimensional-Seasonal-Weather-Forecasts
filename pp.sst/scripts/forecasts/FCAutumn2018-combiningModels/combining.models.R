

##### setting up ######

rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)


name_abbr_1 = "Aut_2018_small" 
save_dir_1= paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr_1,"/")

name_abbr_2 = "Aut_2018_small/GCFS1" 
save_dir_2 = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr_2,"/")


### load data ###

load(file = paste0(save_dir_1,"setup.RData"))

DT_1 = DT
sc_ema_1 = sc_ema
sc_sma_1 = sc_sma
opt_par_1 = opt_par

# get RMSE of the model

if(opt_par_1[1] == "ema")
{
  RMSE_1 = sc_ema_1[a == opt_par_1[2],sqrt(MSE)]
} else{
  RMSE_1 = sc_sma_1[win_length == opt_par_1[2],sqrt(MSE)]
}


# get second model

load(file = paste0(save_dir_2,"setup.RData"))

DT_2 = DT
sc_ema_2 = sc_ema
sc_sma_2 = sc_sma
opt_par_2 = opt_par

# get RMSE of the model

if(opt_par_2[1] == "ema")
{
  RMSE_2 = sc_ema_2[a == opt_par_2[2],sqrt(MSE)]
} else{
  RMSE_2 = sc_sma_2[win_length == opt_par_2[2],sqrt(MSE)]
}


# get weights:

# heuristically:
lambda_1 = RMSE_1^(-2)/(RMSE_1^(-2)+RMSE_2^(-2))
lambda_2 = 1- lambda_1

mixed_fc = lambda_1*DT_1[,SST_hat]+lambda_2 * DT_2[,SST_hat]

DT_m = DT_1[,.(year,month,Lon,Lat,grid_id,YM,SST_bar)]
DT_m[,SST_hat := mixed_fc]

RMSE_mixed = DT_m[year %in% validation_years, sqrt(mean((SST_bar - SST_hat)^2,na.rm = TRUE))]

RMSEs_MM = data.table(method = "heuristically",lambda_1 = round(lambda_1,3),lambda_2 = round(lambda_2,3),RMSE_m = round(RMSE_mixed,4),RMSE_NCPM = round(RMSE_1,4),RMSE_GCFS1 = round(RMSE_2,4))


#tryout:

n = 100
lambda_vec = seq(0,1,length.out = 100)

RMSE_vec = c()

for(lambda_1 in lambda_vec)
{
  lambda_2 = 1-lambda_1
  
  # get a data table with mixed forecasts:
  mixed_fc = lambda_1*DT_1[,SST_hat]+lambda_2 * DT_2[,SST_hat]
  
  DT_m = DT_1[,.(year,month,Lon,Lat,grid_id,YM,SST_bar)]
  DT_m[,SST_hat := mixed_fc]
  
  RMSE_mixed = DT_m[year %in% validation_years, sqrt(mean((SST_bar - SST_hat)^2,na.rm = TRUE))]
  RMSE_vec = c(RMSE_vec,RMSE_mixed)
}

plot(lambda_vec,RMSE_vec)
abline(h=RMSE_1)

lambda_1 = lambda_vec[which(RMSE_vec == min(RMSE_vec))]
lambda_2 = 1-lambda_1
RMSE_mixed = RMSE_vec[which(RMSE_vec == min(RMSE_vec))]

RMSEs_MM = rbindlist(list(RMSEs_MM,
                          data.table(method = "tryout",lambda_1 = round(lambda_1,3),lambda_2 = round(lambda_2,3),RMSE_m = round(RMSE_mixed,4),RMSE_NCPM = round(RMSE_1,4),RMSE_GCFS1 = round(RMSE_2,4))))
            

#tryout out of sample for comparison

validation_years_2 = 2003:2013

n = 100
lambda_vec = seq(0,1,length.out = 100)

RMSE_vec = c()

for(lambda_1 in lambda_vec)
{
  lambda_2 = 1-lambda_1
  
  # get a data table with mixed forecasts:
  mixed_fc = lambda_1 * DT_1[,SST_hat] + lambda_2 * DT_2[,SST_hat]
  
  DT_m = DT_1[,.(year,month,Lon,Lat,grid_id,YM,SST_bar)]
  DT_m[,SST_hat := mixed_fc]
  
  RMSE_mixed = DT_m[year %in% validation_years_2, sqrt(mean((SST_bar - SST_hat)^2,na.rm = TRUE))]
  RMSE_vec = c(RMSE_vec,RMSE_mixed)
}


lambda_1 = lambda_vec[which(RMSE_vec == min(RMSE_vec))]
lambda_2 = 1-lambda_1

mixed_fc = lambda_1 * DT_1[,SST_hat] + lambda_2 * DT_2[,SST_hat]

DT_m = DT_1[,.(year,month,Lon,Lat,grid_id,YM,SST_bar)]
DT_m[,SST_hat := mixed_fc]

RMSE_mixed = DT_m[year %in% validation_years, sqrt(mean((SST_bar - SST_hat)^2,na.rm = TRUE))]

RMSEs_MM = rbindlist(list(RMSEs_MM,
                          data.table(method = "tryout_oos",lambda_1 = round(lambda_1,3),lambda_2 = round(lambda_2,3),RMSE_m = round(RMSE_mixed,4),RMSE_NCPM = round(RMSE_1,4),RMSE_GCFS1 = round(RMSE_2,4))))



# try to figure out th minimum more formally

MSE_1 = DT_1[year %in% validation_years, mean((SST_hat - SST_bar)^2,na.rm = TRUE)]
MSE_2 = DT_2[year %in% validation_years, mean((SST_hat - SST_bar)^2,na.rm = TRUE)]
MSE_cross = mean(DT_1[year %in% validation_years, SST_hat - SST_bar] * DT_2[year %in% validation_years, SST_hat - SST_bar],na.rm = TRUE)

f = function(lambda){return(lambda^2*MSE_1 + (1-lambda)^2*MSE_2 + 2*lambda*(1-lambda)*MSE_cross)}


lambda_opt = optimize(interval = c(0,1),f = function(lambda){return(lambda^2*MSE_1 + (1-lambda)^2*MSE_2 + 2*lambda*(1-lambda)*MSE_cross)})

lambda_1 = lambda_opt$minimum
lambda_2 = 1-lambda_1

mixed_fc = lambda_1*DT_1[,SST_hat]+lambda_2 * DT_2[,SST_hat]

DT_m = DT_1[,.(year,month,Lon,Lat,grid_id,YM,SST_bar)]
DT_m[,SST_hat := mixed_fc]

RMSE_mixed = DT_m[year %in% validation_years, sqrt(mean((SST_bar - SST_hat)^2,na.rm = TRUE))]

RMSEs_MM = rbindlist(list(RMSEs_MM,
                          data.table(method = "estimated",lambda_1 = round(lambda_1,3),lambda_2 = round(lambda_2,3),RMSE_m = round(RMSE_mixed,4),RMSE_NCPM = round(RMSE_1,4),RMSE_GCFS1 = round(RMSE_2,4))))



# out of sample estimation for comparison:


MSE_1 = DT_1[year %in% validation_years_2, mean((SST_hat - SST_bar)^2,na.rm = TRUE)]
MSE_2 = DT_2[year %in% validation_years_2, mean((SST_hat - SST_bar)^2,na.rm = TRUE)]
MSE_cross = mean(DT_1[year %in% validation_years_2, SST_hat - SST_bar] * DT_2[year %in% validation_years_2, SST_hat - SST_bar],na.rm = TRUE)

f = function(lambda){return(lambda^2*MSE_1 + (1-lambda)^2*MSE_2 + 2*lambda*(1-lambda)*MSE_cross)}

lambda_opt = optimize(interval = c(0,1),f = function(lambda){return(lambda^2*MSE_1 + (1-lambda)^2*MSE_2 + 2*lambda*(1-lambda)*MSE_cross)})

lambda_1 = lambda_opt$minimum
lambda_2 = 1-lambda_1

mixed_fc = lambda_1*DT_1[,SST_hat]+lambda_2 * DT_2[,SST_hat]

DT_m = DT_1[,.(year,month,Lon,Lat,grid_id,YM,SST_bar)]
DT_m[,SST_hat := mixed_fc]

RMSE_mixed = DT_m[year %in% validation_years, sqrt(mean((SST_bar - SST_hat)^2,na.rm = TRUE))]

RMSEs_MM = rbindlist(list(RMSEs_MM,
                          data.table(method = "estimated_oos",lambda_1 = round(lambda_1,3),lambda_2 = round(lambda_2,3),RMSE_m = round(RMSE_mixed,4),RMSE_NCPM = round(RMSE_1,4),RMSE_GCFS1 = round(RMSE_2,4))))


# return results

RMSEs_MM


########### Handpick your preferred lambdas and conduct the mixture #########

# for comparison first estimate lambda out of sample

lambda_1 = RMSEs_MM[method == "tryout",lambda_1]
lambda_2 = 1-lambda_1

mixed_fc = lambda_1*DT_1[,SST_hat]+lambda_2 * DT_2[,SST_hat]

DT = DT_1[,.(year,month,Lon,Lat,grid_id,YM,SST_bar)]
DT[,SST_hat := mixed_fc]

#save

save_dir = paste0(save_dir_1,"mm_oos/")
dir.create(save_dir,showWarnings = FALSE)

save.image(file = paste0(save_dir,"setup.RData"))

# now, finally, choose weights and get mixture model:

lambda_1 = RMSEs_MM[method == "tryout",lambda_1]
lambda_2 = 1-lambda_1


mixed_fc = lambda_1*DT_1[,SST_hat]+lambda_2 * DT_2[,SST_hat]

DT = DT_1[,.(year,month,Lon,Lat,grid_id,YM,SST_bar)]
DT[,SST_hat := mixed_fc]
DT[,SST_NCPM := DT_1[,SST_hat]]
DT[,SST_GCFS1 := DT_2[,SST_hat]]

#save

save_dir = paste0(save_dir_1,"mm/")
dir.create(save_dir,showWarnings = FALSE)

save.image(file = paste0(save_dir,"setup.RData"))



# now, after finding the mixed prediction we can model the variance of it as usual and thereafter get a CRPS-value.
# To do this, just run master script #03 (variance modelling) twice with save_dir_mm and save_dir_mm_oos as the corresponding save_dirs.


############################################
#### get all scores into one data table ####
############################################

name_abbrs = c("Aut_2018_small","Aut_2018_small/GCFS1","Aut_2018_small/mm","Aut_2018_small/mm_oos")
save_dirs =  paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbrs,"/")

methods = c("NorCPM","GCFS1","mixed, i.s.","mixed, o.o.s.")

scores_dt = data.table(method = methods)

scores_dt[method == "mixed, i.s." , RMSE := RMSEs_MM[method == "tryout",RMSE_m]]
scores_dt[method == "mixed, o.o.s." , RMSE := RMSEs_MM[method == "tryout_oos",RMSE_m]]


# get RMSEs

for(i in 1:2)
{
  load(file = paste0(save_dirs[i],"setup.RData"))
  
  load(file = paste0(save_dirs[i],"scores.bc.",opt_par[1],".RData"))
  
  if(opt_par[1] == "ema") sc = sc_ema
  if(opt_par[1] == "sma") sc = sc_sma
  
  scores_dt[method == methods[i], RMSE := sc[,sqrt(min(MSE))]]
}

# get CRPSs

for(i in 1:4)
{
  load(file = paste0(save_dirs[i],"setup.RData"))
  
  a = load(file = paste0(save_dirs[i],"scores.bc.sd.",opt_par_var[1],".Rdata"))
  
  if(opt_par_var[1] == "ema") sc = sc_ema_var
  if(opt_par_var[1] == "sma") sc = sc_sma_var
  
  scores_dt[method == methods[i], CRPS := sc[,min(CRPS)]]
}

# save

save(RMSEs_MM,scores_dt,file = paste0(save_dir_1,"scores.mixed.model.RData"))

