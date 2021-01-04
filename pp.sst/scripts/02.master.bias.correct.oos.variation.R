################################################################################

###################  master script part 2 - bias correct  ######################

################################################################################

# This script conducts the bias correction as described in the paper, by simple and moving averages.
# It complements the data table DT by two new columns, Bias_Est and SST_hat (which is just Ens_bar + Bias_Est, truncated
# at freezing temperature).



##### setting up ######

rm(list = ls())

options(max.print = 1e3)

library(pp.sst)

name_abbr = "pp.sst/Full" 

save_dir = file.path('~','SST','Derived', name_abbr)

load(file = paste0(save_dir,"setup.RData"))

# start timer:

time_s2 = proc.time()


###### run bias analysis for simple moving averages ######

past_0 = 5 # How many years before validation period do we consider for the MSEs:
           # For prediction in year y, the parameters are chosen to minimize the MSE over the years y_0 - past_0, ... ,y - 1,
           # where y_0 is the first year of the validation period.


num_years = max(validation_years) - DT[,min(year)]

win_length = 1 : (num_years-1)

MSE_by_par = function(k)
{
  temp = bias_correct(dt = DT, 
                        method = "sma", 
                        par_1 = k)
  return(temp[year %between% c(min(validation_years - past_0),max(validation_years)),SST_hat])
}


BC = parallel::mclapply(X = win_length, FUN = MSE_by_par,mc.cores = mc_cores)

# restrict to validation years + the past_0 years before 

Bias_est_dt = DT[year %between% c(min(validation_years - past_0),max(validation_years)),
                 .(year,month,Lon,Lat,SST_bar)]

for(k in win_length)
{print(k)
  Bias_est_dt[,paste0('l',k):= BC[[k]]]
  Bias_est_dt[,paste0('err',k):= (SST_bar - eval(parse(text = paste0('l',k))))^2]
}

# get MSE by year and month, for each year based on all previous years contained in Bias_est_dt

sc_sma = NULL

for(y in validation_years)
{print(y)
  mean_sc_bm = function(m)
  {
      temp = data.table(year = y ,month = m,Bias_est_dt[year < y & month == m,
                                                        lapply(X = .SD,FUN = mean,na.rm = TRUE),
                                                        .SDcols = paste0('err',win_length)])
      return(temp)     
  }
  
  MSE_y = rbindlist(parallel::mclapply(X = months,FUN = mean_sc_bm,mc.cores = mc_cores))
  sc_sma = rbindlist(list(sc_sma,MSE_y))
}


save(sc_sma, file = paste0(save_dir,"scores.bc.sma.RData"))


###### run bias analysis for exponential moving averages ######

par_vec = seq(0.05,0.4,length.out = 24) 

MSE_by_par = function(a)
{
  temp = bias_correct(dt = DT, 
                        method = "ema", 
                        par_1 = a)
  return(temp[year %between% c(min(validation_years - past_0),max(validation_years)),SST_hat])
}


BC = parallel::mclapply(X = par_vec, FUN = MSE_by_par,mc.cores = mc_cores)

# restrict to validation years + the past_0 years before 

Bias_est_dt = DT[year %between% c(min(validation_years - past_0),max(validation_years)),
                 .(year,month,Lon,Lat,SST_bar)]

for(a in par_vec)
{ind = which(par_vec == a )
  print(ind)
  ra = round(a,4)
  Bias_est_dt[,paste0('a',ra):= BC[[ind]]]
  Bias_est_dt[,paste0('err',ra):= (SST_bar - eval(parse(text = paste0('a',ra))))^2]
}

#get MSE by year and month, for each year based on all previous years contained in Bias_est_dt

sc_ema = NULL

for(y in validation_years)
{print(y)
  mean_sc_bm = function(m)
  {
    temp = data.table(year = y ,month = m,Bias_est_dt[year < y & month == m,
                                                      lapply(X = .SD,FUN = mean,na.rm = TRUE),
                                                      .SDcols = paste0('err',round(par_vec,4))])
    return(temp)     
  }
  
  MSE_y = rbindlist(parallel::mclapply(X = months,FUN = mean_sc_bm,mc.cores = mc_cores))
  sc_ema = rbindlist(list(sc_ema,MSE_y))
}


save(sc_ema, file = paste0(save_dir,"scores.bc.ema.RData"))


###### getting mean scores for different ways of bias correction ######

load(paste0(save_dir,"scores.bc.sma.RData"))
load(paste0(save_dir,"scores.bc.ema.RData"))

# get means over all months for sma

msc_sma = sc_sma[,lapply(.SD,mean),by = year][,month := NULL]

# get minimum for each row

row_min = NULL
min_l = NULL

for(i in 1:msc_sma[,.N])
{
  row = msc_sma[,-1,with = FALSE][i,]
  row_min = c(row_min,min(row))
  min_l = c(min_l,which.min(row))
}
  
msc_sma[,min_MSE := row_min][,min_l := min_l]


# get means over all months for ema

msc_ema = sc_ema[,lapply(.SD,mean),by = year][,month := NULL]

# get minimum for each row

row_min = NULL
min_a = NULL

for(i in 1:msc_sma[,.N])
{
  row = msc_ema[,-1,with = FALSE][i,]
  row_min = c(row_min,min(row))
  min_a = c(min_a,par_vec[which.min(row)])
}

msc_ema[,min_MSE := row_min][,min_a := min_a]


###### finding optimal way of bias correction for each year in the validation period ######

# exponential moving averages are preferred, as the parameter estimation for them typically is more stable ######

opt_par = data.table(year = validation_years)

opt_par[,mse_sma := msc_sma[,min_MSE]]
opt_par[,par_sma := msc_sma[,min_l]]

opt_par[,mse_ema := msc_ema[,min_MSE]]
opt_par[,par_ema := msc_ema[,min_a]]

# exponential moving averages are more stable, so we give them a slight edge:
opt_par[, par := par_sma * (mse_sma < 0.95 * mse_ema) + par_ema * (mse_sma >= 0.95 * mse_ema)]


opt_par[mse_sma < 0.95 * mse_ema, method :='sma' ]
opt_par[mse_sma > 0.95 * mse_ema, method :='ema' ]

#################################################################
######## bias correction year by year, validation period ########

# sma

for(y in validation_years)
{
  print(y) 
  temp_sma = bias_correct(dt = DT,
                            method = 'sma',
                            par_1 = opt_par[year == y,par_sma])[year == y,]
  DT[year == y,Bias_Est_SMA := temp_sma[,Bias_est]][year == y ,SST_hat_sma := temp_sma[,SST_hat]]
}
rm(temp_sma)

#ema

for(y in validation_years)
{
  print(y) 
  temp_ema = bias_correct(dt = DT,
                            method = 'ema',
                            par_1 = opt_par[year == y,par_ema])[year == y,]
  DT[year == y,Bias_Est_EMA := temp_ema[,Bias_est]][year == y ,SST_hat_ema := temp_ema[,SST_hat]]
}
rm(temp_ema)


# For the training years, the bias correction considers also the future (contained in the training period),
# and uses the parameter estimated for the first validation year:

#sma

opt_sma_par_tr = opt_par[year == min(year),par_sma]

DT_sma_tr = bias_correct(dt = DT[year %in% training_years],
                           method = 'sma',par_1 = opt_sma_par_tr,
                           twosided = TRUE)

DT[year %in% training_years, c('Bias_Est_SMA','SST_hat_sma') := DT_sma_tr[,.(Bias_est,SST_hat)]]

rm(DT_sma_tr,opt_sma_par_tr)

#ema

opt_ema_par_tr = opt_par[year == min(year),par_ema]

DT_ema_tr = bias_correct(dt = DT[year %in% training_years],
                           method = 'ema',par_1 = opt_ema_par_tr,
                           twosided = TRUE)

DT[year %in% training_years, c('Bias_Est_EMA','SST_hat_ema') := DT_ema_tr[,.(Bias_est,SST_hat)]]

rm(DT_ema_tr,opt_ema_par_tr)

######################################################
### pick optimal performance method for each year: ###

ema_years = opt_par[method == 'ema',year]

if(min(validation_years) %in% ema_years) ema_years = c(training_years,ema_years)

DT[year %in% ema_years, Bias_Est := Bias_Est_EMA]
DT[year %in% ema_years, SST_hat := SST_hat_ema]

DT[!(year %in% ema_years), Bias_Est := Bias_Est_SMA]
DT[!(year %in% ema_years), SST_hat := SST_hat_sma]




#### time, update script counter, save ####

time_s2 = proc.time() - time_s2

script_counter = 2

save.image(file = paste0(save_dir,"setup.RData"))
