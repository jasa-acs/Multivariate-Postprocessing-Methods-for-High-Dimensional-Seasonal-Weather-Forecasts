
#######################################################################################

###################  master script part 3 - variance estimation  ######################

#######################################################################################

# This script estimates the marginal variance of the ensemble forecast by 
# simple and exponential moving averages and computes average scores for each method.
# It then selects and applies the optimal way of variance estimation and
# complements the data table DT by a new column SD_hat.
#
# 
##### setting up ######

rm(list = ls())

options(max.print = 1e3)

library(pp.sst)

name_abbr = "pp.sst/Full" 

save_dir = file.path('~','SST','Derived', name_abbr)

load(file = paste0(save_dir,"setup.RData"))

# start timer:

time_s3 = proc.time()

###### getting sample variances of ensemble  ######

DT[,var_bar := (SST_hat - SST_bar)^2]


#####################################################################
###### running variance estimation with simple moving averages ######


past_0 = 5  # How many years before validation period do we consider for the CRPS:
            # For prediction in year y, the parameters are chosen to minimize the MSE over the years y_0 - past_0, ... ,y - 1,
            # where y_0 is the first year of the validation period.

num_years = max(validation_years) - DT[,min(year)]

# get SDs for a range of parameters and all relevant years

SD_by_par = function(k)
  {
  sd_hat = sd_est_2(DT,method = 'sma',par_1 = k)
  return(sd_hat[year %between% c(min(validation_years - past_0),max(validation_years)),SD_hat])
  }
  
win_length = 1 : (num_years-1)

SDs = parallel::mclapply(X = win_length, FUN = SD_by_par,mc.cores = mc_cores)


# restrict to validation years + the past_0 years before 

SD_dt = DT[year %between% c(min(validation_years - past_0),max(validation_years)),
                 .(year,month,Lon,Lat,SST_bar,SST_hat)]

for(k in win_length)
{print(k)
  SD_dt[,paste0('l',k):= SDs[[k]]]
  SD_dt[,paste0('crps',k):= crps_na_rm(y = SST_bar,mean = SST_hat,sd = eval(parse(text = paste0('l',k))))]
}




#get mean CRPS by month, for each year averaged over all previous years contained in SD_dt

sc_sd_sma = NULL

for(y in validation_years)
{print(y)
  mean_sc_sd_bm = function(m)
  {
    temp = data.table(year = y ,month = m,SD_dt[year < y & month == m,
                                                      lapply(X = .SD,FUN = mean,na.rm = TRUE),
                                                      .SDcols = paste0('crps',win_length)])
    return(temp)     
  }
  
  CRPS_y = rbindlist(parallel::mclapply(X = months,FUN = mean_sc_sd_bm,mc.cores = mc_cores))
  sc_sd_sma = rbindlist(list(sc_sd_sma,CRPS_y))
}


save(sc_sd_sma, file = paste0(save_dir,"scores.bc.sd.sma.RData"))


####################################################################
###### run variance estimation by exponential moving averages ######

par_vec = seq(0.05, 0.4, length.out = 24) 

SD_by_par = function(a)
{
  sd_hat = sd_est_2(dt = DT,
                  method = "ema", 
                  par_1 = a)
  return(sd_hat[year %between% c(min(validation_years - past_0),max(validation_years)),SD_hat])
}


SDs = parallel::mclapply(X = par_vec, FUN = SD_by_par,mc.cores = mc_cores)

# restrict to validation years + the past_0 years before 

SD_dt = DT[year %between% c(min(validation_years - past_0),max(validation_years)),
                 .(year,month,Lon,Lat,SST_bar,SST_hat)]


for(a in par_vec)
{ind = which(par_vec == a )
print(ind)
ra = round(a,4)
SD_dt[,paste0('a',ra):= SDs[[ind]]]
SD_dt[,paste0('crps',ra):= crps_na_rm(y = SST_bar,mean = SST_hat,sd = eval(parse(text = paste0('a',ra))))]
}

#get mean crps by month, averaged over all previous years contained in SD_dt

sc_sd_ema = NULL

for(y in validation_years)
{print(y)
  mean_sc_sd_bm = function(m)
  {
    temp = data.table(year = y ,month = m,SD_dt[year < y & month == m,
                                                lapply(X = .SD,FUN = mean,na.rm = TRUE),
                                                .SDcols = paste0('crps',round(par_vec,4))])
    return(temp)     
  }
  
  CRPS_y = rbindlist(parallel::mclapply(X = months,FUN = mean_sc_sd_bm,mc.cores = mc_cores))
  sc_sd_ema = rbindlist(list(sc_sd_ema,CRPS_y))
}


save(sc_sd_ema, file = paste0(save_dir,"scores.bc.sd.ema.RData"))



###### getting mean scores for different ways of variance estimation ######

load(paste0(save_dir,"scores.bc.sd.sma.RData"))
load(paste0(save_dir,"scores.bc.sd.ema.RData"))

# get means over all months for sma

msc_sd_sma = sc_sd_sma[,lapply(.SD,mean),by = year][,month := NULL]

# get minimum for each row

row_min = NULL
min_l = NULL

for(i in 1:msc_sd_sma[,.N])
{
  row = msc_sd_sma[,-1,with = FALSE][i,]
  row_min = c(row_min,min(row))
  min_l = c(min_l,which.min(row))
}

msc_sd_sma[,min_crps := row_min][,min_l := min_l]


# get means over all months for ema

msc_sd_ema = sc_sd_ema[,lapply(.SD,mean),by = year][,month := NULL]

# get minimum for each row

row_min = NULL
min_a = NULL

for(i in 1:msc_sd_sma[,.N])
{
  row = msc_sd_ema[,-1,with = FALSE][i,]
  row_min = c(row_min,min(row))
  min_a = c(min_a,par_vec[which.min(row)])
}

msc_sd_ema[,min_crps := row_min][,min_a := min_a]


###### finding optimal way of variance estimation for each year in the validation period ######

# exponential moving averages are preferred, as the parameter estimation for them typically is more stable ######

opt_par_sd = data.table(year = validation_years,method = NA_character_,par = NA_real_)


opt_par_sd = data.table(year = validation_years)

opt_par_sd[,crps_sma := msc_sd_sma[,min_crps]]
opt_par_sd[,par_sma := msc_sd_sma[,min_l]]

opt_par_sd[,crps_ema := msc_sd_ema[,min_crps]]
opt_par_sd[,par_ema := msc_sd_ema[,min_a]]


# exponential moving averages are more stable, so we give them a slight edge:
opt_par_sd[, par := par_sma * (crps_sma < 0.95 * crps_ema) + par_ema * (crps_sma >= 0.95 * crps_ema)]


opt_par_sd[crps_sma < 0.95 * crps_ema, method :='sma' ]
opt_par_sd[crps_sma > 0.95 * crps_ema, method :='ema' ]


###############################################################
##### variance correction year by year, validation period #####

# sma

for(y in validation_years)
{
  print(y) 
  temp = sd_est_2(dt = DT,
                        method = 'sma',
                        par_1 = opt_par_sd[year == y,par_sma])[year == y,]
  DT[year == y,SD_hat_sma := temp[,SD_hat]]
}

rm(temp)

# ema

for(y in validation_years)
{
  print(y) 
  temp = sd_est_2(dt = DT,
                  method = 'ema',
                  par_1 = opt_par_sd[year == y,par_ema])[year == y,]
  DT[year == y,SD_hat_ema := temp[,SD_hat]]
}

rm(temp)



# For the training years, the variance correction considers also the future (contained in the training period),
# and uses the parameter estimated for the first validation year:

#sma

opt_sma_par_tr_sd = opt_par_sd[year == min(year),par_sma]

DT_sma_tr_sd = sd_est_2(dt = DT[year %in% training_years],
                        method = 'sma', par_1 = opt_sma_par_tr_sd,
                        twosided = TRUE)

DT[year %in% training_years, SD_hat_sma := DT_sma_tr_sd[,.(SD_hat)]]

rm(DT_sma_tr_sd,opt_sma_par_tr_sd)

#ema

opt_ema_par_tr_sd = opt_par[year == min(year),par_ema]

DT_ema_tr_sd = sd_est_2(dt = DT[year %in% training_years],
                         method = 'ema',par_1 = opt_ema_par_tr_sd,
                         twosided = TRUE)

DT[year %in% training_years, SD_hat_ema := DT_ema_tr_sd[,SD_hat]]

rm(DT_ema_tr_sd,opt_ema_par_tr_sd)


######################################################
### pick optimal performance method for each year: ###

ema_years = opt_par_sd[method == 'ema',year]

if(min(validation_years) %in% ema_years) ema_years = c(training_years,ema_years)

DT[year %in% ema_years, SD_hat := SD_hat_ema]

DT[!(year %in% ema_years), SD_hat := SD_hat_sma]





#### time, update script counter, save ####

time_s3 = proc.time() - time_s3

script_counter = 3

save.image(file = paste0(save_dir,"setup.RData"))
