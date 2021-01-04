##### setting up ######

rm(list = ls())

time_s4 = proc.time()

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "test" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))

#########################################

library(scoringRules)

# functions:

mean_crps_sim = function(l1,l2,DT_sm)
{
  DT_temp = copy(DT_sm)
  DT_temp[,"Bias_est_2" := sim_mov_av(l = l1, vec = Ens_bar - SST_bar, years = year,), by = .(Lon,Lat, month)]
  DT_temp[, SST_hat_2 := Ens_bar - Bias_est_2]
  
  
  DT_temp[,"SD_hat_2" := sqrt(sim_mov_av(l = l2, vec = (SST_hat_2 - SST_bar)^2, years = year,)), by = .(Lon,Lat, month)]
  
  DT_temp_val = DT_temp[year %in% validation_years]
  
  CRPS = crps_norm(y = DT_temp_val[, SST_bar],mean = DT_temp_val[,SST_hat_2],sd = DT_temp_val[,SD_hat_2])
  return(mean(CRPS))
}

#######

# maximize separately 

mean_MSE_sma = function(l1,DT_sm)
{
  DT_temp = copy(DT_sm)
  DT_temp[,"Bias_est_2" := sim_mov_av(l = l1, vec = Ens_bar - SST_bar, years = year,), by = .(Lon,Lat, month)]
  DT_temp_val = DT_temp[year %in%validation_years,]
  MSE = DT_temp_val[, mean((SST_bar-(Ens_bar - Bias_est_2))^2)]
  return(MSE)
}

########################################################




test_sep_vs_sim_est = function(DT,n = 3, N = 10,mc_cores = 8)
{

  ret_val = data.table()
  
  for(i in 1:N)
  {
    print(i)
  
    # get n water grid points and a random month:
    
    gids = sample(x = DT[!is.na(SST_hat),unique(grid_id)],n)
    m = sample(1:12,1)
    
    DT_sm = DT[grid_id %in% gids & month == m]
    
    
    # get SST_hat and SD_hat simultaneously
    
    # simple moving averages
    
    
    
    
    ####### get mean crps for all combinations: ########
    
    ls = 1:31
    
    
    
    
    crps_by_l = function(ll1)
    {
      ret_dt = data.table(l2 = ls)
      for(ll2 in ls)
      {
        ret_dt[ ll2 == l2, mCRPS:= mean_crps_sim(ll1,ll2,DT_sm)]
      }
    
      return(ret_dt[,l1:=ll1])
    }
    
    mean_crps_dt = rbindlist(parallel::mclapply(X = ls,FUN = crps_by_l,mc.cores = mc_cores))
    
    
    l1_min =  mean_crps_dt[which.min(mCRPS),l1]
    DT_sm[, Bias_est_2 := sim_mov_av(l = l1_min, vec = Ens_bar - SST_bar, years = year,), by = .(Lon,Lat, month)]
    DT_sm[, SST_hat_2 := Ens_bar - Bias_est_2]
    
    dt_sim_min = mean_crps_dt[which.min(mCRPS)][,MSE := DT_sm[year %in% validation_years,mean((SST_hat_2 - SST_bar)^2)]]
    
    dt_sim_min[,method := 'sim']
    
    
    
    ####### get mean crps for best MSE: ########
    
    ls = 1:31
    
    mean_mse_dt = data.table(l1 = ls)
    
    for(ll1 in ls)
    {
      mean_mse_dt[l1 == ll1, mse:= mean_MSE_sma(ll1,DT_sm)]
    }
    
    ll1 = mean_mse_dt[which.min(mse),l1]
    
    DT_sm[, Bias_est_3 := sim_mov_av(l = ll1, vec = Ens_bar - SST_bar, years = year,), by = .(Lon,Lat, month)]
    DT_sm[, SST_hat_3 := Ens_bar - Bias_est_3]
    
    mean_crps_dt_sep = data.table(l2 = ls)
    
    for(ll2 in ls)
    {
      mean_crps_dt_sep[ ll2 == l2, mCRPS:= mean_crps_sim(ll1,ll2,DT_sm)]
    }
    
    
    dt_sep_min = mean_crps_dt_sep[which.min(mCRPS)][,l1 := ll1]
    dt_sep_min[,MSE := DT_sm[year %in% validation_years,mean((SST_hat_3 - SST_bar)^2)]]
    dt_sep_min[,method := 'sep']
    
    
    ind_ret = rbindlist(list(dt_sim_min,dt_sep_min))[,ind:=i]
    
    ret_val = rbindlist(list(ret_val,ind_ret))
    
  }

  return(ret_val)
}

####

min_crps_diff = test_sep_vs_sim_est(DT, n = 5, N = 500 )

save(min_crps_diff,file = paste0(save_dir,'crps_diff.RData'))
