# time univariate methods: 

# bias correction:

##### setting up ######

rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "NAO" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))

#### set number of cores ###

mc_cores = 1

#start timer:

time_sma = proc.time()

###### run bias analysis for simple moving averages ######

num_years = DT[,range(year)][2] - DT[,range(year)][1] + 1


dummy_function = function(k){
  temp = bias_correct(dt = DT, method = "sma", par_1 = k)
  
  temp[,SST := DT[,SST_bar]]
  
  MSE = temp[, mean((SST_hat - SST)^2,na.rm = TRUE),by = .(year,month)]
  setnames(MSE,old = 3,new = 'MSE')
  
  MSE[,l:=k]
  return(MSE)
}

sc_sma = parallel::mclapply(X = 1:(num_years-1), FUN = dummy_function, mc.cores = mc_cores)
sc_sma = rbindlist(sc_sma)

############## get optimal parameter and estimate bias ################

### sma ###

temp = sc_sma[,mean(MSE),by = l]
setnames(temp,2,'MSE')
sma_par = temp[,which.min(MSE)]
sma_mse = temp[,min(MSE)]

temp_sma = bias_correct(dt = DT,
                        method = 'sma',
                        par_1 = sma_par)
DT[,Bias_Est_SMA := temp_sma[,Bias_est]][,SST_hat_sma := temp_sma[,SST_hat]]

### stop timer ###

#start timer:

time_sma = proc.time() - time_sma
print(time_sma)

##################################################################################

### exponential moving averages ###

#start timer:

time_ema = proc.time()


###### run bias analysis for simple moving averages ######

###### run bias analysis for exponential moving averages ######

par_vec = seq(0.05,0.4,length.out = 24) 


dummy_function = function(k){
  temp = bias_correct(dt = DT, method = "ema", par_1 = par_vec[k])
  
  temp[,SST := DT[,SST_bar]]
  
  MSE = temp[,mean((SST_hat - SST)^2,na.rm = TRUE),by = .(year,month)]
  setnames(MSE,old = 3,new = 'MSE')
  
  MSE[,a:=par_vec[k]]
  return(MSE)
}

sc_ema = parallel::mclapply(X = 1:length(par_vec), FUN = dummy_function,mc.cores = mc_cores)
sc_ema = rbindlist(sc_ema)


############## get optimal parameter and estimate bias ################

### ema ###

temp = sc_ema[,mean(MSE),by = a]
setnames(temp,2,'MSE')
ema_par = temp[,which.min(MSE)]
ema_mse = temp[,min(MSE)]

temp_ema = bias_correct(dt = DT,
                        method = 'ema',
                        par_1 = ema_par)
DT[,Bias_Est_EMA := temp_ema[,Bias_est]][,SST_hat_ema := temp_ema[,SST_hat]]


#stop timer:

time_ema = proc.time() - time_ema
print(time_ema)

############# bias estimation linear regression by both ###############

time_bb = proc.time()

if(mc_cores == 1)
{
  DT = bias_lr_bb(DT,validation_years = 2016)
} else{
  DT = bias_lr_bb_par(DT,months = months,validation_years = 2016,mc_cores = mc_cores)
}

time_bb = proc.time() - time_bb
print(time_bb)


############# bias estimation linear regression by month ###############

time_bm = proc.time()

if(mc_cores == 1)
{
  DT = bias_lr_bm(DT,validation_years = 2016)
} 

time_bm = proc.time() - time_bm
print(time_bm)



time_bl = proc.time()

if(mc_cores == 1)
{
  DT = bias_lr_bl(DT,validation_years = 2016)
} 

time_bl = proc.time() - time_bl
print(time_bl)
####################################################

############## variance estimation ###################


##### getting scores for simple moving averages ######

#start timer:

time_sma_var = proc.time()

###

num_years = DT[,range(year)][2] - DT[,range(year)][1] + 1

dummy_function = function(k){
  temp = sd_est_2(dt = DT, method = "sma", par_1 = k)
  
  temp[,SST_hat := DT[,SST_hat]]
  temp[,SST_bar := DT[,SST_bar]]
  
  CRPS = temp[, crps_na_rm(SST_bar,mean = SST_hat,sd = SD_hat),by = .(year,month)]
  setnames(CRPS,old = 3,new = 'CRPS')
  
  CRPS[,l:=k]
  return(CRPS)
}


sc_sma_var = parallel::mclapply(X = 1:(num_years-1), FUN = dummy_function, mc.cores = mc_cores)
sc_sma_var = rbindlist(sc_sma_var)



### sma ###

temp = sc_sma_var[,mean(CRPS,na.rm = TRUE),by = l]
setnames(temp,2,'CRPS')
sma_par_var = temp[,which.min(CRPS)]
sma_crps_var = temp[,min(CRPS)]

temp_sma_var = sd_est_2(dt = DT,
                        method = 'sma',
                        par_1 = sma_par_var)
DT[,SD_hat_SMA := temp_sma_var[,SD_hat]]

#stop timer:

time_sma_var = proc.time() - time_sma_var


########## exponential moving averages ##########

#start timer:

time_ema_var = proc.time() 


###### getting scores for exponential moving averages ######

par_vec = seq(0.01,0.4,length.out = 24)


dummy_function = function(k){
  temp = sd_est_2(dt = DT, method = "ema", par_1 = k)
  
  temp[,SST_hat := DT[,SST_hat]]
  temp[,SST_bar := DT[,SST_bar]]
  
  CRPS = temp[, crps_na_rm(SST_bar,mean = SST_hat,sd = SD_hat),by = .(year,month)]
  setnames(CRPS,old = 3,new = 'CRPS')
  
  CRPS[,l:=k]
  return(CRPS)
}

sc_ema_var = parallel::mclapply(X = 1:length(par_vec), FUN = dummy_function,mc.cores = mc_cores)
sc_ema_var = rbindlist(sc_ema_var)


temp = sc_ema_var[,mean(CRPS,na.rm = TRUE),by = l]
setnames(temp,2,'CRPS')
ema_par_var = temp[,which.min(CRPS)]
ema_crps_var = temp[,min(CRPS)]

temp_ema_var = sd_est_2(dt = DT,
                        method = 'ema',
                        par_1 = ema_par_var)
DT[,SD_hat_EMA := temp_ema_var[,SD_hat]]

#stop timer:

time_ema_var = proc.time() - time_ema_var

################### NGR ####################

#stop timer:

time_bb_var = proc.time()

DT = var_est_NGR_bb(DT,validation_years = 2016,mc.cores = mc_cores)

time_bb_var = proc.time() - time_bb_var

# by month:

time_bm_var = proc.time()

DT = var_est_NGR_bm(DT,validation_years = 2016)

time_bm_var = proc.time() - time_bm_var

# by location:

time_bl_var = proc.time()

DT = var_est_NGR_bl(DT,validation_years = 2016,mc.cores = mc_cores)

time_bl_var = proc.time() - time_bl_var



############

mod_names = c('sma','ema','bb','bm','bl')

var_names = paste0(mod_names,'_var')

names = c(mod_names,var_names)


times = data.table(model = names)

for(nn in names)
{
  times[model == nn, c('user','system','elapsed') := as.list(get(paste0('time_',nn))[1:3])]
  
}



save(times,file = paste0(save_dir,'times.univ.methods.RData'))
