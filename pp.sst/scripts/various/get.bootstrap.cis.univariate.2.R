

# get permutation tests etc for the Full/lv run. To this end we run the oos script number 3 again, since in the original run SD_hat was not computed 
# with SMAs. The first script ran for a long time and then ran into an error. It si corrected now and we don't run everything again but try 
# kicking off where we started.

rm(list = ls())

time_s3oos = proc.time()

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "Full/lv" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))

# add this now, for the future

comment = 'This is the full globe, lv stands for long validation (years 2001-2016).'

###### getting sample variances of ensemble  ######

DT[,var_bar := (SST_hat - SST_bar)^2]



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

opt_ema_par_tr_sd = opt_par_sd[year == min(year),par_ema]

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

script_counter = 3


save.image(file = paste0(save_dir,'setup.RData'))

##############################################

### correct colnames ###

DT[,Bias_Est_EMA := Bias_Est]
DT[,SST_hat_EMA := SST_hat]

######################################################################################
### permutation tests and bootstrap confidence intervals for al model combinations ###
######################################################################################

models = c('lr_m','lr_loc','lr_both','sma','ema')

scores = c('MSE','CRPS')

perm_test_dt = as.data.table(expand.grid(model1 = models,model2 = models,score = scores))

################

N=5000 # number of permutations for permutation test and resamples for bootstrap

q_probs = c(0.025,0.975) # probabilities for bootstrap quantiles of the score differences

# MSE:

for(mod1 in models)
{
  print(mod1)
  
  dat1 = MSE_linear_models[,get(mod1)]
  
  for(mod2 in models)
  {
    dat2 = MSE_linear_models[,get(mod2)]
    
    p_val = permutation_test_difference(dat1,dat2,N)$p_val
    
    perm_test_dt[model1 == mod1 & model2 == mod2 & score == 'MSE',p_value := p_val]
    
    bt = bootstrap_difference(dat1,dat2,N = N,q_prob = q_probs)
    
    perm_test_dt[model1 == mod1 & model2 == mod2 & score == 'MSE',paste0('bt_q',q_probs[1]) := bt$q[1]]
    perm_test_dt[model1 == mod1 & model2 == mod2 & score == 'MSE',paste0('bt_q',q_probs[2]) := bt$q[2]]
  }
}

# CRPS (for variance estimation)


for(mod1 in models)
{
  print(mod1)
  
  dat1 = CRPS_comparison[,get(mod1)]
  
  for(mod2 in models)
  {
    dat2 = CRPS_comparison[,get(mod2)]
    
    p_val = permutation_test_difference(dat1,dat2,N)$p_val
    
    perm_test_dt[model1 == mod1 & model2 == mod2 & score == 'CRPS',p_value := p_val]
    
    bt = bootstrap_difference(dat1,dat2,N = N,q_prob = q_probs)
    
    perm_test_dt[model1 == mod1 & model2 == mod2 & score == 'CRPS',paste0('bt_q',q_probs[1]) := bt$q[1]]
    perm_test_dt[model1 == mod1 & model2 == mod2 & score == 'CRPS',paste0('bt_q',q_probs[2]) := bt$q[2]]
  }
}


###############33

save(perm_test_dt,file = paste0(save_dir,'perm_test_univ_mod.RData'))


save.image(file = paste0(save_dir,"setup.RData"))






