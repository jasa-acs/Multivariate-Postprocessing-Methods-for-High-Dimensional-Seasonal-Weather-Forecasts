
##################################################################§§§######################

#############  side script 3.3 - permutation tests for univariate methods  ################

#####################################################################§§§###################

# This script conducts permutation tests in order to assess significance of the score differences for the univariate methods (both MSE and CRPS).


##### setting up ######

rm(list = ls())

options(max.print = 1e3)

library(pp.sst)

name_abbr = "Full" 

save_dir = file.path('~','SST','Derived', name_abbr)

load(file = paste0(save_dir,"setup.RData"))

#############

time_s33 = proc.time()

load(file = paste0(save_dir,'univ_scores_comp_NGR.RData'))

load(file = paste0(save_dir,'bc_scores_comp_NGR.RData'))


######################################################################################
### permutation tests and bootstrap confidence intervals for al model combinations ###
######################################################################################

models_mse = c('lr_bm','lr_bl','lr_bb','lr_bm_la','sma','ema')

perm_test_dt_mse = as.data.table(expand.grid(model1 = models_mse,model2 = models_mse))

################

N=5000 # number of permutations for permutation test and resamples for bootstrap

q_probs = c(0.025,0.975) # probabilities for bootstrap quantiles of the score differences

# MSE:

for(mod1 in models_mse)
{
  print(mod1)
  
  dat1 = MSE_linear_models[,get(paste0('MSE_',mod1))]
  
  for(mod2 in models_mse)
  {
    dat2 = MSE_linear_models[,get(paste0('MSE_',mod2))]
    
    p_val = permutation_test_difference(dat1,dat2,N)$p_val
    
    perm_test_dt_mse[model1 == mod1 & model2 == mod2,p_value := p_val]
    
    bt = bootstrap_difference(dat1,dat2,N = N,q_prob = q_probs)
    
    perm_test_dt_mse[model1 == mod1 & model2 == mod2 ,paste0('bt_q',q_probs[1]) := bt$q[1]]
    perm_test_dt_mse[model1 == mod1 & model2 == mod2 ,paste0('bt_q',q_probs[2]) := bt$q[2]]
  }
}

perm_test_dt_mse[,score:= 'MSE']

# CRPS (for variance estimation)

models_crps = c('lr_bm','lr_bl','lr_bb','sma','ema')

perm_test_dt_crps = as.data.table(expand.grid(model1 = models_crps,model2 = models_crps))

for(mod1 in models_crps)
{
  print(mod1)
  
  dat1 = CRPS_comparison[,get(mod1)]
  
  for(mod2 in models_crps)
  {
    dat2 = CRPS_comparison[,get(mod2)]
    
    p_val = permutation_test_difference(dat1,dat2,N)$p_val
    
    perm_test_dt_crps[model1 == mod1 & model2 == mod2 ,p_value := p_val]
    
    bt = bootstrap_difference(dat1,dat2,N = N,q_prob = q_probs)
    
    perm_test_dt_crps[model1 == mod1 & model2 == mod2 ,paste0('bt_q',q_probs[1]) := bt$q[1]]
    perm_test_dt_crps[model1 == mod1 & model2 == mod2 ,paste0('bt_q',q_probs[2]) := bt$q[2]]
  }
}

perm_test_dt_crps[,score:= 'CRPS']


perm_test_dt = rbindlist(list(perm_test_dt_mse,perm_test_dt_crps))

#### save stuff ####

time_s33 = proc.time() - time_s33

save.image(file = paste0(save_dir,"setup.RData"))

