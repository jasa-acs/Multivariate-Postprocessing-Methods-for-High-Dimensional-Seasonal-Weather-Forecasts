
##################################################################§§§############################

#############  side script 3.4 - month-by-month analysis for univariate methods  ################

#####################################################################§§§#########################

# This script analyses the univatiate postprocessing methods month by month, see supplementary material.


##### setting up ######

rm(list = ls())

options(max.print = 1e3)

library(pp.sst)
library(boot)
library(ggplot2)
library(latex2exp)

name_abbr = "Full" 

save_dir = file.path('~','SST','Derived', name_abbr)

load(file = paste0(save_dir,"setup.RData"))


time_s34 = proc.time()


load(file = paste0(save_dir,'univ_scores_comp_NGR.RData'))

load(file = paste0(save_dir,'bc_scores_comp_NGR.RData'))


#################################

#### month by month analysis ####

#################################

models_mse = c('lr_bm','lr_bl','lr_bb','lr_bm_la','sma','ema')

models_crps = c('lr_bm','lr_bl','lr_bb','sma','ema')

N=5000 # number of permutations for permutation test and resamples for bootstrap

q_probs = c(0.025,0.975) # probabilities for bootstrap quantiles of the score differences



perm_test_dt_mse_bm = as.data.table(expand.grid(model1 = models_mse, model2 = models_mse, month = 1:12))

################

for(m in 1:12)
{
  print(m)
  
  N=5000 # number of permutations for permutation test and resamples for bootstrap
  
  q_probs = c(0.025,0.975) # probabilities for bootstrap quantiles of the score differences
  
  # MSE:
  
  for(mod1 in models_mse)
  {
    print(mod1)
    
    dat1 = MSE_linear_models[month == m,get(paste0('MSE_',mod1))]
    
    for(mod2 in models_mse)
    {
      dat2 = MSE_linear_models[month == m,get(paste0('MSE_',mod2))]
      
      p_val = permutation_test_difference(dat1,dat2,N)$p_val
      
      perm_test_dt_mse_bm[model1 == mod1 & model2 == mod2 & month == m, p_value := p_val]
      
      bt = bootstrap_difference(dat1,dat2,N = N,q_prob = q_probs)
      
      perm_test_dt_mse_bm[model1 == mod1 & model2 == mod2 & month == m,paste0('bt_q',q_probs[1]) := bt$q[1]]
      perm_test_dt_mse_bm[model1 == mod1 & model2 == mod2 & month == m ,paste0('bt_q',q_probs[2]) := bt$q[2]]
    }
  }
  
}

perm_test_dt_mse_bm[,score:= 'MSE']

# CRPS (for variance estimation)

models_crps = c('lr_bm','lr_bl','lr_bb','sma','ema')

perm_test_dt_crps_bm = as.data.table(expand.grid(model1 = models_crps,model2 = models_crps,month = 1:12))

for(m in 1:12)
{
  print(m)
  for(mod1 in models_crps)
  {
    print(mod1)
    
    dat1 = CRPS_comparison[month == m,get(mod1)]
    
    for(mod2 in models_crps)
    {
      dat2 = CRPS_comparison[month == m,get(mod2)]
      
      p_val = permutation_test_difference(dat1,dat2,N)$p_val
      
      perm_test_dt_crps_bm[model1 == mod1 & model2 == mod2 & month == m, p_value := p_val]
      
      bt = bootstrap_difference(dat1,dat2,N = N,q_prob = q_probs)
      
      perm_test_dt_crps_bm[model1 == mod1 & model2 == mod2 & month == m, paste0('bt_q',q_probs[1]) := bt$q[1]]
      perm_test_dt_crps_bm[model1 == mod1 & model2 == mod2 & month == m, paste0('bt_q',q_probs[2]) := bt$q[2]]
    }
  }
}
perm_test_dt_crps_bm[,score:= 'CRPS']


perm_test_dt_bm= rbindlist(list(perm_test_dt_mse_bm,perm_test_dt_crps_bm))


###############################

#### get bootstrap samples ####

###############################


### get skill scores ###

for(mod in models_mse)
{
  MSE_linear_models[,paste0('SS_',mod,'_ema') := 1 - get(paste0('MSE_',mod))/MSE_ema]  
}


for(mod in models_crps)
{
  CRPS_comparison[,paste0('SS_',mod,'_ema') := 1 - get(mod)/ema]  
}

####

# get bootstrap samples for skill scores:

qq1 = 0.025
qq2 = 0.975

N = 500

boot_dt_mse = as.data.table(expand.grid(model = models_mse,month = 1:12))

for(m in 1:12)
{
  print(m)
  for(mod in models_mse)
  {
    data = MSE_linear_models[month == m, get(paste0('SS_',mod,'_ema'))]
    
    bt_fct = function(x,inds){return(mean(x[inds]))}
    
    bt = boot(data = data, bt_fct, R = N)
    
    boot_dt_mse[model == mod & month == m, SS := bt$t0]
    boot_dt_mse[model == mod & month == m, q1 := quantile(bt$t,qq1)]
    boot_dt_mse[model == mod & month == m, q2 := quantile(bt$t,qq2)]
  }
}


boot_dt_crps = as.data.table(expand.grid(model = models_crps,month = 1:12))

for(m in 1:12)
{
  print(m)
  for(mod in models_crps)
  {
    data = CRPS_comparison[month == m, get(paste0('SS_',mod,'_ema'))]
    
    bt_fct = function(x,inds){return(mean(x[inds]))}
    
    bt = boot(data = data, bt_fct, R = N)
    
    boot_dt_crps[model == mod & month == m, SS := bt$t0]
    boot_dt_crps[model == mod & month == m, q1 := quantile(bt$t,qq1)]
    boot_dt_crps[model == mod & month == m, q2 := quantile(bt$t,qq2)]
  }
}



##################
#### plotting ####

#####

title_size = 22
labels_size = 18
ticks_size = 14


mse_labs = list(expression('NGR'[s]),
                expression('NGR'['s,m']),
                expression('NGR'[m]^'la'),
                expression('sma'),
                expression('ema'))

crps_labs = list(expression('NGR'[s]),
                 expression('NGR'['s,m']),
                 expression('sma'),
                 expression('ema'))
#####

p_mse = ggplot(data = boot_dt_mse[!(model %in% c('lr_bm','lr_bb_la','lr_bl_la'))],aes(x = month,y = SS,color = model)) + geom_line() + geom_point()

# add errorbars

p_mse = p_mse + geom_errorbar(aes(ymin = q1,ymax = q2,color = model,width = 0.4),alpha = 0.5)

# add custom title

p_mse = p_mse + labs(title = 'Skill score: MSE by month',
                     y = 'skill score') +
  theme(plot.title = element_text(size = title_size),
        axis.title = element_text(size = labels_size),
        axis.text = element_text(size = ticks_size)) +
  scale_x_continuous(breaks = seq(2,12,2)) +
  scale_color_discrete(labels = mse_labs)


pdf(paste0(plot_dir,'MSE_by_month.pdf'))
print(p_mse)
dev.off()

##### CRPS #####

p_crps = ggplot(data = boot_dt_crps[model != 'lr_bm'],aes(x = month,y = SS,color = model)) + geom_line() + geom_point()

# add errorbars

p_crps = p_crps + geom_errorbar(aes(ymin = q1,ymax = q2,color = model,width = 0.4),alpha = 0.5)

# add custom title

p_crps = p_crps + labs(title = 'Skill score: CRPS by month',
                       y = 'skill score') +
  theme(plot.title = element_text(size = title_size),
        axis.title = element_text(size = labels_size),
        axis.text = element_text(size = ticks_size)) +
  scale_x_continuous(breaks = seq(2,12,2)) +
  scale_color_discrete(labels = crps_labs)


pdf(paste0(plot_dir,'CRPS_by_month.pdf'))
print(p_crps)
dev.off()

#### save stuff ####

time_s34 = proc.time() - time_s34

save.image(file = paste0(save_dir,"setup.RData"))

