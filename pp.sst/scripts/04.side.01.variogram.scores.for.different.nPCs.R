##############################################################################################################################

###################  side script 4.1 - compute variogram scores for a range of different numbers of PCs ######################

##############################################################################################################################

# This script computes variogram scores for a range of different number of PCs. This is often not conclusive, takes a lot of time,
# and is not recommended (see discussion section of the paper) 



##### setting up ######

rm(list = ls())


setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)
library(ggplot2)

name_abbr = "test" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))

# start timer:

time_s41 = proc.time()

# parameters

mc_cores = 3

nPCs = seq(5,30,by = 5) # number of PCs for which to compute the variogram score

mms = 1:6 # months (restrict this to make computation time bearable)

N_sim = 500 # number of simulations/size of predictive ensemble for variogram evaluation 

###########################################################################
# get for a range of nPCs the setup to simulate with n principle components

valPC_dir = paste0(PCA_dir,'valnPC/')
dir.create(valPC_dir,showWarnings = FALSE)


for(mm in mms)
{

  covs_mult_nPCs = function(nPC)
  {
  PCA_cov(DT,weight_mat = wm,
          Y = c(training_years,validation_years),
          M = mm,
          nPCs = nPC,
          save_years = validation_years,
          save_dir = valPC_dir,
          file_name = paste0('PCA_cov_n',nPC))
  
  }
  
  parallel::mclapply(nPCs,covs_mult_nPCs,mc.cores = mc_cores)

}

###########################################################################

dt_vs = data.table()

for(nnPC in nPCs)
{
  print(paste0('n = ',nnPC))
  # get forecast: (parallelized along months)
  
  print('forecasting...')
  fcs = forecast_PCA_mult_corr(DT, Y = validation_years, M = mms,
                               n = N_sim, nPCs = rep(nnPC,max(mms)),mc_cores = mc_cores,
                               cov_dir = valPC_dir,
                               cov_file_name = paste0('PCA_cov_n',nnPC))
  
  # get variogram score: (parallelized along months)
  
  print('computing variogram score...')
  dt_temp = var_sc_par(dt_fc = fcs,
                       years = validation_years,
                       ms = mms,
                       n = N_sim,
                       mc_cores = mc_cores)
  
  
  dt_vs = rbindlist(list(dt_vs,dt_temp[,nPC:= nnPC]))

}



#############################################################################

# get mean variogram score

dt_m_vs = dt_vs[,mean(vs),by = .(nPC,month)]

setnames(dt_m_vs,3,'vs')

p = ggplot() + geom_line(data = dt_m_vs,aes(x = nPC,y = vs,color = as.factor(month)))

print(p)

save(dt_vs,file = paste0(save_dir,'var_sc_different_nPCs.RData'))
