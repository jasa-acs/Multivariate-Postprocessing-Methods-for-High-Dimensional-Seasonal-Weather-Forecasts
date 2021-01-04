

#tryout, do somethings get better when we take more principal components and move the distance up to 4000 km?


##################################################################################################

###################  master script part 4 - setting up multivariate models  ######################

##################################################################################################

# This script sets up the different methods of multivariate post-processing, by computing and saving covariance estimates etc.
# As part of this script you should choose whether to consider centered or even standardized variables.
#
# Files generated:
#   
# Data files: PCA/sam_cov_m+_y++.RData, SE/cov_est_SE_m+_y++.RData, 
#             where + labels the months considered and ++ the years in the validation period.
#
# Requires previous run of 03.master.bias.correct 
# with the same value of name_abbr as below.

##### setting up ######

rm(list = ls())

time_s4 = proc.time()

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "NAO/lv" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))

# reset name abbr

name_abbr = "NAO/lv/tryout" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")
dir.create(save_dir,showWarnings = FALSE)

plot_dir = paste0('./figures/',name_abbr,'/')
dir.create(save_dir,showWarnings = FALSE)




##### get weight matrix for tapering and svds of data matrices over the training period #####


Tap_dir = paste0(save_dir,"Tap/")
dir.create(Tap_dir,showWarnings = FALSE)

wm = weight_mat(DT,L = 5000) 

num_loc = DT[year == min(year) & month == min(month)][!(is.na(SST_bar) | is.na(SST_hat)) ,.N]


for(y in validation_years)
{
  
  print(y)
  
  svd_by_month = function(m)
  {
    print(paste0("month =",m))  
    
    train_years = DT[month == m][year < y & year > min(year),][,unique(year)]
    
    data_mat = matrix(DT[month == m][!(is.na(SST_bar) | is.na(SST_hat)) & year %in% train_years,
                                     SST_bar - SST_hat],
                      nrow = num_loc)
    
    sam_cov_mat = 1/length(train_years) * data_mat %*% t(data_mat) 
    
    tap_cov_mat = wm * sam_cov_mat
    
    return(svd(tap_cov_mat))
  }
  sin_val_dec = parallel::mclapply(X = months,FUN = svd_by_month,mc.cores = mc_cores)
  
  save(sin_val_dec,file = paste0(Tap_dir,'svd_y',y,'.RData'))
}




###################################################
###################### PCA  #######################
###################################################

PCA_dir = paste0(save_dir,"PCA/")
dir.create(PCA_dir,showWarnings = FALSE)


# how many PCs should we use?

nPCs = c()

for(m in months)
{
  #plot(sin_val_dec[[m]]$d, main = paste0('month = ',m))
  
  sum_vec = cumsum(sin_val_dec[[m]]$d)
  sum_tot = sum_vec[length(sum_vec)]
  
  nPCs = c(nPCs,which(sum_vec > 0.95*sum_tot)[1])
  
}


PCA_cov(DT,weight_mat = wm, 
        Y = 1985:2016,
        M = months,
        nPCs = nPCs,
        save_years = validation_years,
        save_dir = PCA_dir)


###################################################
################## geostationary ##################
###################################################

GS_dir = paste0(save_dir, "GS/")
dir.create(GS_dir, showWarnings = FALSE)

geostationary_training(dt = DT, 
                       training_years = training_years,
                       m = months,
                       save_dir = GS_dir,mc_cores)

########################################
################ ECC  ##################
########################################

ECC_dir = paste0(save_dir, "ECC/")
dir.create(ECC_dir, showWarnings = FALSE)

#####

####### generate forecasts ###############

# specifications for the desired forecasts:

fc_years = validation_years
fc_months = months
fc_ens_size = 500

mod_vec = c('PCA_mc','PCA_ac','GS','ECC')


###################################################
###################### PCA  #######################
###################################################

PCA_fc_mc = forecast_PCA_mult_corr(DT, 
                                   Y = fc_years,
                                   M = fc_months,
                                   n = fc_ens_size,
                                   nPCs = nPCs,
                                   cov_dir = PCA_dir)

save(PCA_fc_mc,file = paste0(PCA_dir,"fc_mc.RData"))

rm(PCA_fc_mc)
gc()

PCA_fc_ac = forecast_PCA_add_corr(DT, 
                                  Y = fc_years,
                                  M = fc_months,
                                  n = fc_ens_size,
                                  nPCs = nPCs,
                                  cov_dir = PCA_dir)

save(PCA_fc_ac,file = paste0(PCA_dir,"fc_ac.RData"))

rm(PCA_fc_ac)

gc()

###################################################
################## geostationary ##################
###################################################

GS_fc = forecast_GS(DT,
                    Y = validation_years,
                    M = months,
                    n = fc_ens_size,
                    var_dir = GS_dir,
                    mc_cores = mc_cores)

save(GS_fc,file = paste0(GS_dir,"fc.RData"))
rm(GS_fc)

gc()

########################################
################ ECC  ##################
########################################

ECC_fc = forecast_ECC(DT,
                      Y = validation_years,
                      M = months,
                      ens_size = ens_size)

save(ECC_fc,file = paste0(ECC_dir,"fc.RData"))
rm(ECC_fc)

save.image(file = paste0(save_dir,"setup.RData"))


