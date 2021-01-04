
##################################################################################################

###################  master script part 4 - setting up multivariate models  ######################

##################################################################################################

# This script sets up the different methods of multivariate post-processing, by computing and saving covariance estimates etc.
# As part of this script you should choose whether to consider centered or even standardized variables (in the paper this is not done, 
# it does not affect the results).


#### setting up ######

rm(list = ls())

options(max.print = 1e3)

library(pp.sst)

name_abbr = "pp.sst/Full" 

save_dir = file.path('~','SST','Derived', name_abbr)

load(file = paste0(save_dir,"setup.RData"))


time_s4 = proc.time()

##### centered variables? #####

# decide whether you want to work with SST, or with SST centered around climatology, or with SST standardized w.r.t. climatology 

SST = ""  # takes 'centered','standardized', or ''

clim_years = training_years


if(SST == "centered")
{
  DT = dt_transform_center(DT,clim_years)

  name_abbr = paste0(name_abbr,"/centered" )
  
  save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")
  dir.create(save_dir,showWarnings = FALSE)
  
  plot_dir = paste0("./figures/", name_abbr,"/")
  dir.create(plot_dir, showWarnings = FALSE)
  
}
if(SST == "standardized")
{
  DT = dt_transform_stan(DT,clim_years)

  name_abbr = paste0(name_abbr,"/standardized")
  
  save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")
  dir.create(save_dir,showWarnings = FALSE)
  
  plot_dir = paste0("./figures/", name_abbr,"/")
  dir.create(plot_dir, showWarnings = FALSE)
  
}


##### get weight matrix for tapering and svds of data matrices over the training period #####


Tap_dir = paste0(save_dir,"Tap/")
dir.create(Tap_dir,showWarnings = FALSE)

wm = weight_mat(DT,L = 2500) 

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

percentage = 90 # percentage of variance recovered, see discussion section in the paper

nPCs = c()

for(m in months)
{
  #plot(sin_val_dec[[m]]$d, main = paste0('month = ',m))
  load(paste0(Tap_dir,'svd_y',2010,'.RData'))
  
  sum_vec = cumsum(sin_val_dec[[m]]$d)
  sum_tot = sum_vec[length(sum_vec)]
  
  nPCs = c(nPCs,which(sum_vec > 0.01*percentage*sum_tot)[1])
  
}


PCA_cov(DT, weight_mat = wm, 
        Y = c(training_years,validation_years),
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
                       save_dir = GS_dir,mc_cores = mc_cores)

#########################################################
################ ECC & Schaake shuffle ##################
#########################################################

# they don't need setup, we just create directories:

ECC_dir = paste0(save_dir, "ECC/")
dir.create(ECC_dir, showWarnings = FALSE)

Schaake_dir = paste0(save_dir, "Schaake/")
dir.create(Schaake_dir, showWarnings = FALSE)

#### time, update script counter, save ####

time_s4 = proc.time() - time_s4

script_counter = 4

save.image(file = paste0(save_dir,"setup.RData"))


