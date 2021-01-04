
###############################################################################

#############  side script 5.1 - multivariate rank histograms  ################

###############################################################################

# This script compares the multivariate post-processing methods by means of their multivariate rank histograms.
# The rank histograms are for the entire region contained in DT (not shown in the paper).


#### setting up ######

rm(list = ls())

options(max.print = 1e3)

library(pp.sst)

name_abbr = "Full" 

save_dir = file.path('~','SST','Derived', name_abbr)

load(file = paste0(save_dir,"setup.RData"))

time_s51 = proc.time()

##################

brks = 6 # number of bins for the rank histogram - don't change that

### PCA ###

# get forecast:
load(paste0(PCA_dir,"fc_mc.RData"))

rks_PCA_mc = mv_rank_hist_new(PCA_fc_mc, fc_ens_size = fc_ens_size,
                              mc_cores = mc_cores,
                              breaks = brks, mn = "PCA_mc rank histograms",
                              save_pdf = TRUE, plot_dir = plot_dir, file_name = "rank_histo_PCA_mc")



load(paste0(PCA_dir,"fc_ac.RData"))

rks_PCA_ac = mv_rank_hist_new(PCA_fc_ac, fc_ens_size = fc_ens_size, 
                              mc_cores = mc_cores,
                              breaks = brks, mn = "PCA_ac rank histograms",
                              save_pdf = TRUE, plot_dir = plot_dir, file_name = "rank_histo_PCA_ac")


### geostationary ###

# get forecast:
load(paste0(GS_dir,"fc.RData"))

rks_GS = mv_rank_hist_new(GS_fc, fc_ens_size = fc_ens_size, 
                          mc_cores = mc_cores,
                          breaks = brks, mn = "GS rank histograms",
                          save_pdf = TRUE, plot_dir = plot_dir, file_name = "rank_histo_GS")


### ECC ###

# get forecast:
load(paste0(ECC_dir,"fc.RData"))

rks_ECC = mv_rank_hist_new(ECC_fc, fc_ens_size = ens_size, 
                           mc_cores = mc_cores,
                           breaks = brks, mn = "ECC rank histograms",
                           save_pdf = TRUE, plot_dir = plot_dir, file_name = "rank_histo_ECC")

### Schaake ###

# get forecast:
load(paste0(Schaake_dir,"fc.RData"))

# subsample for fair comparison: draw maximum ensemble size n such that n+1 is divisible by brks-1

rem = length(training_years) %% (brks-1)

ens_size_Schaake = length(training_years) - (rem + 1)
ens_mem_Schaake = sample(length(training_years),ens_size_Schaake)

fc_Schaake = paste0('fc',ens_mem_Schaake)

Schaake_fc[,paste0('fc',1:ens_size_Schaake,'temp'):=.SD,.SDcols = fc_Schaake]

Schaake_fc[,paste0('fc',1:length(training_years)):= NULL,]

setnames(Schaake_fc,old = paste0('fc',1:ens_size_Schaake,'temp'),new = paste0('fc',1:ens_size_Schaake))

#####################

rks_Schaake = mv_rank_hist_new(Schaake_fc, fc_ens_size = ens_size_Schaake, 
                           mc_cores = mc_cores,
                           breaks = brks, mn = "Schaake shuffle rank histograms",
                           save_pdf = TRUE, plot_dir = plot_dir, file_name = "rank_histo_Schaake")



time_s51 = proc.time() - time_s51

save.image(file = paste0(save_dir,"setup.RData"))
