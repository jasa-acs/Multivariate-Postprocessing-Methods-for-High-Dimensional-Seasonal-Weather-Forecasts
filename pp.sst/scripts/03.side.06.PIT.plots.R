#######################################################################################

###################  side script 3.6 - PIT plots  ######################

#######################################################################################

# This script generates plots of the mean and standard deviation of the marginal probability integral transforms 
# for the marginally calibrated forecast, both averaged over the entire year (paper) and month-by-month (supplementary material).
#
# This requires a previous run of the side script 03.side.01.compare.univariate.bias.correction.R
# 
##### setting up ######

rm(list = ls())

options(max.print = 1e3)

library(pp.sst)

name_abbr = "Full" 

save_dir = file.path('~','SST','Derived', name_abbr)

load(file = paste0(save_dir,"setup.RData"))


######## get the distribution fct. of a censored normal #############

# value, mean and sd need to be vectors of equal length.
# returns F(value), where F is dist. fct. of a normal with parameters mean and sd, censored at trc_value

dist_fun_tn = function(value, mean, sd, trc_value = -1.79){ 
  a=rep(0,times = length(value))
  na_loc = which(is.na(value) | is.na(sd) | sd == 0)
  trc_loc = which(value <= trc_value & sd > 0)
  nor_loc = which(value > trc_value & sd > 0)
  
  a[na_loc] = NA
  a[trc_loc] = runif(length(trc_loc), max = pnorm(trc_value, mean = mean[trc_loc], sd = sd[trc_loc]))
  a[nor_loc] = pnorm(value[nor_loc], mean = mean[nor_loc], sd = sd[nor_loc])
    
  return(a)
}

########### get PITs ###############

DT_calib_1 = DT[year %in% validation_years,.(Lon,Lat,year,month,grid_id,SST_bar,SST_hat,SD_hat,Ens_bar,Ens_sd)]

# PIT for marginally corrected forecast:
DT_calib_1[,"PIT_mc" := dist_fun_tn(SST_bar, mean = SST_hat, sd = SD_hat)]

DT_calib_1[,"PIT_mc_mean" := mean(PIT_mc), by = grid_id]
DT_calib_1[,"PIT_mc_sd" := sd(PIT_mc), by = grid_id]

# by month

DT_calib_1[,"PIT_mc_mean_bm" := mean(PIT_mc), by = .(grid_id,month)]
DT_calib_1[,"PIT_mc_sd_bm" := sd(PIT_mc), by = .(grid_id,month)]

########### plot ################
Lat_res = c(-75,80)

pdf(paste0(plot_dir,'PIT_mean.pdf'),width = 14)

  par('cex' = 2, 'cex.axis' = 0.75)
  plot_diagnostic(DT_calib_1[year == min(year) & Lat %between% Lat_res & month == min(month), .(Lon,Lat,PIT_mc_mean)],
                  rr = c(0,1),
                  mn = paste0("PIT mean"))
dev.off()


unif_sd = sqrt(1/12) #standard deviation of uniform distribution

pdf(paste0(plot_dir,'PIT_sd.pdf'),width = 14)
  
  par('cex' = 2, 'cex.axis' = 0.75)
  plot_diagnostic(DT_calib_1[year == min(year) & Lat %between% Lat_res & month == min(month), .(Lon,Lat,PIT_mc_sd)],
                  rr = c(0,2*unif_sd),set_white = unif_sd,
                  mn = paste0("PIT standard deviation"))
dev.off()

# month-by-month plot

for(mm in months)
{
pdf(paste0(plot_dir,'PIT_mean',mm,'.pdf'),width = 14)

par('cex' = 2, 'cex.axis' = 0.75)
plot_diagnostic(DT_calib_1[year == min(year) &  Lat %between% Lat_res & month == mm, .(Lon,Lat,PIT_mc_mean_bm)],
                rr = c(0,1),
                mn = paste0("PIT mean, month ",mm))
dev.off()


unif_sd = sqrt(1/12) #standard deviation of uniform distribution

pdf(paste0(plot_dir,'PIT_sd',mm,'.pdf'),width = 14)

par('cex' = 2, 'cex.axis' = 0.75)
plot_diagnostic(DT_calib_1[year == min(year) & Lat %between% Lat_res & month == mm, .(Lon,Lat,PIT_mc_sd_bm)],
                rr = c(0,2*unif_sd),set_white = unif_sd,
                mn = paste0("PIT standard deviation, month ",mm))
dev.off()
}




###### for estimated mean by linear regression ######

DT_calib_3 = DT[year %in% validation_years,.(Lon,Lat,year,month,grid_id,SST_bar,T_hat_lr_bb,SD_hat,Ens_bar,Ens_sd)]

# PIT for marginally corrected forecast:
DT_calib_3[,"PIT_mc" := dist_fun_tn(SST_bar, mean = T_hat_lr_bb, sd = SD_hat)]

DT_calib_3[,"PIT_mc_mean" := mean(PIT_mc), by = grid_id]
DT_calib_3[,"PIT_mc_sd" := sd(PIT_mc), by = grid_id]

# by month

DT_calib_3[,"PIT_mc_mean_bm" := mean(PIT_mc), by = .(grid_id,month)]
DT_calib_3[,"PIT_mc_sd_bm" := sd(PIT_mc), by = .(grid_id,month)]

########### plot ################
Lat_res = c(-75,80)

pdf(paste0(plot_dir,'PIT_mean_NGR.pdf'),width = 14)

par('cex' = 2, 'cex.axis' = 0.75)
plot_diagnostic(DT_calib_3[year == min(year) & Lat %between% Lat_res & month == min(month), .(Lon,Lat,PIT_mc_mean)],
                rr = c(0,1),
                mn = expression("PIT mean, NGR"['m,s']))
dev.off()

# month-by-month plot

for(mm in months)
{
  pdf(paste0(plot_dir,'PIT_mean_NGR',mm,'.pdf'),width = 14)
  
  par('cex' = 2, 'cex.axis' = 0.75)
  plot_diagnostic(DT_calib_3[year == min(year) &  Lat %between% Lat_res & month == mm, .(Lon,Lat,PIT_mc_mean_bm)],
                  rr = c(0,1),
                  mn = paste0("PIT mean NGR, month ", mm))
  dev.off()
  
}


save.image(file = paste0(save_dir,"setup.RData"))

