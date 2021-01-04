

####################################################
### This script generates the plots of section 2 ###
####################################################

# this script is not optimized to be computationally efficient. For any plot, it loads a saved workspace image to retain the corresponding data set.
# It is adviced not to run the entire script if possible, but to use this to reconstruct and or modify single plots by running the corresponding section below.



rm(list = ls())

time_s4 = proc.time()

setwd("~/NR/SFE")
options(max.print = 1e3)

library(data.table)

devtools::load_all()

########### set parameters for plots #######################

par('cex' = 0.75, 'cex.lab' = 0.6,'cex.axis' = 0.6)

plot_dir0 = './figures/paper/'
dir.create(plot_dir0,showWarnings = FALSE)


Lat_res = c(-75,80) # Latitude restrictions for area plots in order to exclude the polar regions


########################################################


############# Plot of normalized SST #######################


name_abbr = "Full"

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))

# reset plot directory

plot_dir = plot_dir0

# year and month
y = 2016
m = 5


DT = DT[Lat %between% Lat_res]

dt_ym = DT[year == y & month == m]

# get sample mean and sample variance

test_dt = DT[month == m & year < y & !is.na(SST_hat),][,.('clim' = mean(SST_bar,na.rm = TRUE),'clim_sd' = sd(SST_bar,na.rm = TRUE)),by = grid_id]

n_years = dt[month == m & year < y & !is.na(SST_hat),][grid_id == min(grid_id),.N]

temp = dt_ym[!is.na(SST_hat),][,'clim' := test_dt[,clim]][,'clim_sd' := test_dt[, clim_sd]]

dt_ym2 = merge(temp,dt_ym,by = colnames(dt_ym),all.y = TRUE)


# bound standard deviations away from 0

dt_ym2[clim_sd < 0.1,clim_sd := 0.1]

# get standardized SST and plot
dt_ym2[ ,'SST_stan' := (SST_bar - clim)/clim_sd]

plot_smooth(dt_ym2,'SST_stan',rr = c(-5,5),
            mn = 'normalized sea surface temperature in May 2016',
            pixels = 512,
            save_pdf = TRUE, save_dir = plot_dir,file_name = 'sst_nor_m05_y2016')





############################################################


#### plot of the first principal component ####


# get data

name_abbr = "NAO"

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))


# reset plot directory

plot_dir = plot_dir0


# year and month

y = 2016
m = 7

# get the singular value decomposition of weight matrix * sample covariance matrix

wm = weight_mat(DT,L = 7000)


num_loc = DT[year == min(year) & month == min(month)][!(is.na(SST_bar) | is.na(SST_hat)) ,.N]

train_years = DT[month == m][year < y & year > min(year),][,unique(year)]

data_mat = matrix(DT[month == m][!(is.na(SST_bar) | is.na(SST_hat)) & year %in% train_years,
                                 SST_bar - SST_hat],
                  nrow = num_loc)

sam_cov_mat = 1/length(train_years) * data_mat %*% t(data_mat)


sin_val_dec_1 = svd(wm * sam_cov_mat)

# some data preparation

dt_water = DT[!(is.na(Ens_bar) | is.na(SST_bar))]

SD_cols = c("Lon","Lat","grid_id","month","year","YM",
            "SST_hat","SST_bar","Ens_bar","Bias_Est","var_bar","SD_hat")
SD_cols = SD_cols[which(SD_cols %in% colnames(dt))]
fc_water <- na.omit( dt_water[,.SD,.SDcols = SD_cols])

dt_ym = fc_water[month == m & year == y,]


# complement dt_ym by principal components

for(i in 1:30)
{
  temp =  sin_val_dec_1$u[,i]

  dt_ym[,paste0('PC',i):= temp]

  dt_test = rbindlist(list(dt_ym,dt[year == y & month ==m][is.na(Ens_bar) | is.na(SST_bar),.SD,.SDcols = SD_cols]), fill = TRUE)

}


# generate plot for restricted area

dt_test_new = dt_test[Lon >= -20 & Lat >= 50,]

rr = dt_test_new[,range(PC1,na.rm = TRUE)]
rr = c(-(max(abs(rr))),(max(abs(rr))))

plot_smooth(dt_test_new,paste0('PC1'),mn = '1st PC of covariance matrix for June',save_pdf = TRUE,file_name = '1stPCtapered',save_dir = plot_dir,xlab = '',ylab = '' )


