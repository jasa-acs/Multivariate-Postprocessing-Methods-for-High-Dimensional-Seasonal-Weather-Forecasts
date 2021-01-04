

##### setting up ######

rm(list = ls())

time_s31 = proc.time()

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = 'Full/lv' 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))


######################

which.min(sc_sma[,lapply(.SD,mean),.SDcols = paste0('err',1:30)])


###### finding optimal way of SMA bias correction for each year in the validation period ######

opt_par_sma = data.table(year = validation_years,method = NA_character_,par = NA_real_)


for(y in validation_years)
{
    opt_par_sma[year == y, method := 'sma']
    opt_par_sma[year == y,par := msc_sma[year == y,min_l]]
}


# bias correction year by year by SMAs

for(y in validation_years)
{
  print(y) 
  temp = bias_correct(dt = DT,
                      method = opt_par_sma[year == y, method],
                      par_1 = opt_par_sma[year == y,par])[year == y,]
  DT[year == y,][,Bias_Est_sma := temp[,Bias_Est]][,SST_hat_sma := temp[,SST_hat]]
}

rm(temp)

### conduct permutation test SMA vs EMA ###



#####################################################
### permutation tests for MSE EMA vs SMA ###
#####################################################

perm_test_dt_MA = DT[year %in% validation_years & month %in% months,.(year,month,SST_bar,SST_hat,SST_hat_sma)]

# getting MSEs

perm_test_dt_MA[,MSE_ema := (SST_bar - SST_hat)^2]
perm_test_dt_MA[,MSE_sma := (SST_bar - SST_hat_sma)^2]


### permutation test for MSE_ma ~ MSE_lr_bb ###

pt_MAs = permutation_test_difference(na.omit(perm_test_dt_MA[,MSE_ema]),na.omit(perm_test_dt[,MSE_sma]), N = N )

pdf(paste0(plot_dir,'Perm_test_MSE.pdf'))

rr = max(abs(1.1*pt_MSE$d_bar),abs(1.1*pt_MSE$D))
rr = c(-rr,rr)

hist(pt_MSE$D, xlim = rr,breaks = 10,
     xlab = '', main = latex2exp::TeX('MSE permutation test for EMA vs. $NGR_{m,s}$'))

abline(v = pt_MSE$d_bar,col = 'red')


qq = quantile(pt_MSE$D,c(0.05))
abline(v = qq,lty = 2)

qq_2 = quantile(pt_MSE$D,c(0.01))
abline(v = qq_2,lty = 3)


dev.off()

# permutation test for MSE_ma vs MSE_lr_bb, averaged over the globe

ptbm = perm_test_dt_MSE[,.('MSE_ma' = mean(MSE_ma,na.rm = TRUE),'MSE_lr_bb' = mean(MSE_lr_bb,na.rm = TRUE)),by = .(year,month)]

pt_MSE_bm = permutation_test_difference(ptbm[,MSE_ma],ptbm[,MSE_lr_bb], N = 10000  )

pdf(paste0(plot_dir,'Perm_test_glob_mean_MSE.pdf'))

rr = max(abs(1.1*pt_MSE_bm$d_bar),abs(1.1*pt_MSE_bm$D))
rr = c(-rr,rr)

hist(pt_MSE_bm$D, xlim = rr,breaks = 20,
     xlab = '', main = latex2exp::TeX('globally averaged MSE permutation test for EMA vs. $NGR_{m,s}$'))

abline(v = pt_MSE_bm$d_bar,col = 'red')

qq = quantile(pt_MSE_bm$D,c(0.05))
abline(v = qq,lty = 2)

qq_2 = quantile(pt_MSE_bm$D,c(0.01))
abline(v = qq_2,lty = 3)


dev.off()





