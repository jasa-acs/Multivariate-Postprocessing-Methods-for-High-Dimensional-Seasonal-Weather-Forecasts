

##################################################################
### This script generates all the plots contained in the paper ###
##################################################################

# This script generates all the plots shown in the paper, subject to having only the example data available that is provided with this R package.
# It requires previous runs of the following scripts
#   - all master scripts
#   - 03.side.01.compare.univariate.bias.correction.linear.methods.R
#   - 05.side.04.route.MVRHs.R

rm(list = ls())

options(max.print = 1e3)

library(pp.sst)

name_abbr = "pp.sst/Full"

save_dir = file.path('~','SST','Derived', name_abbr) # where did you store your run?

plot_dir0 = '/nr/project/stat/SFE/figures/test/' # where do you wnat to save your Figures? Adjust!
  dir.create(plot_dir0,recursive = T,showWarnings = FALSE)

load(file = paste0(save_dir,"setup.RData"))

########### set parameters for plots #######################

par('cex' = 1, 'cex.lab' = 1,'cex.axis' = 1)




Lat_res = c(-75,80) # Latitude restrictions for area plots in order to exclude the polar regions

########################################################
############# Figure 1: Plot of normalized SST #######################

# year and month, these are different than in the paper in order to run on the example data provided in this package
y = 2016
m = 2

# get out-of-sample mean and sd

DT = DT[Lat %between% Lat_res]

mu = DT[year < y & month == m,mean(SST_bar),.(Lon,Lat)][[3]]
sd = DT[year < y & month == m,sd(SST_bar),.(Lon,Lat)][[3]]

dt_ym = DT[year == y & month == m][,SST_stan := (SST_bar-mu)/sd]

# 'proper' scaling of the plot
scale_factor = (dt_ym[,max(Lon)]-dt_ym[,min(Lon)])/(dt_ym[,max(Lat)]-dt_ym[,min(Lat)])

pdf(paste0(plot_dir0,'Fig1.pdf'),height = 7,width = scale_factor * 7)
  plot_smooth(dt_ym,'SST_stan',mn = paste0('Standardized SST anomaly, month ',m))
dev.off()

################################################################################
############# Figure 2: Plot of MSE for SMA and EMA by parameter ###############

# more or less identical to script 02.side.01.RMSEplots.R

# get data for last year
y = max(validation_years)


row_sma = msc_sma[year == y,-1,with = FALSE]
row_ema = msc_ema[year == y,-1,with = FALSE]

y_range = range(list(row_sma[,-'min_l',with = FALSE],row_ema[,-'min_a',with = FALSE]))

y_range = sqrt(y_range)

## plot for sma ##

pdf(paste0(plot_dir0,"Fig2_1.pdf"))
  plot(x = win_length,
       y = as.vector(sqrt(row_sma[,-c('min_MSE','min_l'),with = FALSE])),
       ylim = y_range,
       type = "b",
       col = "blue",
       main = paste0("RMSE for bias correction by SMA"),
       xlab = "window length",
       ylab = "RMSE"
  )

  # highlight minimum and add minimum reference line
  abline(h = row_sma[,sqrt(min_MSE)], lty = "dashed", col = adjustcolor("blue",alpha = .5))

  points(x = row_sma[,min_l],
         y = row_sma[,sqrt(min_MSE)],
         col = "blue",
         bg = "blue",
         pch = 21)
dev.off()

## plot for ema ##

pdf(paste0(plot_dir0,"Fig2_2.pdf"))
plot(x = par_vec,
     y = sqrt(row_ema[,-c('min_MSE','min_a'),with = FALSE]),
     ylim = y_range,
     type = "b",
     col = "blue",
     main = paste0("RMSE for bias correction by EMA"),
     xlab = "scale parameter",
     ylab = "RMSE"
)

# highlight minimum and add minimum reference line
abline(h = row_ema[,sqrt(min_MSE)], lty = "dashed", col = adjustcolor("blue",alpha = .5))

points(x = row_ema[,min_a],
       y = row_ema[,sqrt(min_MSE)],
       col = "blue",
       bg = "blue",
       pch = 21)
dev.off()


########################################################
############# Figure 3: Plot of PIT mean ###############

# more or less identical to script 03.side.06.PITplots.R

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

# plot

pdf(paste0(plot_dir0,'Fig3_1.pdf'),width = 14)

  par('cex' = 2, 'cex.axis' = 0.75)
  plot_diagnostic(DT_calib_1[year == min(year) & month == min(month), .(Lon,Lat,PIT_mc_mean)],
                  rr = c(0,1),
                  mn = paste0("PIT mean"))
dev.off()


###### for estimated mean by linear regression ######
# this part requires a previous run of side script 03.01.

DT_calib_3 = DT[year %in% validation_years,.(Lon,Lat,year,month,grid_id,SST_bar,T_hat_lr_bb,SD_hat,Ens_bar,Ens_sd)]

# PIT for marginally corrected forecast:
DT_calib_3[,"PIT_mc" := dist_fun_tn(SST_bar, mean = T_hat_lr_bb, sd = SD_hat)]

DT_calib_3[,"PIT_mc_mean" := mean(PIT_mc), by = grid_id]
DT_calib_3[,"PIT_mc_sd" := sd(PIT_mc), by = grid_id]

# by month

DT_calib_3[,"PIT_mc_mean_bm" := mean(PIT_mc), by = .(grid_id,month)]
DT_calib_3[,"PIT_mc_sd_bm" := sd(PIT_mc), by = .(grid_id,month)]

########### plot ################

pdf(paste0(plot_dir0,'Fig3_2.pdf'),width = 14)
  par('cex' = 2, 'cex.axis' = 0.75)
  plot_diagnostic(DT_calib_3[year == min(year) & month == min(month), .(Lon,Lat,PIT_mc_mean)],
                  rr = c(0,1),
                  mn = expression("PIT mean, NGR"['m,s']))
dev.off()

#######################################
#### Fig 4: PIT standard deviation ####

unif_sd = sqrt(1/12) #standard deviation of uniform distribution

pdf(paste0(plot_dir0,'Fig4.pdf'),width = 14)
  par('cex' = 2, 'cex.axis' = 0.75)
  plot_diagnostic(DT_calib_1[year == min(year) & month == min(month), .(Lon,Lat,PIT_mc_sd)],
                  rr = c(0,2*unif_sd),set_white = unif_sd,
                  mn = paste0("PIT standard deviation"))
dev.off()

####################################################
#### Fig 5: Example residual and shipping route ####

# choose month and year for observed residual
m = 2
y = 2016

#in order to add route, we need to set the stretch parameter manually:

pdf(paste0(plot_dir0,'Fig5.pdf'),width = 7 * scale_factor,height = 7)

  plot_smooth(DT[year == y & month == m,.(Lon,Lat,SST_bar - SST_hat)],
              mn = paste0('Forecast residual, ',m,'/',y),
              rr = c(-3,3),
              brks = c(-3,-1.5,0,1.5,3))

  # add Bordeaux - Norfolk shipping route

  Bordeaux = c(-0.57,44.8)
  Norfolk = c(-76.3,36.9)

  p1 = data.table(Lon = Bordeaux[1], Lat = Bordeaux[2], Loc = 'Bordeaux')
  p2 = data.table(Lon = Norfolk[1], Lat = Norfolk[2], Loc = "Norfolk")

  cities = rbindlist(list(p1,p2))

  points(cities[,Lon],cities[,Lat],col="black", cex=2, pch=20)

  par('cex' = save_cex)

  # Connection between Bordeaux and Norfolk
  inter <- geosphere::gcIntermediate(Bordeaux, Norfolk, n=100, addStartEnd=TRUE, breakAtDateLine=F)
  lines(inter, col="black", lwd=2)

dev.off()


################################
#### Fig 6: Rank histograms ####


brks = 6 # how many bins do you want?

# range
rr = c(0,20)

titles = c(latex2exp::TeX('$\\widehat{\\Sigma}^{mc}$'),
           latex2exp::TeX('$\\widehat{\\Sigma}^{ac}$'),
           'GS',
           'ECC',
           'Schaake')

### average rank histograms ###

pdf(paste0(plot_dir0,'Fig6_1.pdf'),width = 35,height = 8)

  par(oma = c(1,1,1,1), mfrow=c(1,5), mar=c(2,1,2,1) )
  par('cex' = 3,'cex.axis' = 0.75)

  ind = 1
  for(mod in mod_vec)
  {

    ens_size = fc_ens_size
    if(mod == 'ECC') ens_size = 9
    if(mod == 'Schaake') ens_size = length(training_years)


    temp = get(paste0('rks_',mod,'_route'))

    # for the Schaake shuffle the number of bins does not divide the number of ranks:
    if(mod == 'Schaake')
    {
      temp = data.table(YM = rks_Schaake_route[,YM],av.rk.obs = hist_transform(rks_Schaake_route[,av.rk.obs]))
    }

    rhist_dt(temp[,.(YM,av.rk.obs)],
             max_rk = ens_size +1,
             breaks = brks,
             hist_xlab = "",
             hist_ylim = rr)
    title(main = titles[ind],line = -1)

    ind = ind + 1
  }

  par('cex' = 3)
  title(main = 'Average rank',outer = TRUE,line = -1)

dev.off()


### band depth rank histograms ###

pdf(paste0(plot_dir0,'Fig6_2.pdf'),width = 35,height = 8)

  par(oma = c(1,1,1,1), mfrow=c(1,5), mar=c(2,1,2,1) )
  par('cex' = 3,'cex.axis' = 0.75)

  ind = 1
  for(mod in mod_vec)
  {

    ens_size = fc_ens_size
    if(mod == 'ECC') ens_size = 9
    if(mod == 'Schaake') ens_size = length(training_years)


    temp = get(paste0('rks_',mod,'_route'))

    # for the Schaake shuffle the number of bins does not divide the number of ranks:
    if(mod == 'Schaake')
    {
      temp = data.table(YM = rks_Schaake_route[,YM],bd.rk.obs = hist_transform(rks_Schaake_route[,bd.rk.obs]))
    }

    rhist_dt(temp[,.(YM,bd.rk.obs)],
             max_rk = ens_size +1,
             breaks = brks,
             hist_xlab = "",
             hist_ylim = rr)
    title(main = titles[ind],line = -1)

    ind = ind + 1
  }

  par('cex' = 3)
  title(main = 'Band depth rank',outer = TRUE,line = -1)

dev.off()

