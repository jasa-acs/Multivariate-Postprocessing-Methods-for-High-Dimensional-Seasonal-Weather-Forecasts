

###################################################################################
### This script generates the plots for the supplementary material of the paper ###
###################################################################################



rm(list = ls())

time_s4 = proc.time()

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)


########### set parameters for plots #######################

par('cex' = 0.75, 'cex.lab' = 0.6,'cex.axis' = 0.6)

plot_dir0 = './figures/paper/'
dir.create(plot_dir0,showWarnings = FALSE)


Lat_res = c(-75,80) # Latitude restrictions for area plots in order to exclude the polar regions


#######################################################

#########  Permutation test CRPS and MSE  #############

# for the spatially averaged permutation test:
ptbm = perm_test_dt_MSE[,.('MSE_ma' = mean(MSE_ma,na.rm = TRUE),'MSE_lr_bb' = mean(MSE_lr_bb,na.rm = TRUE)),by = .(year,month)]

pt_MSE_bm = permutation_test_difference(ptbm[,MSE_ma],ptbm[,MSE_lr_bb], N = 10000  )


# permutation test for CRPS_ma ~ CRPS_lr_bb, averaged over the globe

ptbm = perm_test_dt_CRPS[,.('CRPS_ma' = mean(CRPS_ma,na.rm = TRUE),'CRPS_lr_bb' = mean(CRPS_lr_bb,na.rm = TRUE)),by = .(year,month)]

pt_CRPS_bm = permutation_test_difference(ptbm[,CRPS_ma],ptbm[,CRPS_lr_bb], N = 10000  )



pdf(paste0(plot_dir,'Perm_test_glob_mean_MSE_CRPS.pdf'),width = 15,height = 7)

par('mfrow' = c(1,2))

#MSE

rr = max(abs(1.1*pt_MSE_bm$d_bar),abs(1.1*pt_MSE_bm$D))
rr = c(-rr,rr)

hist(pt_MSE_bm$D, xlim = rr,breaks = 20,
     main = latex2exp::TeX(paste0('mean estimation: $NGR_{m,s}$ vs. EMA$')),
     col = rgb(t(col2rgb('gray')/255),alpha = 0.5,maxColorValue = 1),
     border = rgb(t(col2rgb('gray')/255),alpha = 0.5,maxColorValue = 1),
     xlab = '',
     probability = TRUE,
     ylab = '',
     axes = FALSE)
axis(side = 1)

abline(v = pt_MSE_bm$d_bar,col = 'red')

qq = quantile(pt_MSE_bm$D,0.05)

abline(v = qq,lty = 2)

# CRPS

rr = max(abs(1.1*pt_CRPS_bm$d_bar),abs(1.1*pt_CRPS_bm$D))
rr = c(-rr,rr)

hist(pt_CRPS_bm$D, xlim = rr,breaks = 20,
     main = latex2exp::TeX(paste0('variance estimation: $NGR_{m,s}$ vs. EMA$')),
     col = rgb(t(col2rgb('gray')/255),alpha = 0.5,maxColorValue = 1),
     border = rgb(t(col2rgb('gray')/255),alpha = 0.5,maxColorValue = 1),
     xlab = '',
     probability = TRUE,
     ylab = '',
     axes = FALSE)
axis(side = 1)

abline(v = pt_CRPS_bm$d_bar,col = 'red')

qq = quantile(pt_CRPS_bm$D,0.05)

abline(v = qq,lty = 2)


dev.off()

#############################################################

################## CRPS by parameter ########################

# get data for last year
y = max(validation_years)


row_sma = msc_sd_sma[year == y,-1,with = FALSE]
row_ema = msc_sd_ema[year == y,-1,with = FALSE]

y_range = range(list(row_sma[,-'min_l',with = FALSE],row_ema[,-'min_a',with = FALSE]))  

pdf(paste0(plot_dir,"/mean_CRPS_sd.pdf"),width = 15)

par(oma = c(1,1,1,1), mfrow=c(1,2), mar=c(2,1,2,1) )
par('cex' = 1.4,'cex.axis' = 0.75)

## plot for sma ##

plot(x = win_length, 
     y = as.vector(row_sma[,-c('min_crps','min_l'),with = FALSE]),
     ylim = y_range,
     type = "b",
     col = "blue",
     main = paste0("CRPS for variance estimation by SMA"),
     xlab = "window length",
     ylab = "CRPS")

# highlight minimum and add minimum reference line 
abline(h = row_sma[,min_crps], lty = "dashed", col = adjustcolor("blue",alpha = .5))

points(x = row_sma[,min_l],
       y = row_sma[,min_crps],
       col = "blue",
       bg = "blue",
       pch = 21)


## plot for ema ##

plot(x = par_vec, 
     y = row_ema[,-c('min_crps','min_a'),with = FALSE],
     ylim = y_range,
     type = "b",
     col = "blue",
     main = paste0("CRPS for variance estimation by EMA"),
     xlab = "scale parameter",
     ylab = "CRPS"
)

# highlight minimum and add minimum reference line 
abline(h = row_ema[,min_crps], lty = "dashed", col = adjustcolor("blue",alpha = .5))

points(x = row_ema[,min_a],
       y = row_ema[,min_crps],
       col = "blue",
       bg = "blue",
       pch = 21)

dev.off()


####################################################

######## Plots for multivariate calibration ########

####################################################

name_abbr = "NAO/lv/2" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))


# reset plot directory

plot_dir = plot_dir0


####################################################################

##################### plot example residuals #######################

m = 6
y = 2016

#in order to add route, we need to set the stretch parameter manually:

temp = DT[year == y & month == m,.(Lon,Lat,SST_bar - SST_hat)]

Lons = unique(temp[,Lon])
Lats = unique(temp[,Lat])

n_lon = length(Lons)
n_lat = length(Lats)

save_cex = par('cex')


pdf(paste0(plot_dir,'Example_res.pdf'),width = 7,height = 7 * n_lat/n_lon)

par('cex' = save_cex)
plot_smooth(DT[year == y & month == m,.(Lon,Lat,SST_bar - SST_hat)], 
            mn = latex2exp::TeX('observed forecast residual June 2016'), 
            rr = c(-3,3),
            pixels = 512,
            save_pdf = FALSE,
            save_dir = plot_dir,
            file_name = 'Example_res'
)  

#add points and shipping route

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


####################################################################

################### plot simulated residuals #######################

# get data:

load(file = paste0(PCA_dir,"fc_mc.RData"))

load(file = paste0(PCA_dir,"fc_ac.RData"))

load(file = paste0(GS_dir,"fc.RData"))

load(file = paste0(ECC_dir,"fc.RData"))

# get plotting function that doesn't add labels, 
# just a few tweeks from plot_smooth, but plot_smooth wraps image.plot, which always adds label, and this wraps image:

plot_smooth_nl = function( dt, var = colnames(dt)[3], mn = var, rr = NULL,...,
                           theta = 0.5, pixels = 256,
                           col_scheme = "bwr", set_white = NULL,
                           xlab = "", ylab = "",
                           save_pdf = FALSE, save_dir = "./figures/", file_name = "diag_plot", stretch_par = NULL)
{
  # prepare data table
  
  if("year" %in% colnames(dt))
  {
    if("month" %in% colnames(dt))
    {
      dt = dt[year == min(year) & month == min(month),.SD,.SDcols = c('Lon','Lat',var)][order(Lat,Lon)]
    } else {
      dt = dt[month == min(month),.SD,.SDcols = c('Lon','Lat',var)][order(Lat,Lon)]
    }
  } else {
    dt = dt[,.SD,.SDcols = c('Lon','Lat',var)][order(Lat,Lon)]
  }
  
  
  #--- create image ---
  
  x = dt[,.(Lon,Lat)]
  setnames(x,c("Lon","Lat"), c("lon","lat"))
  
  Lons = unique(dt[,Lon])
  Lats = unique(dt[,Lat])
  
  n_lon = length(Lons)
  n_lat = length(Lats)
  
  A = matrix(dt[[3]],  n_lon, n_lat)
  
  im_0 = fields::image.smooth(fields::as.image(A,x = x,nx = pixels,ny = pixels),theta = theta)
  
  ## Find the points that fall over land
  
  if(!exists("wrld_simpl")) data(wrld_simpl, package = 'maptools') 
  
  all_loc = expand.grid(lat = im_0$x,lon = im_0$y)
  pts <- sp::SpatialPoints(all_loc, proj4string=sp::CRS(sp::proj4string(wrld_simpl)))
  ii <- !is.na(sp::over(pts, wrld_simpl)$FIPS)
  im_0$z[ii] = NA
  
  # --- fix range of plot and fill in values for points out of range ---
  
  if(is.null(rr))  rr = range(im_0$z,na.rm=TRUE)
  if(!is.null(rr)){
    im_0$z[im_0$z< min(rr)] = min(rr)
    im_0$z[im_0$z> max(rr)] = max(rr)
  }
  
  # --- scaling and colors ---
  
  brk = seq(rr[1],rr[2],length = 500)
  brk.ind = round(seq(1,length(brk),length = 10))
  brk.lab = round(brk[brk.ind],2)
  brk.at = brk[brk.ind]
  
  if(col_scheme == "bwr"){
    if(is.null(set_white)){
      color <- fields::designer.colors(n=length(brk)-1, col = c("darkblue","white","darkred"))
    }else{
      zero.ind = min(which(brk > set_white))/length(brk)
      color <- fields::designer.colors(n=length(brk)-1, col = c("darkblue","white","darkred"), x = c(0,zero.ind,1))
    }
  }
  if(col_scheme == "wr"){
    color <- fields::designer.colors(n=length(brk)-1, col = c("white","darkred"))
  }
  if(col_scheme == "wb"){
    color <- fields::designer.colors(n=length(brk)-1, col = c("white","blue"))
  }
  
  # color NAs grey
  
  newz.na <- rr[2]+(rr[2]-rr[1])/length(color) # new z for NA
  im_0$z[which(is.na(im_0$z))] <- newz.na 
  rr[2] <- newz.na # extend the range to include the new value 
  color <- c(color, 'gray') # extend the color range by gray
  brk = c(brk,rr[2]) # extend the vector of breaks
  
  #--- plotting ---
  
  if(save_pdf) 
  {
    if (is.null(stretch_par)) stretch_par = n_lat/n_lon
    
    par_0 = par() # allow to set par manually before calling the function
    
    pdf(paste0(save_dir,file_name,".pdf"),width = 7,height = stretch_par * 7)
    
    suppressWarnings(par(par_0))
  }
  
  par(mar = c(2,2,2,2))
  
  image(im_0,
        zlim=rr, main = mn,...,
        xlim = range(Lons), xlab=xlab,
        ylim = range(Lats), ylab=ylab,
        col=color)
  
  # add world map
  
  maps::map("world", add = TRUE)
  
  if(save_pdf) dev.off()
  
}

### plotting the slightly hand-picked examples ###

pdf(paste0(plot_dir,'Example_res_sim.pdf'),width = 15)

  par('mfrow' = c(2,2),'cex' = 1.25)
  
  plot_smooth_nl(PCA_fc_mc[year == y & month == m,.(Lon,Lat,fc6 - SST_hat)],
                 mn = latex2exp::TeX('simulated forecast residual, $\\widehat{\\Sigma}^{mc}$ '),
                 pixels = 512,
                 axes = FALSE,
                 rr = c(-3,3))
  
  
  plot_smooth_nl(PCA_fc_ac[year == y & month == m,.(Lon,Lat,fc11 - SST_hat)],
                 mn = latex2exp::TeX('simulated forecast residual, $\\widehat{\\Sigma}^{ac}$ '),
                 pixels = 512,
                 axes = FALSE,
                 rr = c(-3,3))
  
  plot_smooth_nl(GS_fc[year == y & month == m,.(Lon,Lat,fc1 - SST_hat)],
                 mn = latex2exp::TeX('simulated forecast residual, GS '),
                 pixels = 512,
                 axes = FALSE,
                 rr = c(-3,3))
  
  
  plot_smooth_nl(ECC_fc[year == y & month == m,.(Lon,Lat,fc1 - SST_hat)],
                 mn = latex2exp::TeX('simulated forecast residual, ECC '),
                 pixels = 512,
                 axes = FALSE,
                 rr = c(-3,3))


dev.off()

###### plotting the first four examples for  each method - no hand-picking ######

pdf(paste0(plot_dir,'Example_res_sim_PCA_mc.pdf'),width = 15)

  par('mfrow' = c(2,2),'cex' = 1.25, 'oma' = c(0,0,2,0))
  
  
  for(i in 1:4)
  {
    fc_name = paste0('fc',i)
      plot_smooth_nl(PCA_fc_mc[year == y & month == m,.(Lon,Lat,.SD - SST_hat),.SDcols = fc_name],
                   mn = "",
                   pixels = 256,
                   axes = FALSE,
                   rr = c(-3,3))
  
  }
  title(latex2exp::TeX('simulated forecast residuals, $\\widehat{\\Sigma}^{mc}$ '),outer=TRUE)
dev.off()


pdf(paste0(plot_dir,'Example_res_sim_PCA_ac.pdf'),width = 15)
  
  par('mfrow' = c(2,2),'cex' = 1.25, 'oma' = c(0,0,2,0))
  
  
  for(i in 1:4)
  {
    fc_name = paste0('fc',i)
    plot_smooth_nl(PCA_fc_ac[year == y & month == m,.(Lon,Lat,.SD - SST_hat),.SDcols = fc_name],
                   mn = "",
                   pixels = 256,
                   axes = FALSE,
                   rr = c(-3,3))
    
  }
  title(latex2exp::TeX('simulated forecast residuals, $\\widehat{\\Sigma}^{ac}$ '),outer=TRUE)
dev.off()



pdf(paste0(plot_dir,'Example_res_sim_GS.pdf'),width = 15)
  
  par('mfrow' = c(2,2),'cex' = 1.25, 'oma' = c(0,0,2,0))
  
  
  for(i in 1:4)
  {
    fc_name = paste0('fc',i)
    plot_smooth_nl(GS_fc[year == y & month == m,.(Lon,Lat,.SD - SST_hat),.SDcols = fc_name],
                   mn = "",
                   pixels = 256,
                   axes = FALSE,
                   rr = c(-3,3))
    
  }
  title(latex2exp::TeX('simulated forecast residuals, GS '),outer=TRUE)
dev.off()


pdf(paste0(plot_dir,'Example_res_sim_ECC.pdf'),width = 15)
  
  par('mfrow' = c(2,2),'cex' = 1.25, 'oma' = c(0,0,2,0))
  
  
  for(i in 1:4)
  {
    fc_name = paste0('fc',i)
    plot_smooth_nl(ECC_fc[year == y & month == m,.(Lon,Lat,.SD - SST_hat),.SDcols = fc_name],
                   mn = "",
                   pixels = 256,
                   axes = FALSE,
                   rr = c(-3,3))
    
  }
  title(latex2exp::TeX('simulated forecast residuals, ECC '),outer=TRUE)
dev.off()


# make room again:
rm(PCA_ac_fc,PCA_mc_fc,GS_fc,ECC_fc)
gc()


#####################################################

#### plot permutation tests for variogram scores ####


### permutation tests ###

N = 20000

mod_com_ls = list(c('PCA_mc' , 'PCA_ac'),
                  c('PCA_mc' , 'GS'),
                  c('PCA_mc' , 'ECC'),
                  c('PCA_ac' , 'GS'),
                  c('PCA_ac' , 'ECC'),
                  c('ECC'     , 'GS'))

mod_com_names = list(c('$\\widehat{\\Sigma}^{mc}$ ', '$\\widehat{\\Sigma}^{ac}$ '),
                     c('$\\widehat{\\Sigma}^{mc}$ ', 'GS'),
                     c('$\\widehat{\\Sigma}^{mc}$ ', 'ECC'),
                     c('$\\widehat{\\Sigma}^{ac}$ ', 'GS'),
                     c('$\\widehat{\\Sigma}^{ac}$ ', 'ECC'),
                     c('ECC','GS'))

pdf(paste0(plot_dir,'perm_test_vs.pdf'),width = 24,height = 15)

par(oma = c(1,1,1,1), mfrow=c(2,3), mar=c(2,1,2,1) )

par('cex' = 2.5,'cex.axis' = 0.75)

for(ind in 1:6)
{ mod_com = mod_com_ls[[ind]]
mod_names = mod_com_names[[ind]]  
mod1 = mod_com[1]
mod2 = mod_com[2]

vs_1 = get(paste0('vs_',mod1))
vs_2 = get(paste0('vs_',mod2))

perm_test_dt = merge(vs_1, vs_2, by=c('year', 'month'))
setnames(perm_test_dt,c('vs.x', 'vs.y'), c(mod1, mod2))

# permutation test for mod1 ~ mod2
pt_vs = permutation_test_difference(perm_test_dt[,get(mod1)],perm_test_dt[,get(mod2)], N=N)


#x_lim for the plot:
x_lim_max = 1.1*max(abs(c(pt_vs$d_bar,pt_vs$D)))

#make ticks at:
ticks_pos = c(-0.6*x_lim_max, 0, 0.6*x_lim_max)
ticks_vals = round(ticks_pos,4)


hist(pt_vs$D, xlim = c(-x_lim_max,x_lim_max),xlab = '',
     main = latex2exp::TeX(paste0(mod_names[1],' vs. ',mod_names[2])),breaks = 20,
     col = rgb(t(col2rgb('gray')/255),alpha = 0.5,maxColorValue = 1),
     border = rgb(t(col2rgb('gray')/255),alpha = 0.5,maxColorValue = 1),
     probability = TRUE,     
     axes = FALSE)


# axis

# prefer scientific notation 
options(scipen = -2)

axis(side = 1, cex.axis = 1, at = ticks_pos, labels = ticks_vals)


#actual value

abline(v = pt_vs$d_bar,col = 'red')

#quantiles

qq_1 = quantile(pt_vs$D,0.05)

abline(v = qq_1,lty = 2)

# qq_2 = quantile(pt_vs$D,0.01)
# 
# abline(v = qq_2,lty = 3)


}

dev.off()


##############################################

### multivariate rank histograms for route ###

# range of probability the probability plots
rr = c(0,100)

### average rank histograms ###

pdf(paste0(plot_dir,'avg_rhs_route.pdf'))

par(oma = c(1,1,1,1), mfrow=c(2,2), mar=c(2,1,2,1) )
par('cex' = 1.4,'cex.axis' = 0.75)

rhist_dt(rks_PCA_mc_route[,.(YM,av.rk.obs)], 
         max_rk = fc_ens_size +1,
         breaks = breaks, 
         hist_xlab = "average",
         hist_ylim = rr,
         mn = latex2exp::TeX('average rank, $\\widehat{\\Sigma}^{mc}$')
)

rhist_dt(rks_PCA_ac_route[,.(YM,av.rk.obs)], 
         max_rk = fc_ens_size +1,
         breaks = breaks, 
         hist_xlab = "average",
         hist_ylim = rr,
         mn = latex2exp::TeX('average rank, $\\widehat{\\Sigma}^{ac}$')
)

rhist_dt(rks_GS_route[,.(YM,av.rk.obs)], 
         max_rk = fc_ens_size +1,
         breaks = breaks, 
         hist_xlab = "average",
         hist_ylim = rr,
         mn = latex2exp::TeX('average rank, GS')
)

rhist_dt(rks_ECC_route[,.(YM,av.rk.obs)], 
         max_rk = ens_size +1,
         breaks = breaks, 
         hist_xlab = "average",
         hist_ylim = rr,
         mn = latex2exp::TeX('average rank, ECC')
)



dev.off()

### band depth rank ###

### average rank histograms ###

pdf(paste0(plot_dir,'bd_rhs_route.pdf'))

par(oma = c(1,1,1,1), mfrow=c(2,2), mar=c(2,1,2,1) )
par('cex' = 1.4,'cex.axis' = 0.75)

rhist_dt(rks_PCA_mc_route[,.(YM,bd.rk.obs)], 
         max_rk = fc_ens_size +1,
         breaks = breaks, 
         hist_xlab = "average",
         hist_ylim = rr,
         mn = latex2exp::TeX('band depth rank, $\\widehat{\\Sigma}^{mc}$')
)

rhist_dt(rks_PCA_ac_route[,.(YM,bd.rk.obs)], 
         max_rk = fc_ens_size +1,
         breaks = breaks, 
         hist_xlab = "average",
         hist_ylim = rr,
         mn = latex2exp::TeX('band depth rank, $\\widehat{\\Sigma}^{ac}$')
)

rhist_dt(rks_GS_route[,.(YM,bd.rk.obs)], 
         max_rk = fc_ens_size +1,
         breaks = breaks, 
         hist_xlab = "average",
         hist_ylim = rr,
         mn = latex2exp::TeX('band depth rank, GS')
)

rhist_dt(rks_ECC_route[,.(YM,bd.rk.obs)], 
         max_rk = ens_size +1,
         breaks = breaks, 
         hist_xlab = "average",
         hist_ylim = rr,
         mn = latex2exp::TeX('band depth rank, ECC')
)



dev.off()


