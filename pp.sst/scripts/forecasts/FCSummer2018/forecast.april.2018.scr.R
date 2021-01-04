
###########################################################################

############## script for issuing the summer forecast 2018 ################

###########################################################################


#Issues the forecast for SST for summer 2018 in the northern atlantic ocean.

# This script is based on the master.NAO script with some stuff left out and some minor modifications.
# The model is trained on April vintage instead of the most recent vintage and assumes a gap in the data of 7 years 
# prior to the forecast. PCA uses all data available until 2010.


rm(list = ls())

library(PostProcessing)
library(data.table)
library(parallel)
library(fExtremes)
library(fields)
library(vegan)

setwd("~/NR/SFE")
options(max.print = 1e3)

# box for analysis, months for the forecast, and name abbreviation for saved files and directories

lat_box = c(40,85)
lon_box = c(-30,60)

months = 4:9

name_abbr = "Apr18" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr)
dir.create(save_dir, showWarnings = FALSE)

plot_dir = paste0("./figures/", name_abbr)
dir.create(plot_dir, showWarnings = FALSE)


########################################
### construct and load wide data set ###
########################################



y_start = 1986 # starting at 1985 gives an issue with choosing the April vintage
y_stop = 2010
vintage = "Apr"
data_dir = "~/PostClimDataNoBackup/SFE/FcApr2018/"
grid_mapping_loc = "~/PostClimDataNoBackup/SFE/Derived/"
output_loc = save_dir 
output_name = paste0("dt_combine_",name_abbr,"_wide.RData")


# --- here comes a variation on the function make_combined_wide() ---
  
  ##----- Load Grid Mapping ---
  ff = paste0(grid_mapping_loc,"dt_map.RData")
  if(file.exists(ff))
  {
    load(ff)
    names(dt_map)= c("Lon_Obs","Lat_Obs","Lon_Ens","Lat_Ens") ## Get rid of this eventually.
  }else{
    stop("Could not find grid mapping info")
  }
  ##--------------------------
  
  ##------ Loop ----------
  dt_combine_all = list()
  k = 1
  for(y in y_start:y_stop)
  {
    for(m in months)
    {
      print(c(y,m))
      dt_ens = load_ensemble(y,m,vintage,data_dir = data_dir)
      dt_obs = load_observations(y,m)
      dt_combine_all[[k]] = combine_data_wide(dt_ens, dt_obs, dt_map)
      dt_combine_all[[k]][,year:=y]
      dt_combine_all[[k]][,month:=m]
      k = k + 1
    }
  }
  
  # attach 2018 forecast
  
  y = 2018
  for(m in months)
  {
    print(c(y,m))
    dt_ens = load_ensemble(y,m,vintage,data_dir = data_dir)
    dt_obs = load_observations(2010,1) # hacked - we don't have an observation for 2018, so we just attach any observation
    dt_combine_all[[k]] = combine_data_wide(dt_ens, dt_obs, dt_map)
    dt_combine_all[[k]][,year:=y]
    dt_combine_all[[k]][,month:=m]
    k = k + 1
  }
  
  ##------------------------
  
  ##--------- Combine -----
  DT = rbindlist(dt_combine_all)
  DT[, YM := year * 12 + month]
  setkey(DT, "YM", "Lon", "Lat")
  ##------------------------
  
  ##---- Restrict ---
  DT = DT[ (Lon >= lon_box[1]) & (Lon <= lon_box[2] )& (Lat >= lat_box[1]) & (Lat <= lat_box[2])]
  
  ##----- Should I save or should I go? -----
  
    if(is.null(output_name))
    {
      output_name = paste0("dt_combine_",vintage,"_wide.RData")
    }
    f_name = paste0(output_loc,"/",output_name)
    save(DT,file = f_name)
    

  ##-------------------------------------------
  
  

 DT = load_combined_wide(data_dir = save_dir, output_name = paste0("dt_combine_",name_abbr,"_wide_bc.RData"))



########################################
####### univariate modelling ###########
########################################


############### bias ###################


### run bias analysis - train moving averages subject on the 7 most recent years of data missing  ###

validation_years = 2001:2010

# for simple moving averages

num.years = DT[,range(year)][2] - DT[,range(year)][1] + 1

sc_sma = list()
dummy_function = function(k){
  temp = bias_correct(dt = DT, 
                      method = "sma", 
                      par_1 = k,
                      scores = TRUE,
                      eval_years = validation_years,
                      saveorgo = FALSE,
                      skip = 7)
  sc_sma[[k]] = temp[,"win_length" := k]
}

sc_sma = mclapply(X = 1:(num.years-1), FUN = dummy_function, mc.cores = 8)
sc_sma = rbindlist(sc_sma)

save(sc_sma, file = paste0(save_dir,"/scores.bc.sma.Rdata"))

# for exponential moving averages

par.vec = seq(0.05,0.4,length.out = 24)

sc_ema = list()

dummy_function = function(k){
  temp = bias_correct(dt = DT, 
                      method = "ema",
                      par_1 = par.vec[k], 
                      scores = TRUE,
                      eval_years = validation_years,
                      saveorgo = FALSE,
                      skip = 7)
  sc_ema[[k]] = temp[,"a" := par.vec[k]]
}

sc_ema = mclapply(X = 1:length(par.vec), FUN = dummy_function,mc.cores = 8)
sc_ema = rbindlist(sc_ema)

save(sc_ema, file = paste0(save_dir,"/scores.bc.ema.Rdata"))


#### plotting scores for different ways of bias correction ####

load(paste0(save_dir,"/scores.bc.sma.Rdata"))
load(paste0(save_dir,"/scores.bc.ema.Rdata"))

# ensure that they are plotted on the the same range
y_range = range(c(sc_sma[,sqrt(MSE)],sc_ema[,sqrt(MSE)]))  

## plot for sma ##

pdf(paste0(plot_dir,"/mean_scores_sma.pdf"))
plot(x = sc_sma[,win_length],
     y = sc_sma[,sqrt(MSE)],
     ylim = y_range,
     type = "b",
     col = "blue",
     main = paste0("RMSE for ",name_abbr," bias correction by SMA"),
     xlab = "window length",
     ylab = "RMSE"
)

# highlight minimum and add minimum reference line 
abline(h = sc_sma[,min(sqrt(MSE))], lty = "dashed", col = adjustcolor("blue",alpha = .5))

min_loc_RMSE = sc_sma[,which.min(MSE)]

points(x = sc_sma[,win_length][min_loc_RMSE],
       y = sc_sma[,sqrt(MSE)][min_loc_RMSE],
       col = "blue",
       bg = "blue",
       pch = 21)

dev.off()

## plot for ema ##

pdf(paste0(plot_dir,"/mean_scores_ema.pdf"))
plot(x = sc_ema[,a],
     y = sc_ema[,sqrt(MSE)],
     ylim = y_range,
     type = "b",
     col = "blue",
     main = paste0("RMSE for ",name_abbr," bias correction by EMA"),
     xlab = "weight parameter a",
     ylab = "RMSE"
)

# highlight minimum and add minimum reference line 
abline(h = sc_ema[,min(sqrt(MSE))], lty = "dashed", col = adjustcolor("blue",alpha = .5))

min_loc_RMSE = sc_ema[,which.min(MSE)]

points(x = sc_ema[,a][min_loc_RMSE],
       y = sc_ema[,sqrt(MSE)][min_loc_RMSE],
       col = "blue",
       bg = "blue",
       pch = 21)

dev.off()

### finding optimal way of bias correcion ###

if(sc_sma[,min(MSE)] < sc_ema[,min(MSE)]){
  print(paste0("optimal bias correction uses simple moving averages with window length of ",sc_sma[,which.min(MSE)], " years, and achieves a RMSE of ",round(sc_sma[,sqrt(min(MSE))],3),
               ". Best correction with exponential moving averages achieves a RMSE of ",round(sc_ema[,sqrt(min(MSE))],3),"."))
  opt_par = c("sma",sc_sma[,which.min(MSE)])
} else{
  print(paste0("optimal bias correction uses exponential moving averages with parameter a = ",round(sc_ema[,a][sc_sma[,which.min(MSE)]],3),
               ", and achieves a RMSE of ",round(sc_ema[,sqrt(min(MSE))],3),". Best correction with simple moving averages achieves a RMSE of ",round(sc_sma[,sqrt(min(MSE))],3),"."))
  opt_par = c("ema",sc_ema[,a][sc_sma[,which.min(MSE)]])
}


### bias correction ###

bias_correct(dt = DT,
             method = opt_par[1],
             par_1 = as.double(opt_par[2]),
             save_dir = paste0(save_dir,"/"),
             file_name = paste0("dt_combine_",name_abbr,"_wide_bc.RData"),
             skip = 7)

DT = load_combined_wide(data_dir = save_dir, output_name = paste0("dt_combine_",name_abbr,"_wide_bc.RData"))




#################################################
########## multivariate modelling ###############
#################################################


###################  PCA  #######################

training_years = 1986:2010
eval_years = 2018

ens_size = 9

##### compute and save principal components #####

cov_dir = paste0(save_dir,"/PCACov")
dir.create(cov_dir, showWarnings = FALSE)

for_res_cov(Y = training_years,dt = DT, save_dir = cov_dir,ens_size = ens_size)

############## set up PCA #######################

setup_PCA(dt = DT, y = eval_years,m = months, cov_dir = cov_dir)


###### Example plots of forecasted SST and anomalies w.r.t climatology #######

opt_num_PCs = 17 # This is an educated guess: for 16 years training the results are best for 11 PCs and very good for 9 to 17 PCs, now we have 26 years of training.

ex_depth = opt_num_PCs
ex_months = months
ex_month_names = c("April","May","June","July","August","September")
ex_year = 2018

clim_years = 1986:2010 # the years to compute the climatology from

MC_sample_size = 30   # number of plots with independently generated noise

PCs = opt_num_PCs      # number of considered principal components


#compute climatology
climatology = DT[year %in% clim_years, clim := mean(SST_bar),by = .(grid_id,month)][year == min(year) ,.(Lon,Lat,month,clim)]
#print climatology every year
DT[year == 2018, clim := DT[year == min(clim_years),clim]]



for(m in ex_months){
  print(paste0("month = ",m))
  
  #generate noise:
  no.dt = list()
  for(i in 1:MC_sample_size){
    no.dt[[i]] = forecast_PCA(m = m, y = ex_year, PCA_depth = PCs, saveorgo = FALSE)[,.(Lon,Lat,noise), keyby = .(Lon,Lat)][,noise]
  }
  no.dt = as.data.table(no.dt)
  setnames(no.dt,paste0("no",1:MC_sample_size))
  
  DT_pca_plot = no.dt[,c("year","month","YM","Lat","Lon","SST_bar","Bias_Est",paste0("Ens",1:ens_size)) := 
                        DT[year == ex_year & month == m, c("year","month","YM","Lat","Lon","SST_bar","Bias_Est",paste0("Ens",1:ens_size)),
                           with = FALSE]]
  DT_pca_plot[,clim := climatology[month == m, clim]]
  
  # choose random ensemble members (REM) and generate forecast as REM + bias + noise
  
  ens_mem = sample.int(ens_size,MC_sample_size,replace = TRUE)
  for(i in 1:MC_sample_size){
    dummy_dt = DT_pca_plot[,.SD,.SDcols = c(paste0("no",i),paste0("Ens",ens_mem[i]),"Bias_Est")]
    forecast = dummy_dt[[1]] + dummy_dt[[2]] + dummy_dt[[3]]
    DT_pca_plot = DT_pca_plot[,paste0("fc",i) := forecast]
  }
  
  #forecast plots:

  trc = function (x){ 
    truncation.value = -1.769995
    x = truncation.value * (x < truncation.value) + x * (x >= truncation.value)
    return(x)}
  
  for(i in 1:MC_sample_size){
    plot_diagnostic(DT_pca_plot[,.(Lon,Lat,trc(eval(parse(text = paste0("fc",i)))))],
                    #rr = rr_sst,
                    mn = paste0("SST forecast for ",ex_month_names[which(ex_months == m)]),
                    save_pdf = TRUE, 
                    save_dir = paste0(plot_dir,"/"),
                    file_name = paste0("m",m,"_fc",i),
                    stretch_par = .8)
  }
  
  ### anomaly plots ###
  
  for(i in 1:MC_sample_size){
    plot_diagnostic(DT_pca_plot[,.(Lon,Lat,trc(eval(parse(text = paste0("fc",i))))-clim)],
                    rr = c(-3,3),
                    mn = paste0("Anomaly forecast for ",ex_month_names[which(ex_months == m)]),
                    save_pdf = TRUE, 
                    save_dir = paste0(plot_dir,"/"),
                    set_white = 0,
                    file_name = paste0("m",m,"_afc",i),
                    stretch_par = .8)
  }
  
}

################## plots that Erik actually wants ######################

month_vec = c("Apr.","May","June","July","Aug.","Sept.")

# mean SST forecasts

for (m in months){
  DT_plot = DT[year == 2018 & month == m ,.(Lon,Lat,trc(SST_hat))]
  plot_diagnostic(DT_plot,
                  mn = paste0("SST forecast for ",month_vec[m-3]," 2018"),
                  save_pdf = TRUE,
                  save_dir = paste0(plot_dir,"/"),
                  file_name = paste0("mean_SST_forecast_m",m),
                  stretch_par = .8
                  )
}

# mean anomaly forecasts

for (m in months){
  DT_plot = DT[year == 2018 & month == m ,.(Lon,Lat,trc(SST_hat)-clim)]
  plot_diagnostic(DT_plot,
                  mn = paste0("SST anomaly fc for ",month_vec[m-3]," 2018"),
                  rr = c(-3,3),
                  save_pdf = TRUE,
                  save_dir = paste0(plot_dir,"/"),
                  file_name = paste0("mean_ano_forecast_m",m),
                  set_white = 0,
                  stretch_par = .8
  )
}


###### getting terciles of the empirical anomaly distribution ######

q1 = 1/3
q2 = 2/3

# plot the percentage of ensemble members below q1 and above q2

n=100

DT_new = forecast_PCA_new(dt = DT,
                          y = 2018,
                          m = months,
                          n=n,
                          PCA_depth = opt_num_PCs,
                          ens_member = FALSE,
                          saveorgo = FALSE,
                          cov_dir = cov_dir)

rel_col_names = c("Lon","Lat","month",paste0("fc",1:n,"PC",opt_num_PCs)) 
DT_new = DT_new[,.SD,.SDcols = rel_col_names]

# get observations for climatology

obs_years = 1986:2010

prep_DT = DT[year < 2017 & month > 3,.(Lon,Lat,year,month,SST_bar)][,"name" := paste0("SST",year)]
prep_DT[,SST_SD := sd(SST_bar),by = .(month,Lon,Lat)]

SST_obs = dcast(prep_DT, Lon + Lat + month + SST_SD ~ name,value.var = "SST_bar")

rowrank = function(A){# takes m x n matrix and returns vector of length m containing the ranks of the entries of the first column 
  dA = dim(A)
  a = c()
  for( i  in 1:dA[1]){
    a = c(a,rank(A[i,],na.last = "keep",ties.method = "random")[1])
  }
  return(a)
  }

setkey(SST_obs,month,Lon,Lat)
setkey(DT_new,month,Lon,Lat)

res_dt = merge(SST_obs,DT_new,all = TRUE)

for(i in 1:n){
  print(i)
  col_names = c(paste0("fc",i,"PC",opt_num_PCs),paste0("SST",obs_years))
  vec = c()
  for(m in months){
    vec = c(vec,rowrank(res_dt[month == m,.SD,.SDcols = col_names]))
  }
  res_dt[,paste0("rk",i):= vec]
}



# get the percentages
threshold1 = q1 *  26 # 25 years of obs and 1 forecast
threshold2 = q2 *  26 #

results = res_dt[,sq1 := sum(.SD < threshold1),.SDcols = paste0("fc",1:n,"PC",opt_num_PCs),by = .(month,Lon,Lat)]
results = res_dt[,sq2 := sum(.SD > threshold2),.SDcols = paste0("fc",1:n,"PC",opt_num_PCs),by = .(month,Lon,Lat)]



below_t1 = res_dt[,rk1<threshold1]
for(i in 2:n){below_t1 = below_t1 + res_dt[,eval(parse(text = paste0("rk",i)))<threshold1]
}

above_t2 = res_dt[,rk1>threshold2]
for(i in 2:n){above_t2 = above_t2 + res_dt[,eval(parse(text = paste0("rk",i)))>threshold2]
}

results[,below_t1:= below_t1][,above_t2:= above_t2]

results[SST_SD < 0.2 & !is.na(below_t1),below_t1 := (threshold1+threshold2)/2]
results[SST_SD < 0.2 & !is.na(above_t2),above_t2 := (threshold1+threshold2)/2]



for(m in months){
  plot_diagnostic(results[month == m,.(Lon,Lat,below_t1)],
                mn = paste0("percent below first tercile, ",month_vec[m-3]," 2018"),
                save_pdf = TRUE,
                col_scheme = "wb",
                save_dir = paste0(plot_dir,"/"),
                file_name = paste0("below_first_tercile_m",m),
                stretch_par = .8)
  
  plot_diagnostic(results[month == m,.(Lon,Lat,above_t2)],
                  mn = paste0("percent above second tercile, ",month_vec[m-3]," 2018"),
                  save_pdf = TRUE,
                  save_dir = paste0(plot_dir,"/"),
                  col_scheme = "wr",
                  file_name = paste0("above_second_tercile_m",m),
                  stretch_par = .8)
  
}


 
