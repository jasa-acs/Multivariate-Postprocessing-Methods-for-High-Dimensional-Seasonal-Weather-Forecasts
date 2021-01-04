

##### setting up ######

rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "Full" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))


ex_plot_dir = paste0(plot_dir,"Examples/")
dir.create(ex_plot_dir, showWarnings = FALSE)


### example plots ###

ex_month = 9

cex = 1.1

# get climatology

training_years = 1985:2000

DT[year %in% training_years,clim := mean(SST_bar),by =.(month,grid_id)]
clim_vec = rep(DT[year == min(training_years),clim],length(DT[,unique(year)]))
DT[,clim := clim_vec]
DT[,ano := SST_bar - clim]



######################

rr = c(-5,5)

ex_years = 2001:2003

mn = ""
add_title = TRUE # do we want titles on the plots or do we want to play the guess-my-method game?
repetitions = 1



# SST residuals

r_plus = max(abs(range(DT[year %in% ex_years & month %in% ex_month,SST_bar - SST_hat],na.rm = TRUE)))

rr = c(-r_plus,r_plus)

for(m in ex_month)
{
  for(y in ex_years)
  {
    if(add_title)
    {
      mn = paste0("Forecast residual, ",m," / ",y)
    }
    for(k in 1:repetitions)
    {
      plot_diagnostic(DT[year == y & month == m,.(Lon,Lat,SST_bar - SST_hat)],
                      mn = mn,
                      ylab = "", xlab = "",
                      rr = rr,
                      cex = cex,
                      save_pdf = TRUE,
                      save_dir = ex_plot_dir,
                      file_name = paste0("res",k,"_m",m,"_y",y))  
    }
  }
}

# raw ensemble member and forecast

for(m in ex_month)
{
  for(y in ex_years)
  {
    ens = sample(ens_size,repetitions)
    
    rr = range(DT[year == y & month == m,.(SST_bar,Ens_bar)],na.rm = TRUE)
      
      if(add_title)
      {
        mn = paste0("raw forecast, ",m," / ",y)
      }
      plot_diagnostic(DT[year == y & month == m,.(Lon,Lat,Ens_bar)],
                      mn = mn,
                      rr = rr,
                      cex = cex,
                      ylab = "", xlab = "",
                      save_pdf = TRUE,
                      save_dir = ex_plot_dir,
                      file_name = paste0("rfc",k,"_m",m,"_y",y))  
      
      if(add_title)
      {
        mn = paste0("observed SST, ",m," / ",y)
      }
      plot_diagnostic(DT[year == y & month == m,.(Lon,Lat,SST_bar)],
                      mn = mn,
                      rr = rr,
                      ylab = "", xlab = "",
                      save_pdf = TRUE,
                      save_dir = ex_plot_dir,
                      file_name = paste0("SST",k,"_m",m,"_y",y))  
      
  }
}

# show principal components:
ex_years = 2002

PCs = 2

prin_comp_dt = get_PCs(dt = DT, y = ex_years, m = ex_month, PCA_depth = PCs, cov_dir = PCA_dir)

#without marginal correction

for(y in ex_years){
  rr_PC_raw = range(prin_comp_dt[year == y,.SD,.SDcols = c(paste0("PC",1:PCs))],na.rm = TRUE)
  rr_PC_raw = c(-max(abs(rr_PC_raw)), max(abs(rr_PC_raw)))
  
  for(d in 1:PCs)
  {
    plot_diagnostic(prin_comp_dt[year == y,.SD,.SDcols = c("Lon","Lat",paste0("PC",d))],
                    rr = rr_PC_raw,
                    mn = paste0("PC ",d,", nmc, 0",ex_month," / ",y),
                    xlab = "", ylab = "",
                    cex = cex,
                    save_pdf = TRUE,
                    save_dir = ex_plot_dir,
                    file_name = paste0("PC",d,"_y",y,"_raw"))
    
  }
}

#with marginal correction:

for(y in ex_years){
  rr_PC_raw = range(prin_comp_dt[year == y,.SD,.SDcols = c(paste0("PC_marcor_",1:PCs))],na.rm = TRUE)
  rr_PC_raw = c(-max(abs(rr_PC_raw)), max(abs(rr_PC_raw)))
  
  for(d in 1:PCs)
  {
    plot_diagnostic(prin_comp_dt[year == y,.SD,.SDcols = c("Lon","Lat",paste0("PC_marcor_",d))],
                    rr = rr_PC_raw,
                    mn = paste0("PC ",d,", 0",ex_month," / ",y),
                    save_pdf = TRUE,
                    xlab = "", ylab = "",
                    cex = cex,
                    save_dir = ex_plot_dir,
                    file_name = paste0("PC",d,"_y",y,"_mc"))
    
  }
}




########## multivariate methods #################

name_abbr = "Pres_Bergen" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))


ex_plot_dir = paste0(plot_dir,"Examples/")
dir.create(ex_plot_dir, showWarnings = FALSE)

ex_month = 9

# get climatology

training_years = 1985:2000

DT[year %in% training_years,clim := mean(SST_bar),by =.(month,grid_id)]
clim_vec = rep(DT[year == min(training_years),clim],length(DT[,unique(year)]))
DT[,clim := clim_vec]
DT[,ano := SST_bar - clim]



######################

rr = c(-5,5)


ex_years = 2001:2003
repetitions = 1
mn = ""
add_title = TRUE # do we want titles on the plots or do we want to play the guess-my-method game?

# observed SST

for(m in ex_month)
{
  for(y in ex_years)
  {
    if(add_title)
    {
      mn = paste0("SST anomaly ",m," / ",y)
    }
    for(k in 1:repetitions)
    {
      plot_diagnostic(DT[year == y & month == m,.(Lon,Lat,ano)],
                      mn = mn,
                      rr = rr,
                      ylab = "", xlab = "",
                      save_pdf = TRUE,
                      save_dir = ex_plot_dir,
                      file_name = paste0("Ano",k,"_m",m,"_y",y))  
    }
  }
}




ex_PCs = 50
ex_years = 2001:2003
repetitions = 1

DT_PCA = forecast_PCA_new(dt = DT,y = ex_years, m = ex_month, n = repetitions, PCA_depth = ex_PCs,saveorgo = FALSE, cov_dir = PCA_dir)
setkeyv(DT_PCA,key(DT))
DT_PCA[,clim:= DT[month %in% ex_month & year %in% ex_years,clim]]

# PCA corrected forecast

add_title = TRUE
mn = ""

for(PC in ex_PCs)
{
  for(m in ex_month)
  {
    for(y in ex_years)
    {
      for(k in 1:repetitions)
      {
        if(add_title)
        {
          mn = paste0("PCA forecast ",m," / ",y)
        }
        plot_diagnostic(DT_PCA[year == y & month == m, .(Lon,Lat,eval(parse(text = paste0("fc",k,"PC",PC)))-clim)],
                        mn = mn,
                        save_pdf = TRUE,
                        rr=rr,
                        ylab = "", xlab = "",
                        save_dir = ex_plot_dir,
                        file_name = paste0("PCA_fc",k,"_m",m,"_y",y))
        
      }
    }
  }
}

ex_years = 2001:2003
# geostationary:
DT_geostat = forecast_geostat(dt = DT, n = repetitions, y = ex_years, m = ex_month,
                              ens_size = ens_size, saveorgo = FALSE, data_dir = geostat_dir)
setkey(DT_geostat,year,Lon,Lat)
DT_geostat[,clim := DT[year %in% ex_years & month %in% ex_month,clim]]

add_title = TRUE
mn = ""

for(y in ex_years)
{
  for(m in ex_month)
  {
    if(add_title)
    {
      mn = paste0("geostationary forecast ",m," / ",y)
    }
    for(k in 1:repetitions)
    {
      plot_diagnostic(DT_geostat[year == y & month == m,.(Lon,Lat,.SD-clim),
                                 .SDcols = paste0("fc",k)],
                      mn = mn,
                      rr=rr,
                      ylab = "", xlab = "",
                      save_pdf = TRUE,
                      save_dir = ex_plot_dir,
                      file_name = paste0("fc_geostat",k,"_m",m,"_y",y))
    }
  }
}

# ECC

ex_years = 2001:2003

DT_ECC = forecast_ECC(dt = DT[year %in% ex_years & month %in% ex_month,],saveorgo = FALSE)
DT_ECC[,clim := DT[year %in% ex_years & month %in% ex_month,clim]]

add_title = TRUE
mn = ""

for(y in ex_years)
{
  for(m in ex_month)
  {
    if(add_title)
    {
      mn = paste0("ECC forecast ",m," / ",y)
    }
    ens = sample(ens_size,repetitions)
    for(k in 1:repetitions)
    {
      plot_diagnostic(DT_ECC[year == y,.(Lon,Lat,.SD - clim),
                             .SDcols = paste0("ecc_fc",ens[k])],
                      mn = mn,
                      rr = rr,
                      save_pdf = TRUE,
                      save_dir = ex_plot_dir,
                      file_name = paste0("fc_ecc",k,"_m",m,"_y",y))
    }
  }
}



#####################




########### marginal SDs and Ensemble spread #############

# Ensemble spread

for(m in ex_month)
{
  for(y in ex_years)
  {
    
    for(k in 1:repetitions)
    {
      
      DT_temp = DT[year %in% training_years, Ens_spread := sqrt(mean(Ens_sd^2)), by = .(Lon,Lat,month)]
      rr_sd = range(c(0,DT_temp[,Ens_spread],DT[year == y & month == m,SD_hat]),na.rm = TRUE)
      
      if(add_title)
      {
        mn = paste0("July ensemble spread")
      }
      plot_diagnostic(DT_temp[year == min(year) & month == m,
                              .(Lon,Lat,Ens_spread)],
                      mn = mn,
                      rr = rr_sd,
                      save_pdf = TRUE,
                      col_scheme = "wr",
                      save_dir = ex_plot_dir,
                      file_name = paste0("ens_spread",k,"_m",m,"_y",y))
      
      if(add_title)
      {
        mn = paste0("July SD, estimated from data")
      }
      plot_diagnostic(DT[year == y & month == m,
                         .(Lon,Lat,SD_hat)],
                      mn = mn,
                      rr = rr_sd,
                      save_pdf = TRUE,
                      col_scheme = "wr",
                      save_dir = ex_plot_dir,
                      file_name = paste0("sd_hat",k,"_m",m,"_y",y))
      
    }
  }
}




# plotting marginal SDs

ex_years_sd = 2002:2007
prin_comp_dt = get_PCs(dt = DT, y = ex_years_sd, m = ex_month, PCA_depth = PCs, cov_dir = PCA_dir)

for(y in ex_years_sd){
  rr = range(prin_comp_dt[,SD_hat],na.rm = TRUE)
  rr = c(0,max(rr))
  plot_diagnostic(prin_comp_dt[year == y,.(Lon,Lat,SD_hat)],
                  mn = paste0("SD 0",ex_month," / ",y),
                  rr = rr,
                  col_scheme = "wr",
                  save_pdf = TRUE,
                  save_dir = paste0(plot_dir,"/"),
                  file_name = paste0("SD_hat",y))
  
}

# plotting difference of SDs

ex_years_sd = 2002:2007
prin_comp_dt = get_PCs(dt = DT, y = ex_years_sd, m = ex_month, PCA_depth = PCs, cov_dir = PCA_dir)

y = min(ex_years_sd)

rr = range(prin_comp_dt[,SD_hat],na.rm = TRUE)
rr = c(0,max(rr))
plot_diagnostic(prin_comp_dt[year == y,.(Lon,Lat,SD_hat)],
                mn = paste0("SD 0",ex_month," / ",y),
                rr = rr,
                col_scheme = "wr",
                save_pdf = TRUE,
                save_dir = paste0(plot_dir,"/"),
                file_name = paste0("SD_hat_diff",y))


for(y in ex_years_sd[2:length(ex_years_sd)]){
  rr = range(prin_comp_dt[,SD_hat],na.rm = TRUE)
  rr = c(0,max(rr))
  plot_diagnostic(prin_comp_dt[year == y,.(Lon,Lat,SD_hat)],
                  mn = paste0("SD 0",ex_month," / ",y),
                  rr = rr,
                  col_scheme = "wr",
                  save_pdf = TRUE,
                  save_dir = paste0(plot_dir,"/"),
                  file_name = paste0("SD_hat",y))
  
}







############################################################


###### Example plots of forecasted SST and anomalies w.r.t climatology #######

ex_depth = opt_num_PCs
ex_months = 4:9
ex_month_names = c("April","May","June","July","August","September")
ex_year = 2010

clim_years = 1985:2009 # the years to compute the climatology from

MC_sample_size = 10   # number of plots with independently generated noise

ens_size = 9   #size of forecast ensemble

PCs = opt_num_PCs      # number of considered principal components


#compute climatology
climatology = DT[year %in% clim_years, clim := mean(SST_bar),by = .(grid_id,month)][year == min(year) ,.(Lon,Lat,month,clim)]



for(m in ex_months){
  print(paste0("month = ",m))
    
    #generate noise:
    no_dt = list()
    for(i in 1:MC_sample_size){
      no_dt[[i]] = forecast_PCA(m = m, y = ex_year, PCA_depth = PCs, saveorgo = FALSE)[,.(Lon,Lat,noise), keyby = .(Lon,Lat)][,noise]
    }
    no_dt = as.data.table(no_dt)
    setnames(no_dt,paste0("no",1:MC_sample_size))
    
    DT_pca_plot = no_dt[,c("year","month","YM","Lat","Lon","SST_bar","Bias_Est",paste0("Ens",1:ens_size)) := 
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
    
    rr_sst = range(na.omit(DT_pca_plot[,.SD,.SDcols = c(paste0("fc",1:MC_sample_size))]))
    
    for(i in 1:MC_sample_size){
      plot_diagnostic(DT_pca_plot[,.(Lon,Lat,eval(parse(text = paste0("fc",i))))],
                      rr = rr_sst,
                      mn = paste0("SST forecast for ",ex_month_names[which(ex_months == m)]),
                      save_pdf = TRUE, 
                      save_dir = paste0(plot_dir,"/"),
                      file_name = paste0("m",m,"_fc",i),
                      stretch_par = .8)
    }
    
    #anomaly plot 
    rr_clim = range(na.omit(DT_pca_plot[,.SD - clim,.SDcols = c(paste0("fc",1:MC_sample_size))]))
    
    for(i in 1:MC_sample_size){
      plot_diagnostic(DT_pca_plot[,.(Lon,Lat,eval(parse(text = paste0("fc",i)))-clim)],
                      rr = rr_clim,
                      mn = paste0("Anomaly forecast for ",ex_month_names[which(ex_months == m)]),
                      save_pdf = TRUE, 
                      save_dir = paste0(plot_dir,"/"),
                      file_name = paste0("m",m,"_afc",i),
                      stretch_par = .8)
    }
    
  }



