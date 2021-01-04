rm(list=ls())

library(fExtremes)
library(fields)
library(vegan)
library(PostProcessing)

setwd("~/NR/SFE")
options(max.print = 1e3)


#########################################################################
  
  
  ##### rank histograms for raw ensemble forecast
  
  # ---- create data table containing the (multivariate) ranks of the observation in the raw forecast ensemble ---
  
  Lon_min = -60
  Lon_max = 15
  Lat_min = 30
  Lat_max = 70
  
  DT = load_combined_wide(bias = TRUE)
  DT = DT[Lon >= Lon_min & Lon <= Lon_max & Lat >= Lat_min & Lat <= Lat_max]
  
  y = 2001:2010  #validation years
  m = 1:12       #validation months
  ens.size = 9
  
  DT_raw = DT[year %in% y & month %in% m,]
  
  ym = unique(DT_raw[,YM])
  
  ranks.matrix = matrix(ym,nrow = length(ym),ncol = 1+4*(ens.size +1)) 
        #ncol: 1 col for YM, 4 methods of ranking, for each we get ranks for observation and each ensemble member 
  
  drm = dim(ranks.matrix)
  
  YM.ind = 0
  
  for(yearmonth in ym){
    print(paste0("YM = ",yearmonth,"/",ym[length(ym)]))
    YM.ind = YM.ind + 1
    fc_obs_mat = na.omit(DT_raw[YM == yearmonth,.SD, .SDcols = c("SST_bar",paste0("Ens",1:9))])
    
    # get ranks
    ranks.matrix[YM.ind,2:drm[2]] = c(mst.rank(as.matrix(fc_obs_mat)),
                                      mv.rank(as.matrix(fc_obs_mat)), 
                                      avg.rank(as.matrix(fc_obs_mat)),
                                      bd.rank(as.matrix(fc_obs_mat)))
  }
  
  names.vec = c("YM","mst.rk.obs",paste0("mst.r.",1:ens.size),
                "mv.rk.obs",paste0("mv.r.",1:ens.size),
                "av.rk.obs",paste0("av.r.",1:ens.size),
                "bd.rk.obs",paste0("bd.rk",1:ens.size))
  
  ranks = data.table(ranks.matrix)
  
  setnames(ranks, names.vec)
  
  # --- save ---
  
  save.dir = "~/PostClimDataNoBackup/SFE/Derived/Calibration/"
  
  save(ranks,file = paste0(save.dir,"ranks_raw.Rdata"))
  
  # ---- plotting ----
  pdf(file="./figures/rks_raw.pdf",width=8,height=2,points=12)
  par(mfrow=c(1,4),mex=0.5,mar=c(2.5,2.5,2.5,2.5)+0.1,mgp=c(0.5,0,0))
  rhist.dt(ranks[,.(YM,mst.rk.obs)],hist_xlab = "minimum spanning tree")
  rhist.dt(ranks[,.(YM,mv.rk.obs)],hist_xlab = "multivariate")
  rhist.dt(ranks[,.(YM,av.rk.obs)],hist_xlab = "average")
  rhist.dt(ranks[,.(YM,bd.rk.obs)],hist_xlab = "band depth")
  dev.off()
  
  ###################
  
  
  ##### rank histograms for bias corrected ensemble forecast
  
  # ---- create data table containing the (multivariate) ranks of the observation in the raw forecast ensemble ---
  
  y = 2001:2010  #validation years
  m = 1:12       #validation months
  ens.size = 9
  
  #bias correct the raw ensemble
  
  DT_bc = DT[year %in% y & month %in% m,][,
             c(paste0("Ens",1:9)) := .SD + Bias_Est,
             .SDcols = c(paste0("Ens",1:9))]
  DT_bc[,c("Lon","Lat",paste0("SST",1:10),"SST_sd","Ens_bar","Ens_sd","grid_id","Bias_Est","SST_hat","oos_clim","oos_clim_past") := NULL]
  
  
  ym = unique(DT_bc[,YM])
  
  ranks.matrix = matrix(ym,nrow = length(ym),ncol = 1+4*(ens.size +1)) 
  #ncol: 1 col for YM, 4 methods of ranking, for each we get ranks for observation and each ensemble member 
  
  drm = dim(ranks.matrix)
  
  YM.ind = 0
  
  for(yearmonth in ym){
    print(paste0("YM = ",yearmonth,"/",ym[length(ym)]))
    YM.ind = YM.ind + 1
    fc_obs_mat = na.omit(DT_bc[YM == yearmonth,.SD, .SDcols = c("SST_bar",paste0("Ens",1:9))])
    
    # get ranks
    ranks.matrix[YM.ind,2:drm[2]] = c(mst.rank(as.matrix(fc_obs_mat)),
                                      mv.rank(as.matrix(fc_obs_mat)), 
                                      avg.rank(as.matrix(fc_obs_mat)),
                                      bd.rank(as.matrix(fc_obs_mat)))
  }
  
  names.vec = c("YM","mst.rk.obs",paste0("mst.r.",1:ens.size),
                "mv.rk.obs",paste0("mv.r.",1:ens.size),
                "av.rk.obs",paste0("av.r.",1:ens.size),
                "bd.rk.obs",paste0("bd.rk",1:ens.size))
  
  ranks = data.table(ranks.matrix)
  
  setnames(ranks, names.vec)
  
  # --- save ---
  
  save.dir = "~/PostClimDataNoBackup/SFE/Derived/Calibration/"
  
  save(ranks,file = paste0(save.dir,"ranks_bc.Rdata"))
  
  # ---- plotting ----
  pdf(file="./figures/rks_bc.pdf",width=8,height=2,points=12)
  par(mfrow=c(1,4),mex=0.5,oma = c(0,0,2.5,0),mar=c(2.5,2.5,2.5,2.5)+0.1,mgp=c(0.5,0,0))
  rhist.dt(ranks[,.(YM,mst.rk.obs)],hist_xlab = "minimum spanning tree")
  rhist.dt(ranks[,.(YM,mv.rk.obs)],hist_xlab = "multivariate")
  rhist.dt(ranks[,.(YM,av.rk.obs)],hist_xlab = "average")
  rhist.dt(ranks[,.(YM,bd.rk.obs)],hist_xlab = "band depth")
  
  title("Rank Histograms for bias corrected forecast",outer = TRUE)
  dev.off()
  
  ###################
  
  
  
  ##### rank histograms for PCA forecast
  
  DT = load_combined_wide(bias = TRUE)
  
  #DT = DT[Lon >= Lon_min & Lon <= Lon_max & Lat >= Lat_min & Lat <= Lat_max]
  
  
  y = 2001:2010  # validation years
  m = 1:12       # validation months
  ens.size = 9   # how often we generate PCA noise
  PCs = 15     # number of considered principal components
  
  setup_PCA(dt = DT, m = m, y = y, max_PCA_depth = 75)
  
  fc.dt = list()
  for(i in 1:ens.size){
    print(paste0(i,"/",ens.size))
    fc.dt[[i]] = forecast_PCA(m = m, y = y, PCA_depth = PCs, saveorgo = FALSE)[,forecast]
  }
  fc.dt = as.data.table(fc.dt)
  
  setnames(fc.dt,paste0("fc",1:ens.size))
  
  DT_pca = fc.dt[,c("YM","SST_bar","Lon","Lat"):= DT[year %in% y & month %in% m,.(YM,SST_bar,Lon,Lat)]]
  DT_pca = DT_pca[Lon >= Lon_min & Lon <= Lon_max & Lat >= Lat_min & Lat <= Lat_max]
  
  
  rm(fc.dt)
  
  ym = unique(DT_pca[,YM])
  
  ranks.matrix = matrix(ym,nrow = length(ym),ncol = 1+3*(ens.size +1)) 
  #ncol: 1 col for YM, 3 methods of ranking, for each we get ranks for observation and each ensemble member 
  
  drm = dim(ranks.matrix)
  
  YM.ind = 0
  
  for(yearmonth in ym){
    print(paste0("YM = ",yearmonth,"/",ym[length(ym)]))
    YM.ind = YM.ind + 1
    fc_obs_mat = na.omit(DT_pca[YM == yearmonth,.SD, .SDcols = c("SST_bar",paste0("fc",1:9))])
    
    # get ranks
    ranks.matrix[YM.ind,2:drm[2]] = c(mst.rank(as.matrix(fc_obs_mat)),
                                      avg.rank(as.matrix(fc_obs_mat)),
                                      bd.rank(as.matrix(fc_obs_mat)))
    rm(fc_obs_mat)
  }
  
  names.vec = c("YM","mst.rk.obs",paste0("mst.r.",1:ens.size),
                "av.rk.obs",paste0("av.r.",1:ens.size),
                "bd.rk.obs",paste0("bd.rk",1:ens.size))
  
  ranks = data.table(ranks.matrix)
  
  setnames(ranks, names.vec)
  
  # --- save ---
  
  save.dir = "~/PostClimDataNoBackup/SFE/Derived/Calibration/"
  
  save(ranks,file = paste0(save.dir,"ranks_pca_",PCs,"pcs.Rdata"))
  
  # ---- plotting ----
  pdf(file=paste0("./figures/rks_pca_",PCs,"pcs.pdf"),width=8,height=2,points=12)
  par(mfrow=c(1,3),mex=0.5,oma = c(0,0,2.5,0),mar=c(2.5,2.5,2.5,2.5)+0.1,mgp=c(0.5,0,0))
  rhist.dt(ranks[,.(YM,mst.rk.obs)],hist_xlab = "minimum spanning tree")
  rhist.dt(ranks[,.(YM,av.rk.obs)],hist_xlab = "average")
  rhist.dt(ranks[,.(YM,bd.rk.obs)],hist_xlab = "band depth")
  
  title(paste0("Rank Histograms for forecast with ",PCs," pcs"),outer = TRUE)
  dev.off()
  
  ###################
  
  # PCA forecasts with perturbed ensemble member rather than ensemble mean
  
  DT = load_combined_wide(bias = TRUE)
  
  y = 2001:2010  # validation years
  m = 1:12       # validation months
  ens.size = 9   # how often we generate PCA noise
  PCs = 5      # number of considered principal components
  
  setup_PCA(dt = DT, m = m, y = y, max_PCA_depth = 100)
  
  no.dt = list()
  for(i in 1:ens.size){
    print(paste0(i,"/",ens.size))
    no.dt[[i]] = forecast_PCA(m = m, y = y, PCA_depth = PCs, saveorgo = FALSE)[,noise]
  }
  no.dt = as.data.table(no.dt)
  
  setnames(no.dt,paste0("no",1:ens.size))
  
  DT_pca = no.dt[,c("YM","SST_bar","Ens1","Bias_Est") := DT[year %in% y & month %in% m,.(YM,SST_bar,Ens1,Bias_Est)]]
  
  DT_pca = DT_pca[,paste0("fc",1:ens.size) := .SD + Ens1 + Bias_Est,.SDcols = paste0("no",1:ens.size)]
  
  ym = unique(DT_pca[,YM])
  
  ranks.matrix = matrix(ym,nrow = length(ym),ncol = 1+3*(ens.size +1)) 
  #ncol: 1 col for YM, 3 methods of ranking, for each we get ranks for observation and each ensemble member 
  
  drm = dim(ranks.matrix)
  
  YM.ind = 0
  
  for(yearmonth in ym){
    print(paste0("YM = ",yearmonth,"/",ym[length(ym)]))
    YM.ind = YM.ind + 1
    fc_obs_mat = na.omit(DT_pca[YM == yearmonth,.SD, .SDcols = c("SST_bar",paste0("fc",1:9))])
    
    # get ranks
    ranks.matrix[YM.ind,2:drm[2]] = c(mst.rank(as.matrix(fc_obs_mat)),
                                      avg.rank(as.matrix(fc_obs_mat)),
                                      bd.rank(as.matrix(fc_obs_mat)))
    rm(fc_obs_mat)
  }
  
  names.vec = c("YM","mst.rk.obs",paste0("mst.r.",1:ens.size),
                "av.rk.obs",paste0("av.r.",1:ens.size),
                "bd.rk.obs",paste0("bd.rk",1:ens.size))
  
  ranks = data.table(ranks.matrix)
  
  setnames(ranks, names.vec)
  
  # --- save ---
  
  save.dir = "~/PostClimDataNoBackup/SFE/Derived/Calibration/"
  
  save(ranks,file = paste0(save.dir,"ranks_pca_em_",PCs,"pcs.Rdata"))
  
  # ---- plotting ----
  pdf(file=paste0("./figures/rks_pca_em_",PCs,"pcs.pdf"),width=8,height=2,points=12)
  par(mfrow=c(1,3),mex=0.5,oma = c(0,0,2.5,0),mar=c(2.5,2.5,2.5,2.5)+0.1,mgp=c(0.5,0,0))
  rhist.dt(ranks[,.(YM,mst.rk.obs)],hist_xlab = "minimum spanning tree")
  rhist.dt(ranks[,.(YM,av.rk.obs)],hist_xlab = "average")
  rhist.dt(ranks[,.(YM,bd.rk.obs)],hist_xlab = "band depth")
  
  title(paste0("RHs for ensemble member perturbed by ",PCs," pcs"),outer = TRUE)
  dev.off()
  
  ###################
  
