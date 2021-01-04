rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "south_afr" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"univariate.model.comparison.RData"))

PCA_dir = paste0(save_dir,"PCA/ssq/")

m = 1

PC_sa = get_PCs(dt = DT[year == 2006 & m == 1,],2006,1,cov_dir = PCA_dir)

ws = 2

La = DT[,unique(Lat)]
La = La[1:(length(La)-ws)]
Lo = DT[,unique(Lon)]
Lo = Lo[1:(length(Lo)-ws)]

Cov_check = as.data.table(expand.grid(La,Lo,La,Lo))
setnames(Cov_check, c("Lat1","Lon1","Lat2","Lon2"))

# strong negative correllation on relative short distance seems to be at Lon 23 to 25 in the windows -43,-40 and -36,-34
mean_cov = c()
for( i in 1:nrow(Cov_check))
{
  if(i %% 1000 == 0)
  {
    print(paste0("i = ",i,"/",nrow(Cov_check)))
  }

  
  la_1 = Cov_check[,Lat1][i]
  lo_1 = Cov_check[,Lon1][i]
  la_2 = Cov_check[,Lat2][i]
  lo_2 = Cov_check[,Lon2][i]
  
  
  lat_1 = c(la_1,la_1+ws)
  lon_1 = c(lo_1,lo_1+ws)
  
  lat_2 = c(la_2,la_2+ws)
  lon_2 = c(lo_2,lo_2+ws)  
  
  
  PC_sa_test_1 = PC_sa[min(lat_1) <= Lat & Lat <= max(lat_1) & min(lon_1) <= Lon & Lon <= max(lon_1) & month == m, ]
  PC_sa_test_2 = PC_sa[min(lat_2) <= Lat & Lat <= max(lat_2) & min(lon_2) <= Lon & Lon <= max(lon_2) & month == m, ]
  
  PC_sa_test = rbindlist(list(PC_sa_test_1[month == m,],PC_sa_test_2[month == m,]))
  
  check_sth = as.matrix(PC_sa_test_1[,.SD,.SDcols = paste0("PC",1:4)]) * as.matrix(PC_sa_test_2[,.SD,.SDcols = paste0("PC",1:4)])
  mean_cov = c(mean_cov, mean(rowSums(check_sth),na.rm = TRUE)) 
}

Cov_check[,mc:= mean_cov]

distances =  geosphere::distHaversine(as.matrix(Cov_check[,.(Lon1,Lat1)]),as.matrix(Cov_check[,.(Lon2,Lat2)]))

Cov_check[,Dist:=distances]

# get number of NAs in the corresponding grid

DT = DT[month == m,]

for(lat in La){
  for(lon in  Lo){
    sDT = DT[Lat >= lat & Lat <= lat + ws & Lon >= lon & Lon <= lon + ws & year == min(year), ]
    DT[Lat == lat & Lon ==lon, num_NAs := sDT[is.na(SST_bar) | is.na(Ens_bar),.N]]    
  }
}

DT_coor = DT[year == min(year) & Lat %in% La & Lon %in% Lo,.(Lon,Lat,num_NAs)]
lDT = DT_coor[,.N]

ind_1 = rep(1:lDT,lDT)
ind_2 = as.vector(t(matrix(ind_1,ncol = lDT)))

Cov_check[,num_NAs := DT_coor[,num_NAs][ind_1] + DT_coor[,num_NAs][ind_2] ]

save(Cov_check,file = paste0(save_dir,"Cov_check.RData"))

#show minimal negative correlation as function of distance

Dist_bins = seq(0,Cov_check[num_NAs == 0,max(Dist)],length.out = 50)

min_cov = c()
for(d in Dist_bins)
{
  print(d)
  min_cov = c(min_cov,Cov_check[num_NAs == 0 & Dist <= d, min(mc)])
}

par(mfrow = c(1,1))

pdf(paste0(plot_dir,"cov_by_dist.pdf"))
plot(Dist_bins,min_cov,main = "minimal covariance by distance", xlab = "Distance", ylab = "min. cov.")
dev.off()


#### handpick a suitable strongly negative covariance that exists at relatively small distance:

picked_indices = c(9,16,26)

par(mfrow = c(1,1))

pdf(paste0(plot_dir,"cov_by_dist_sel.pdf"))
plot(Dist_bins,min_cov,main = "minimal covariance by distance", xlab = "Distance", ylab = "min. cov.")

for(ind in picked_indices)
{
  points(x = Dist_bins[ind],
         y = min_cov[ind],
         col = "black",
         bg = "black",
         pch = 21)
}

dev.off()

##########

ind = 26

min_cov_sel = min_cov[ind]

# find closest locations with this covariance
all_test_locs = Cov_check[mc < -.25 & num_NAs == 0,]

test_loc = Cov_check[mc ==min_cov_sel & num_NAs == 0,][Dist == min(Dist),]


lat_1 = c(test_loc[1,Lat1],test_loc[1,Lat1] + ws)
lat_2 = c(test_loc[1,Lat2],test_loc[1,Lat2] + ws)

lon_1 = c(test_loc[1,Lon1],test_loc[1,Lon1] + ws)
lon_2 = c(test_loc[1,Lon2],test_loc[1,Lon2] + ws)

DT_test_1 = DT[Lat >= min(lat_1) & Lat <= max(lat_1) & Lon >= min(lon_1) & Lon <= max(lon_1) & month == m,]
DT_test_2 = DT[Lat >= min(lat_2) & Lat <= max(lat_2) & Lon >= min(lon_2) & Lon <= max(lon_2) & month == m,]


DT_test = rbindlist(list(DT_test_1,DT_test_2))



#####################################################################
###### generate forecasts by different methods for these areas ######
#####################################################################

Lon_1 = min(lon_1) + ws/2
Lat_1 = min(lat_1) + ws/2

Lon_2 = min(lon_2) + ws/2
Lat_2 = min(lat_2) + ws/2

size_temp = 0

DT_temp_1 = DT[abs(Lon-Lon_1) <= size_temp & abs(Lat-Lat_1) <= size_temp ]
DT_temp_2 = DT[abs(Lon-Lon_2) <= size_temp & abs(Lat-Lat_2) <= size_temp ]

DT_temp = rbindlist(list(DT_temp_1,DT_temp_2))

coords = unique(DT_temp[,.(Lon,Lat)])

cov_training = cov(DT_temp[Lon == min(Lon) & month ==1 & year %in% 1986:2005,
                           SST_bar-SST_hat],DT_temp[Lon == max(Lon) & month == 1 & year %in% 1986:2005,SST_bar-SST_hat])
cov_full = cov(DT_temp[Lon == min(Lon) & month ==1 & year %in% 1986:2010,
                       SST_bar-SST_hat],DT_temp[Lon == max(Lon) & month == 1 & year %in% 1986:2010,SST_bar-SST_hat])

pdf(paste0(plot_dir,"neg_cor_ex3.pdf"))

par(mfrow = c(2,1))

rr = range(DT_temp[month == 1 & year %in% 1986:2010,SST_bar-SST_hat])
plot(DT_temp[Lon == min(Lon) & month ==1 & year %in% 1986:2010,.(year,SST_bar-SST_hat)], ylim = rr,col = "blue", ylab = "residual", main = "example residuals at two locations with strong negative correlation" )
points(DT_temp[Lon == max(Lon) & month ==1 & year %in% 1986:2010,.(year,SST_bar-SST_hat)],col = "darkgreen")
abline(h = 0)
abline(v = 2005.5)

plot(DT_temp[Lon == min(Lon) & month ==1 & year %in% 1986:2010,.(year,SST_bar-SST_hat)], ylim = rr,col = "blue",ylab = "residual",xlab = "", type = "l")
lines(DT_temp[Lon == max(Lon) & month ==1 & year %in% 1986:2010,.(year,SST_bar-SST_hat)],col = "darkgreen")
abline(h = 0)
abline(v = 2005.5)


title(sub = paste0("cov. training period = ",round(cov_training,2),",   cov. all years = ",round(cov_full,2),",   dist = ",round(Dist_bins[ind]/1000,1), " km"))

dev.off()


#get truely strongly negatively correlated locations#



##############################################################

name_abbr = "south_afr" # for south of africa

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"univariate.model.comparison.RData"))

plot_dir = paste0(plot_dir,"upscaling/")
dir.create(plot_dir,showWarnings = FALSE)
save_dir = paste0(save_dir,"upscaling/")
dir.create(save_dir,showWarnings = FALSE)

months = m


ind_m = DT[month == m,which((abs(Lon-Lon_1) <= size_temp & abs(Lat-Lat_1) <= size_temp) | (abs(Lon-Lon_2) <= size_temp &abs(Lat-Lat_2) <= size_temp))]
ind_my = DT[month == m & year == min(year),which((abs(Lon-Lon_1) <= size_temp & abs(Lat-Lat_1) <= size_temp) | (abs(Lon-Lon_2) <= size_temp &abs(Lat-Lat_2) <= size_temp))]


##########################################
##### testing different PCA versions #####
##########################################

# get climatology

DT[year %in% training_years, clim := mean(SST_bar), by = .(Lat,Lon,month)]

for(y in validation_years)
{
  DT[year == y, clim := DT[year == min(training_years),clim]]
}

# get clim_sigma

DT[year %in% training_years, clim_sigma := sd(SST_bar), by = .(Lat,Lon,month)]

for(y in validation_years)
{
  DT[year == y, clim_sigma := DT[year == min(training_years),clim_sigma]]
}



versions = c("aggr_by_season","sum_of_squares","wrt_ens_mean")
version_titles = c("aggr. by season","w.r.t. ensemble members","w.r.t. ensemble mean")
version_abbrs = c("abs","ssq","scm")


PCA_0_dir = paste0(save_dir,"PCA/")
dir.create(PCA_0_dir, showWarnings = FALSE)



 for(ver_ind in 1:length(versions))
 {
  version = versions[ver_ind]
  version_title = version_titles[ver_ind]
  version_abbr = version_abbrs[ver_ind]
  
  print(paste0("version = ",version))
  
  PCA_dir = paste0(PCA_0_dir,version_abbr,"/")
  dir.create(PCA_dir, showWarnings = FALSE)
  
  training_years = DT[!(year %in% validation_years),unique(year)]
  
  for_res_cov(Y = training_years,
              M = months,
              dt = DT, 
              save_dir = PCA_dir,
              ens_size = ens_size,
              version = version)
  
  
  
  ### range of PCs to test ###
  
  
  if(version == "aggr_by_season")
  {
    PCs = c(1:10,15,20)  
  }
  if(version == "sum_of_squares")
  {
    PCs = c(1:10,20,50)  
  }
  if(version == "wrt_ens_mean")
  {
    PCs = c(1:10)  
  }
  
  #### variogram scores computation: ####
  
  
  # setup: get principal components and marginal variances for the given month m:
    
  print(paste0("m = ",m))
  load(file = paste0(PCA_dir,"CovRes_mon",m,".RData"))
  PCA = irlba::irlba(res_cov, nv = max(PCs),fastpath = FALSE)
    
  land_ids = which(DT[year == min(year) & month == min(month),is.na(Ens_bar) | is.na(SST_bar)])
    
    PCA_DT = DT[year == min(year) & month == min(month),][-land_ids,.(Lon,Lat,grid_id)]
    
    for(d in  min(PCs):max(PCs))
    {
      PCA_DT [,paste0("PC",d) := PCA$d[d]*PCA$u[,d]]
    } 
    
    # also get marginal variances
    variances = list()
    d = min(PCs) 
    vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
    variances[[d]] = vec^2
    
    if(length(PCs)>1)
    {
      for(d in (min(PCs)+1):max(PCs)){
        vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
        variances[[d]] = variances[[d-1]] + vec^2
      }
    }
    names(variances) = paste0("var",min(PCs):max(PCs))
    PCA_DT = data.table(PCA_DT,rbindlist(list(variances)))    
    
    # reduce to relevant window
    
    PCA_DT = PCA_DT[(abs(Lon-Lon_1) <= size_temp & abs(Lat-Lat_1) <= size_temp) | (abs(Lon-Lon_2) <= size_temp &abs(Lat-Lat_2) <= size_temp),]
    
    
    # now move to getting the variogram scores:
    
    # without marginal correction:
    
    dummy_fct = function(y)
    {
      var_sc_PCA_standardized(m, y, DT_temp, PCA = PCA, PCA_DT = PCA_DT, dvec = PCs, ens_size = ens_size,
                              weighted = FALSE,
                              file_name = paste0("var_sc_by_PC_stan"),
                              marginal_correction = FALSE, 
                              cov_dir = PCA_dir, save_dir = PCA_dir)
    }
    parallel::mclapply(validation_years,dummy_fct,mc.cores = length(validation_years))
    
    # with marginal correction:
    
    dummy_fct = function(y)
    {
      print(c(paste0("y = ",y),paste0("m = ",m)))
      var_sc_PCA_standardized(m, y, DT_temp, PCA = PCA, PCA_DT = PCA_DT, dvec = PCs, ens_size = ens_size,
                              weighted = FALSE,
                              file_name = paste0("var_sc_by_PC_stan"),
                              cov_dir = PCA_dir, save_dir = PCA_dir)
    }
    parallel::mclapply(validation_years,dummy_fct,mc.cores = length(validation_years))
  

  
  ###### combine: ######
  
  scores_stan_nmc = list()
  k=0
  for(m in months){
    for(y in validation_years)
    {k=k+1
    load(file = paste0(PCA_dir,"var_sc_by_PC_stan_nmc_m",m,"_y",y,".RData"))
    scores_stan_nmc[[k]] = scores
    }
  }
  scores_stan_nmc = rbindlist(scores_stan_nmc)
  
  scores_stan_mc = list()
  k=0
  for(m in months){
    for(y in validation_years)
    {k=k+1
    load(file = paste0(PCA_dir,"var_sc_by_PC_stan_m",m,"_y",y,".RData"))
    scores_stan_mc[[k]] = scores
    }
  }
  scores_stan_mc = rbindlist(scores_stan_mc)
  
  mean_sc = scores_stan_mc[,mean_sc := mean(sc), by = d][,.(d,mean_sc)]
  mean_sc = unique(mean_sc)
  
  
  mean_sc_nmc = scores_stan_nmc[,mean_sc := mean(sc), by = d][,.(d,mean_sc)]
  mean_sc_nmc = unique(mean_sc_nmc)
  
  # save
  
  save(scores_stan_mc,scores_stan_nmc,mean_sc,mean_sc_nmc,file = paste0(PCA_dir,"variogram_scores_stan.RData"))
  
}




##### geostationary #####

geostat_dir = paste0(save_dir, "GeoStat/")
dir.create(geostat_dir, showWarnings = FALSE)

geostationary_training(dt = DT, save_dir = geostat_dir,m = months)

load(paste0(geostat_dir,"variogram_exp_m",m,".RData")) 

sp_sm <- sp::SpatialPoints(cbind(x=coords[, Lon],
                                 y=coords[, Lat]), 
                        proj4string = sp::CRS("+proj=longlat +datum=WGS84"))

Dist_sm <- sp::spDists(sp_sm, longlat = TRUE)

# variogram scores without marginal correction:
  
dummy_fct = function(y)
{
    var_sc_geoStat_new( DT_temp, m, y, Mod = Mod, Dist = Dist_sm,
                        weighted = FALSE,
                        file_name = "var_sc",
                        save_dir = geostat_dir, data_dir = geostat_dir,
                        mar_var_cor = FALSE)
  }
  parallel::mclapply(validation_years,dummy_fct,mc.cores = length(validation_years))
  
  # with marginal correction:
  
  dummy_fct = function(y)
  {
    var_sc_geoStat_new( DT_temp,m, y, Mod = Mod, Dist = Dist_sm,
                        weighted = FALSE,
                        file_name = "var_sc",
                        save_dir = geostat_dir, data_dir = geostat_dir,
                        mar_var_cor = TRUE)
  }
  parallel::mclapply(validation_years,dummy_fct,mc.cores = length(validation_years))



####### combine #########

scores_geostat_nmc = list()
k=0
for(y in validation_years)
  {k=k+1
  load(file = paste0(geostat_dir,"var_sc_nmc_m",m,"_y",y,".RData"))
  scores_geostat_nmc[[k]] = scores
  }

scores_geostat_nmc = rbindlist(scores_geostat_nmc)

scores_geostat_mc = list()
k=0
for(y in validation_years)
  {k=k+1
  load(file = paste0(geostat_dir,"var_sc_m",m,"_y",y,".RData"))
  scores_geostat_mc[[k]] = scores
  }

scores_geostat_mc = rbindlist(scores_geostat_mc)


mean_geostat_sc = scores_geostat_mc[,mean_sc := mean(sc)][,.(mean_sc)]
mean_geostat_sc = unique(mean_geostat_sc)


mean_geostat_sc_nmc = scores_geostat_nmc[,mean_sc := mean(sc)][,.(mean_sc)]
mean_geostat_sc_nmc = unique(mean_geostat_sc_nmc)

# save

save(scores_geostat_mc,scores_geostat_nmc,mean_geostat_sc,mean_geostat_sc_nmc,file = paste0(geostat_dir,"variogram_scores.RData"))


##### the same for standardized variables #####

# setup: get principal components and marginal variances for the given month m:
  
load(paste0(geostat_dir,"variogram_exp_m",m,".RData")) 
  
  # variogram scores without marginal correction:
  
  dummy_fct = function(y)
  {
    var_sc_geoStat_standardized( DT_temp, m, y, Mod = Mod, Dist = Dist_sm,
                                 weighted = FALSE,
                                 save_dir = geostat_dir, data_dir = geostat_dir,
                                 mar_var_cor = FALSE)
  }
  parallel::mclapply(validation_years,dummy_fct,mc.cores = length(validation_years))
  
  # with marginal correction:
  
  dummy_fct = function(y)
  {
    var_sc_geoStat_standardized( DT_temp,m, y, Mod = Mod, Dist = Dist_sm,
                                 weighted = FALSE,
                                 save_dir = geostat_dir, data_dir = geostat_dir,
                                 mar_var_cor = TRUE)
  }
  parallel::mclapply(validation_years,dummy_fct,mc.cores = length(validation_years))



####### combine #########

scores_stan_geostat_nmc = list()
k=0
for(y in validation_years)
  {k=k+1
  load(file = paste0(geostat_dir,"var_sc_stan_nmc_m",m,"_y",y,".RData"))
  scores_stan_geostat_nmc[[k]] = scores
  }

scores_stan_geostat_nmc = rbindlist(scores_stan_geostat_nmc)

scores_stan_geostat_mc = list()
k=0
for(y in validation_years)
  {k=k+1
  load(file = paste0(geostat_dir,"var_sc_stan_m",m,"_y",y,".RData"))
  scores_stan_geostat_mc[[k]] = scores
  }

scores_stan_geostat_mc = rbindlist(scores_stan_geostat_mc)


mean_geostat_sc_stan = scores_stan_geostat_mc[,mean_sc := mean(sc)][,.(mean_sc)]
mean_geostat_sc_stan = unique(mean_geostat_sc_stan)


mean_geostat_sc_stan_nmc = scores_stan_geostat_nmc[,mean_sc := mean(sc)][,.(mean_sc)]
mean_geostat_sc_stan_nmc = unique(mean_geostat_sc_stan_nmc)

# save

save(scores_stan_geostat_mc,scores_stan_geostat_nmc,mean_geostat_sc_stan,mean_geostat_sc_stan_nmc,file = paste0(geostat_dir,"variogram_scores_stan.RData"))







###############################
### plotting and comparison ###
###############################

# for(ver_ind in 1:length(versions))
# {
#   version = versions[ver_ind]
#   version_title = version_titles[ver_ind]
#   version_abbr = version_abbrs[ver_ind]
#   
#   PCA_dir = paste0(save_dir,"PCA/",version_abbr,"/")
#   
#   load(file = paste0(PCA_dir,"variogram_scores.RData"))
#   
#   assign(paste0("mean_sc_",version_abbr),mean_sc)
#   assign(paste0("mean_sc_nmc_",version_abbr),mean_sc_nmc)
# }
# 
# 
# ##################################
# # plot with marginal correction: #
# ##################################
# 
# plot_dir = paste0(plot_dir,"upscaling/")
# 
# pdf(paste0(plot_dir,"/mean_variogram_scores_",size_temp,".pdf"))
# rr = range(c(mean_sc_abs[[2]],mean_sc_ssq[[2]],mean_sc_scm[[2]],mean_geostat_sc))
# 
# plot(x = mean_sc_ssq[[1]],
#      y = mean_sc_ssq[[2]],
#      ylim = rr,
#      type = "b",
#      col = "blue",
#      main = paste0("mean variogram scores"),
#      xlab = "number of principal components",
#      ylab = "mean score")
# 
# #---- add minima: -----
# 
# 
# lines(x = mean_sc_abs[[1]],
#       y = mean_sc_abs[[2]],
#       type = "b",
#       col = "darkred")
# 
# lines(x = mean_sc_scm[[1]],
#       y = mean_sc_scm[[2]],
#       type = "b",
#       col = "darkgreen")
# 
# 
# # --- add geostat value: ----
# 
# abline(h = mean_geostat_sc[,mean_sc], lty = "dashed", col = adjustcolor("black"))
# 
# #abline(h = sc_ECC[,mean(sc)], lty = "dashed", col = adjustcolor("pink"))
# 
# 
# legend("topright",legend = c("w.r.t. ens. mem.","w.r.t ens. mean","aggr. by season","geostat"),
#        col = c("blue","darkgreen","darkred","black"),lty = c(1,1,1,2))
# # legend("topright",legend = c("PCA, m.c.v.","PCA","geostat, m.c.v.","geostat","ECC"),
# #        col = c("blue","darkgreen","darkred","black","pink"),lty = c(1,1,2,2,2))
# dev.off()
# 
# 
# 
# #####################################
# # plot without marginal correction: #
# #####################################
# 
# pdf(paste0(plot_dir,"/mean_variogram_scores_nmc_",size_temp,".pdf"))
# rr = range(c(mean_sc_nmc_abs[[2]],mean_sc_nmc_ssq[[2]],mean_sc_nmc_scm[[2]],mean_geostat_sc_nmc,mean_geostat_sc))
# 
# plot(x = mean_sc_nmc_ssq[[1]],
#      y = mean_sc_nmc_ssq[[2]],
#      ylim = rr,
#      type = "b",
#      col = "blue",
#      main = paste0("mean variogram scores without mar. cor."),
#      xlab = "number of principal components",
#      ylab = "mean score")
# 
# #---- add minima: -----
# 
# 
# lines(x = mean_sc_nmc_abs[[1]],
#       y = mean_sc_nmc_abs[[2]],
#       type = "b",
#       col = "darkred")
# 
# lines(x = mean_sc_nmc_scm[[1]],
#       y = mean_sc_nmc_scm[[2]],
#       type = "b",
#       col = "darkgreen")
# 
# 
# # --- add geostat value: ----
# 
# abline(h = mean_geostat_sc_nmc[,mean_sc], lty = "dashed", col = adjustcolor("black"))
# abline(h = mean_geostat_sc[,mean_sc], lty = "dashed", col = adjustcolor("black"))
# 
# #abline(h = sc_ECC[,mean(sc)], lty = "dashed", col = adjustcolor("pink"))
# 
# 
# legend("topright",legend = c("w.r.t. ens. mem.","w.r.t ens. mean","aggr. by season","geostat"),
#        col = c("blue","darkgreen","darkred","black"),lty = c(1,1,1,2))
# # legend("topright",legend = c("PCA, m.c.v.","PCA","geostat, m.c.v.","geostat","ECC"),
# #        col = c("blue","darkgreen","darkred","black","pink"),lty = c(1,1,2,2,2))
# dev.off()
# 

############# plotting standardized versions ######################

### plotting and comparison ###

for(ver_ind in 1:length(versions))
{
  version = versions[ver_ind]
  version_title = version_titles[ver_ind]
  version_abbr = version_abbrs[ver_ind]
  
  PCA_dir = paste0(save_dir,"PCA/",version_abbr,"/")
  
  load(file = paste0(PCA_dir,"variogram_scores_stan.RData"))
  
  assign(paste0("mean_sc_",version_abbr),mean_sc)
  assign(paste0("mean_sc_nmc_",version_abbr),mean_sc_nmc)
}


##################################
# plot with marginal correction: #
##################################

pdf(paste0(plot_dir,"/mean_variogram_scores_stan_",size_temp,".pdf"))
rr = range(c(mean_sc_abs[[2]],mean_sc_ssq[[2]],mean_sc_scm[[2]],mean_geostat_sc_stan))

plot(x = mean_sc_ssq[[1]],
     y = mean_sc_ssq[[2]],
     ylim = rr,
     type = "b",
     col = "blue",
     main = paste0("mean variogram scores"),
     xlab = "number of principal components",
     ylab = "mean score")

#---- add minima: -----


lines(x = mean_sc_abs[[1]],
      y = mean_sc_abs[[2]],
      type = "b",
      col = "darkred")

lines(x = mean_sc_scm[[1]],
      y = mean_sc_scm[[2]],
      type = "b",
      col = "darkgreen")


# --- add geostat value: ----

abline(h = mean_geostat_sc_stan[,mean_sc], lty = "dashed", col = adjustcolor("black"))

#abline(h = sc_ECC[,mean(sc)], lty = "dashed", col = adjustcolor("pink"))


legend("topright",legend = c("w.r.t. ens. mem.","w.r.t ens. mean","aggr. by season","geostat"),
       col = c("blue","darkgreen","darkred","black"),lty = c(1,1,1,2))
# legend("topright",legend = c("PCA, m.c.v.","PCA","geostat, m.c.v.","geostat","ECC"),
#        col = c("blue","darkgreen","darkred","black","pink"),lty = c(1,1,2,2,2))
dev.off()



#####################################
# plot without marginal correction: #
#####################################

pdf(paste0(plot_dir,"/mean_variogram_scores_stan_nmc_",size_temp,".pdf"))
rr = range(c(mean_sc_nmc_abs[[2]],mean_sc_nmc_ssq[[2]],mean_sc_nmc_scm[[2]],mean_geostat_sc_stan))

plot(x = mean_sc_nmc_ssq[[1]],
     y = mean_sc_nmc_ssq[[2]],
     ylim = rr,
     type = "b",
     col = "blue",
     main = paste0("mean variogram scores without mar. cor."),
     xlab = "number of principal components",
     ylab = "mean score")

#---- add minima: -----


lines(x = mean_sc_nmc_abs[[1]],
      y = mean_sc_nmc_abs[[2]],
      type = "b",
      col = "darkred")

lines(x = mean_sc_nmc_scm[[1]],
      y = mean_sc_nmc_scm[[2]],
      type = "b",
      col = "darkgreen")


# --- add geostat value: ----


abline(h = mean_geostat_sc_stan[,mean_sc], lty = "dashed", col = adjustcolor("black"))

#abline(h = sc_ECC[,mean(sc)], lty = "dashed", col = adjustcolor("pink"))


legend("topright",legend = c("w.r.t. ens. mem.","w.r.t ens. mean","aggr. by season","geostat"),
       col = c("blue","darkgreen","darkred","black"),lty = c(1,1,1,2))
# legend("topright",legend = c("PCA, m.c.v.","PCA","geostat, m.c.v.","geostat","ECC"),
#        col = c("blue","darkgreen","darkred","black","pink"),lty = c(1,1,2,2,2))
dev.off()


pdf(paste0(plot_dir,"values_",size_temp,".pdf"))
 rr = range(DT_temp[month == 1 & year %in% 1986: 2010,SST_bar-SST_hat])
 plot(DT_temp[Lon == 48.5 & month ==1 & year %in% 1986:2010,.(year,SST_bar-SST_hat)], ylim = rr,col = "blue")
 points(DT_temp[Lon == 46.5 & month ==1 & year %in% 1986:2010,.(year,SST_bar-SST_hat)],col = "darkgreen")
 abline(h = 0)
 abline(v = 2005.5)
 dev.off()