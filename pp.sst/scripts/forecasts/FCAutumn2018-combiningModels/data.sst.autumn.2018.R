
# loading the SST data for the autumn forecast in 2018 and bringing it into the standard shape

rm(list = ls())

options(max.print = 1e3)
setwd("~/temp")

dt = readRDS(file = "sst_1985_1989.rds")

dt = as.data.table(dt)

setkeyv(dt, c("year","month","lon","lat"))

setnames(dt,c("lat","lon"),c("Lat","Lon"))

unique(dt_test[,.(year,month,lon,lat)])
dt = dt[order(year,month,Lon,Lat)]


for(y in seq(1995,2015,by = 5))
{
  
  dt_new = readRDS(file = paste0("sst_",y,"_",y+4,".rds"))
  
  dt_new = as.data.table(dt_new)
  
  setkeyv(dt_new, c("year","month","lon","lat"))
  
  setnames(dt_new,c("lat","lon"),c("Lat","Lon"))
  
  dt_new = dt_new[order(year,month,Lon,Lat)]
  
  dt = rbindlist(list(dt,dt_new))
  
}

# get forecasts:

dt_new = readRDS(file = "forecast_sst_2017_2019.rds")

dt_new = as.data.table(dt_new)

setkeyv(dt_new, c("year","month","lon","lat"))

setnames(dt_new,c("lat","lon"),c("Lat","Lon"))

dt_new = dt_new[order(year,month,Lon,Lat)]

dt = rbindlist(list(dt,dt_new))



Lon_trafo = function(x)
{
  x[x>180] = x[x>180]-360
  return(x)
}

dt[,Lon := Lon_trafo(Lon)]
dt = dt[order(year,month,Lon,Lat)]


sf_cols = c(paste0("sst",1:10),"sst_bar","sst_sd")
dt[,eval(sf_cols) := NULL]

names_old = c(paste0("ens",1:9),"ens_bar","ens_sd","noaa")
names_new = c(paste0("Ens",1:9),"Ens_bar","Ens_sd","SST_bar")

setnames(dt,names_old,names_new)

land_ids = dt[,which(is.na(Ens_bar))]

dt[land_ids,SST_bar := NA]

#### retracing grid_ids


lon_all = sort(dt[,unique(Lon)])
n_lon = length(lon_all)
cutoff_lon = c(-Inf,head(lon_all,-1) + diff(lon_all)/2)
f_lon = approxfun(cutoff_lon, 1:n_lon, method="constant", rule = 2)

lat_all = sort(dt[,unique(Lat)])
n_lat = length(lat_all)
cutoff_lat = c(-Inf, head(lat_all,-1) + diff(lat_all)/2)
f_lat = approxfun(cutoff_lat, 0:(n_lat - 1) * n_lon, method="constant", rule = 2)
dt[,grid_id := f_lon(Lon) + f_lat(Lat)]


### saving

save_dir = "~/PostClimDataNoBackup/SFE/Derived/"

save(dt,file = paste0(save_dir,"/dt_combine_mr_wide.RData"))
