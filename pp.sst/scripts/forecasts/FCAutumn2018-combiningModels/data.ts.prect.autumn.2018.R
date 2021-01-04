
# loading the precipitation and surface temperature data for the autumn forecast in 2018 and bringing it into the standard shape

rm(list = ls())

options(max.print = 1e3)
setwd("~/temp")

for(var in c("prect","ts"))
{

  print(var)
  
dt = readRDS(file = paste0(var,"_1985_1989.rds"))

dt = as.data.table(dt)

setkeyv(dt, c("year","month","lon","lat"))

setnames(dt,c("lat","lon"),c("Lat","Lon"))

dt = dt[order(year,month,Lon,Lat)]


for(y in seq(1995,2015,by = 5))
{
  print(y)
  
  dt_new = readRDS(file = paste0(var,"_",y,"_",y+4,".rds"))
  
  dt_new = as.data.table(dt_new)
  
  setkeyv(dt_new, c("year","month","lon","lat"))
  
  setnames(dt_new,c("lat","lon"),c("Lat","Lon"))
  
  dt_new = dt_new[order(year,month,Lon,Lat)]
  
  dt = rbindlist(list(dt,dt_new))
  
}


Lon_trafo = function(x)
{
  x[x>180] = x[x>180]-360
  return(x)
}

dt[,Lon := Lon_trafo(Lon)]
dt = dt[order(year,month,Lon,Lat)]


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

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/",var,"/")
dir.create(save_dir,showWarnings = FALSE)

save(dt,file = paste0(save_dir,"/dt_combine_mr_wide.RData"))

#################################

# get forecasts:

dt_new = readRDS(file = paste0("forecast_",var,"_2017_2019.rds"))

dt_new = as.data.table(dt_new)

setkeyv(dt_new, c("year","month","lon","lat"))

setnames(dt_new,c("lat","lon"),c("Lat","Lon"))

dt_new = dt_new[order(year,month,Lon,Lat)]

dt = dt_new


Lon_trafo = function(x)
{
  x[x>180] = x[x>180]-360
  return(x)
}

dt[,Lon := Lon_trafo(Lon)]
dt = dt[order(year,month,Lon,Lat)]

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

save(dt,file = paste0(save_dir,"/dt_fc_ens30.RData"))

}

# some random Norway plotting...

lon_min = 00
lon_max = 40
lat_min = 55
lat_max = 75

dt_nor = dt[Lon < lon_max & Lon >lon_min & Lat < lat_max & Lat > lat_min,]

rr = range(dt_nor[year == 2018 & month %in% 8:12,ts_bar])

for(m in 8:12)
{
  plot_diagnostic(dt_nor[year == 2018 & month == m,.(Lon,Lat,ts_bar)],rr=rr,mn = paste0("month = ",m))
}

