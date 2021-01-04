

rm(list = ls())

setwd('~/PostClimDataNoBackup/SFE/NORCE')


dt =  readRDS('sst_hindcast.rds')

dt = as.data.table(dt)

setkeyv(dt, c("year","month","lon","lat"))

setnames(dt,c("lat","lon"),c("Lat","Lon"))

setnames(dt,c('noaa_sst',paste0('ens',1:9),'ens_bar','ens_sd'),c('SST_bar',paste0('Ens',1:9),'Ens_bar','Ens_sd'))

dt[,c(paste0('sst',1:10),'sst_bar','sst_sd'):= NULL]

unique(dt[,.(year,month,Lon,Lat)])

setorder(dt,year,month,Lon,Lat)


land_ids = dt[,which(is.na(Ens_bar))]

# convert NaNs to missing values

conv_nan_to_na = function(x)
{
  x[is.nan(x)] <- NA
  return(x)
}


dt[,lapply(.SD,conv_nan_to_na),.SDcols = colnames(dt)]

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


# save
