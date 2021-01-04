
# importing the SST forecasts 1981:2019 and saving them


rm(list = ls())
setwd("~/PostClimDataNoBackup/SFE/GCFS1")

library(ncdf4)
library(data.table)
library(PostProcessing)

# get the first ensemble member for the forecast:

mnum = 1

if(mnum < 10)
{
  mnumstr = paste0("0",mnum)
}else{
  mnumstr = mnum
}

ncname = paste0("./FCAST_082018/Mem",mnumstr,"_GCFS1_temp2_mm_2016-2018_smon07.nc")

#load the data:

nc = nc_open(ncname)

Lons = ncvar_get(nc,"lon")
Lats = ncvar_get(nc,"lat")
t = ncvar_get(nc,"time")
temp2 = ncvar_get(nc,"temp2")


# convert time and get years and months:

t_new = as.POSIXlt(3600*24*t,origin = "1978-6-30 23:50:00",tz = "CEST")
months = as.integer(format(t_new, format = "%m"))
years = as.integer(format(t_new, format = "%Y"))

#pool everything into a data table

DT = list()
for(i in 1: length(t_new))
{
  timeslize = as.data.table(expand.grid(Lons,Lats))
  timeslize[,temp2 := as.vector(temp2[,,i]-273.15)]
  setnames(timeslize,c("Lon","Lat",paste0("Ens",mnum) ))
  timeslize[,"year" := years[i]]
  timeslize[,"month" := months[i]]
  DT = rbindlist(list(DT,timeslize))
}

setkey(DT,year,month,Lon,Lat)

# for some reason there are location duplicates? Perhaps this has to do with measurements at different depths?

DT = unique(DT)


### do the same for the remaining ensemble members: ###

for(mnum in 2:15)
{

print(mnum)
  
  if(mnum < 10)
  {
    mnumstr = paste0("0",mnum)
  }else{
    mnumstr = mnum
  }
  
  ncname = paste0("./FCAST_082018/Mem",mnumstr,"_GCFS1_temp2_mm_2016-2018_smon07.nc")
  
  #load the data:
  
  nc = nc_open(ncname)
  
  Lons = ncvar_get(nc,"lon")
  Lats = ncvar_get(nc,"lat")
  t = ncvar_get(nc,"time")
  temp2 = ncvar_get(nc,"temp2")
  
  # convert time and get years and months:
  
  t_new = as.POSIXlt(24*3600*t,origin = "1978-6-30 23:50:00",tz = "CEST")
  months = as.integer(format(t_new, format = "%m"))
  years = as.integer(format(t_new, format = "%Y"))
  
  #pool everything into a data table
  
  
  DT_temp = list()
  for(i in 1: length(t_new))
  {
    timeslize = as.data.table(expand.grid(Lons,Lats))
    timeslize[,temp2 := as.vector(temp2[,,i]-273.15)]
    setnames(timeslize,c("Lon","Lat",paste0("Ens",mnum) ))
    timeslize[,"year" := years[i]]
    timeslize[,"month" := months[i]]
    DT_temp = rbindlist(list(DT_temp,timeslize))
  }
  
  setkey(DT_temp,year,month,Lon,Lat)
  DT_temp = unique(DT_temp)
  
  DT = merge(DT,DT_temp,by = c("Lon","Lat","year","month"))
  DT = DT[order(year,month,Lon,Lat)]
}

DT_fc = copy(DT)


########### forecasts done, move to hindcasts ############

# get the first ensemble member for the forecast:

mnum = 1

if(mnum < 10)
{
  mnumstr = paste0("0",mnum)
}else{
  mnumstr = mnum
}

ncname = paste0("./HCAST_082018/Mem",mnumstr,"_GCFS1_temp2_mm_1981-2015_smon07.nc")

#load the data:

nc = nc_open(ncname)

Lons = ncvar_get(nc,"lon")
Lats = ncvar_get(nc,"lat")
t = ncvar_get(nc,"time")
temp2 = ncvar_get(nc,"temp2")


# convert time and get years and months:

t_new = as.POSIXlt(3600*24*t,origin = "1978-6-30 23:50:00",tz = "CEST")
months = as.integer(format(t_new, format = "%m"))
years = as.integer(format(t_new, format = "%Y"))

#pool everything into a data table

DT = list()
for(i in 1: length(t_new))
{
  timeslize = as.data.table(expand.grid(Lons,Lats))
  timeslize[,temp2 := as.vector(temp2[,,i]-273.15)]
  setnames(timeslize,c("Lon","Lat",paste0("Ens",mnum) ))
  timeslize[,"year" := years[i]]
  timeslize[,"month" := months[i]]
  DT = rbindlist(list(DT,timeslize))
}

setkey(DT,year,month,Lon,Lat)

# for some reason there are location duplicates? Perhaps this has to do with measurements at different depths?

DT = unique(DT)


### do the same for the remaining ensemble members: ###

for(mnum in 2:15)
{
  
  print(mnum)
  
  if(mnum < 10)
  {
    mnumstr = paste0("0",mnum)
  }else{
    mnumstr = mnum
  }
  
  ncname = paste0("./HCAST_082018/Mem",mnumstr,"_GCFS1_temp2_mm_1981-2015_smon07.nc")
  
  #load the data:
  
  nc = nc_open(ncname)
  
  Lons = ncvar_get(nc,"lon")
  Lats = ncvar_get(nc,"lat")
  t = ncvar_get(nc,"time")
  temp2 = ncvar_get(nc,"temp2")
  
  # convert time and get years and months:
  
  t_new = as.POSIXlt(24*3600*t,origin = "1978-6-30 23:50:00",tz = "CEST")
  months = as.integer(format(t_new, format = "%m"))
  years = as.integer(format(t_new, format = "%Y"))
  
  #pool everything into a data table
  
  
  DT_temp = list()
  for(i in 1: length(t_new))
  {
    timeslize = as.data.table(expand.grid(Lons,Lats))
    timeslize[,temp2 := as.vector(temp2[,,i]-273.15)]
    setnames(timeslize,c("Lon","Lat",paste0("Ens",mnum) ))
    timeslize[,"year" := years[i]]
    timeslize[,"month" := months[i]]
    DT_temp = rbindlist(list(DT_temp,timeslize))
  }
  
  setkey(DT_temp,year,month,Lon,Lat)
  DT_temp = unique(DT_temp)
  
  DT = merge(DT,DT_temp,by = c("Lon","Lat","year","month"))
  DT = DT[order(year,month,Lon,Lat)]
}


DT = DT[order(year,month,Lon,Lat)]


DT = rbindlist(list(DT,DT_fc))




#################################

# get grid map - they use a different for 2mt than for sst

grid_GCFS1_2mt = DT[year == 2017 & month == 7,.(Lon,Lat)]

# get Longitudes in [-180,180]
Lon_trafo = function(x){return((x+ 180) %% 360 -180)}

grid_GCFS1_2mt[,Lon := Lon_trafo(Lon)]

DT[,Lon := Lon_trafo(Lon)]

load(file = "../Derived/dt_combine_mr_wide.RData")

grid_NorCPM = dt[year == 1985 & month == 7,.(Lon,Lat)]

gm = contruct_grid_map(dt_ens = grid_GCFS1_2mt,dt_obs = grid_NorCPM)

setnames(gm,c("Lon2","Lat2","Lon","Lat"))

save(gm, file = "./grid_mapping_2mt.RData")

load( file = "./grid_mapping_2mt.RData")

### project on new grid

# easiest and fastest way is merging data tables:

# get all years in DT
ym = DT[Lon == Lon[1] & Lat == Lat[1],.(year,month)]
years = ym[,year]
months = ym[,month]


#parallelize

ym_ts = function(i)
{
  y = years[i]
  m = months[i]
  print(c(y,m))
  timeslize = merge(gm,DT[year == y & month == m,])
  return(timeslize)
}

DT_new = parallel::mclapply(X = 1:ym[,.N], FUN = ym_ts, mc.cores = 15)
DT_new = rbindlist(DT_new)

#correct names and order

DT_new[,Lon:= NULL]
DT_new[,Lat:= NULL]
setnames(DT_new,c("Lon2","Lat2"), c("Lon","Lat"))

DT_new = DT_new[order(year,month,Lon,Lat)]

DT = DT_new

DT[, Ens_bar := rowMeans(.SD),.SDcols = paste0("Ens",1:15)]
DT[, Ens_sd := rowMeans((.SD - Ens_bar)^2),.SDcols = paste0("Ens",1:15)]


save(DT,file = "./DT_2mt.RData")

for( i in 1:15)
{
  DT[,paste0("Ens",i):=NULL]  
}

save(DT,file = "./DT_2mt_mean.RData")


