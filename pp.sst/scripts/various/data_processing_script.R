
####################################################################################################
########## Import data from netcdf files and store it in a data table of the required format #######
####################################################################################################

# This script shows how to construct a data table of the required type from the original data, provided as NetCDF files.

rm(list = ls())

library(ncdf4)
library(data.table)
library(pp.sst)

# The OISST dataset needs to be downloaded from ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2/sst.mnmean.nc

obs_dir = '~/PostClimDataNoBackup/SFE/OISST/' # the directory you stored the observations in

data = nc_open(paste0(obs_dir,'sst.mnmean.nc'))

# times and locations

times = as.Date(ncvar_get(data,varid = 'time'),origin = '1800-01-01')

lon = ncvar_get(data,varid = 'lon')
lat = ncvar_get(data,varid = 'lat')

# create data table

DT_obs = as.data.table(expand.grid(Lon = lon,Lat = lat,date = times))

DT_obs[,year:= year(date)][,month := month(date)]
DT_obs[,date:=NULL]

DT_obs[,Lon := (Lon <= 180) * Lon + (Lon > 180) * (Lon - 360)]

sst = ncvar_get(data,varid = 'sst')

DT_obs[,SST_bar := sst]

DT_obs = DT_obs[year >= 1985]

# check correct import:

plot_diagnostic(DT_obs,'SST_bar')

##########################################

# get an example forecast. The predictions can be downloaded at  ...

fc_dir = '~/PostClimDataNoBackup/SFE/NorCPM_Ocean/' # where did you store the data?

ens_ex = load_ensemble(year = 1985, month = 1,data_dir = fc_dir)


##########################################

# construct nearest-neighbour grid mapping between the observations grid and the NorCPM prediction grid:

data_dir = '~/PostClimDataNoBackup/SFE/NorCPM_Ocean/' # directory where the ensemble forecast is stored.

if(!file.exists(paste0(data_dir,'dt_map.RData')))
{
  dt_map = construct_grid_map(dt_ens = ens_ex,
                        dt_obs = DT)

  save(dt_map,file = paste0(data_dir,'dt_map.RData'))
} else {load(file = paste0(data_dir,'dt_map.RData'))}

setnames(dt_map,c('Lon_Obs','Lat_Obs','Lon_Ens','Lat_Ens'))

#################################

# construct the full data table

years = 1985:2010
months = 1:12

DT = data.table()

for(yy in years)
{
  for(mm in months)
  {
    print(paste0('year = ',yy,', month = ',mm))
    
    obs = DT_obs[year == yy & month == mm ]
    
    ens_fc = load_ensemble(year = yy, month = mm,data_dir = fc_dir)
    
    dt = combine_data_wide(dt_ens = ens_fc,dt_obs = obs,dt_map)
    
    DT = rbindlist(list(DT,dt))
  }
}

save(DT,file = paste0(fc_dir,'dt.RData'))
