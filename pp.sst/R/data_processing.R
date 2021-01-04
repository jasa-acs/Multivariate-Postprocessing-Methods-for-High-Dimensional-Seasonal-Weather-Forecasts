#' Create a mapping between the grid points on the
#' NorCPM ensemble and the SeNorge observational grid
#'
#' @description This helper function performs a nearest neighbor analysis to find the grid point in the NorCPM ensemble
#' that is closest to each point in the SeNorge dataset.
#' @param dt_ens A data.table of the ensemble members
#' @param dt_obs A data.table of the observation members
#'
#' @return A datatable which maps each lat/lon in the ensemble space to the nearest
#' member in the observation space
#'
#' @export
#' @importFrom geosphere distHaversine
#' @author Alex Lenkoski
#' @examples
#' ##construct_grid_map()
construct_grid_map = function(dt_ens = load_ensemble(1985,1),
                             dt_obs = load_observations(1985,1))
{

  dt_obs_grid = unique(dt_obs[,.(Lon,Lat)])
  dt_ens_grid = unique(dt_ens[,.(Lon,Lat)])
  point_match = NULL
  for(j in 1:dim(dt_obs_grid)[1])
  {
    if(j %% 1e2 == 0)print(paste0(j,"/",dim(dt_obs_grid)[1]))
    ##---------------------------------------
    a = geosphere::distHaversine(as.vector(dt_obs_grid[j,.(Lon,Lat)]),as.matrix(dt_ens_grid[,.(Lon,Lat)]))
    point_match[j] = order(a)[1]
    ##--------------------------------------
  }

  dt_map = data.table(dt_obs_grid[,.(Lon,Lat)], dt_ens_grid[point_match,.(Lon,Lat)])

  return(dt_map)
}

#' Load and postprocess a given NorCPM ensemble member
#'
#' @description This loads a given NorCPM ensemble member given the correct directory structure
#'
#' @param year Integer.  The year you are interested in.
#' @param month Integer.  The month
#' @param vintage One of Jan, Apr, Jul, Oct or mr.  If mr, the vintage closest to the given month is used.
#' @param data_dir String.  The root directory that stores the data.  Data are then assumed stored in \code{SFE/NorCPM_Ocean/}
#'
#' @examples
#' ##load_ensemble(1990,1)
#'
#' @export
load_ensemble = function(year,
                         month,
                         vintage = "mr",
                         data_dir = "~/PostClimDataNoBackup/SFE/NorCPM_Ocean/")
{

  ##------ Setup ----------
  month_num = month;  # store month as an integer
  if(month < 10)month = paste0("0",month)
  
  if(vintage == "Jan")vin_mon = 1
  if(vintage == "Apr")vin_mon = 4
  if(vintage == "Jul")vin_mon = 7
  if(vintage == "Oct")vin_mon = 10
  ##-----------------------

  ##---- Collect Grid Info ------------------
  gridname <- paste0(data_dir,"grid.nc")
  ncgrid <- ncdf4::nc_open(gridname)
  grid_lon_ens <- ncdf4::ncvar_get(ncgrid, "plon")
  grid_lat_ens <- ncdf4::ncvar_get(ncgrid, "plat")
  ncdf4::nc_close(ncgrid)
  ##-----------------------------------------

  ##-------- Find target run ----------------
  filedir <- data_dir
  ff_all = system(paste0("ls ",filedir,"*_mem01.micom.hm.",year,"-",month,".nc"), intern = TRUE)
  
  if (vintage == "mr") {ff_use = tail(ff_all,1)
  } else if (vintage == "2r") {ff_use = rev(ff_all)[2]
  } else if (vintage == "3r") {ff_use = rev(ff_all)[3]
  } else if (vintage == "4r") {ff_use = rev(ff_all)[4]
  } else  {month_num <- (month_num-vin_mon+1)%%12 # start the month labelling at the vintage month
    if(month_num == 0)month_num == 12
    if (month_num < 4) {ff_use = rev(ff_all)[1]
    }else if (month_num < 7) {ff_use = rev(ff_all)[2]
    }else if (month_num < 10) {ff_use = rev(ff_all)[3]
    }else  ff_use = rev(ff_all)[4]
    }
  ##-----------------------------------------
  
  ##----- Get Raw Ensemble ------------------
  sst_ensemble = array(NA,dim = c(dim(grid_lon_ens)[1], dim(grid_lat_ens)[2],9))
  for(j in 1:9)
  {
    ff_new = gsub("mem01",paste0("mem0",j),ff_use)
    ncin <- ncdf4::nc_open(ff_new)
    sst_ensemble[,,j] <- ncdf4::ncvar_get(ncin, "sst")
    ncdf4::nc_close(ncin)
  }
  ##-----------------------------------

  ## --- Make a Data Table ----------------
  dt_ensemble_list = list()
  for(j in 1:9)
  {
    dt_ensemble_list[[j]] = data.table(Lon = as.vector(grid_lon_ens),
                                       Lat = as.vector(grid_lat_ens),
                                       Ens = j,
                                       Forecast = as.vector(sst_ensemble[,,j]))
  }
  dt_ensemble = rbindlist(dt_ensemble_list)
  setkey(dt_ensemble,"Lon","Lat","Ens")
  ##-----------------------------------------

  return(dt_ensemble)
}

#' Load and postprocess a given SeNorge observational member
#'
#' @description This loads a given SeNorge observational member given the correct directory structure
#'
#' @param year Integer.  The year you are interested in.
#' @param month Integer.  The month
#' @param obsdir String.  The root directory that stores the data.
#'
#' @examples
#' ##load_observations(1990,1)
#'
#' @export
load_observations = function(year, month,
                             obsdir = "./Data/HadiSST2/")
{
  ##------ Load Observations -------------------
  if(month < 10)month = paste0("0",month)
  obsname = paste0("SST_ens_", year,"_", month,".nc")
  ncnameobs = paste0(obsdir, obsname)
  ncobs = ncdf4::nc_open(ncnameobs)
  ##--------------------------------------------
  
  ##------- Extract Observations ---------------
  sst_obs = ncdf4::ncvar_get(ncobs, "sst")
  sst_obs = sst_obs - 273.15 ## Convert from Kelvin
  lon_obs = ncdf4::ncvar_get(ncobs, "longitude")
  lat_obs = ncdf4::ncvar_get(ncobs, "latitude")
  n_lon = length(lon_obs)
  n_lat = length(lat_obs)
  grid_lon_obs = matrix(rep(lon_obs,length(lat_obs)),n_lon,n_lat)
  grid_lat_obs = matrix(rep(lat_obs,length(lon_obs)),n_lon,n_lat, byrow=TRUE)
  dt_obs_list = list()
  for(j in 1:dim(sst_obs)[3])
  {
    dt_obs_list[[j]] = data.table(Lon = as.vector(grid_lon_obs),
                                  Lat = as.vector(grid_lat_obs),
                                  Obs = j,
                                  SST = as.vector(sst_obs[,,j]))
  }
  dt_obs = rbindlist(dt_obs_list)
  setkey(dt_obs, "Lon", "Lat")
  ##-------------------------------------------

  ncdf4::nc_close(ncobs)
  
  return(dt_obs)
}

#' Align and combine an ensemble (NorCPM) and an observation (SeNorge) set of data and store this in a wide format
#'
#' @description This aligns a collection of ensemble members with the associated observational member to create a single dataset.  Data are stored in "wide" format, indicating that all 9 ensemble members and all 10 observational outcomes for a given y,m,grid location are stored in one line of data with a whole bunch of columns.
#'
#' @param dt_ens The data.table of ensemble memebers
#' @param dt_obs The data.table of observational members
#' @param dt_map A data.table that maps the two between one another
#'
#' @export
#' @author Alex Lenkoski
#' @examples
#' ##dt_ens = load_ensemble(1990,1)
#' ##dt_obs = load_observations(1990,1)
#' ##dt_map = load_mapping()
#' ##combine_data_wide(dt_ens,dt_obs,dt_map)
combine_data_wide= function(dt_ens, dt_obs, dt_map)
{
  
  ##------- Collapse Ensemble -----------
  dt_ens[,Ens_Name:=paste0("Ens",Ens)]
  dt_ens_wide= dcast(dt_ens, Lon + Lat ~ Ens_Name, value.var = "Forecast")
  dt_ens_wide[,Ens_bar:=rowMeans(dt_ens_wide[,paste0("Ens",1:9)])]
  dt_ens_wide[,Ens_sd:= apply(dt_ens_wide[,paste0("Ens",1:9)],1,"sd")]
  ##-------------------------------------
  
  ##------- First fill out ensemble with obs lookup ---
  setkey(dt_map, "Lon_Ens","Lat_Ens")
  setkey(dt_ens_wide, "Lon", "Lat")
  dt_ens_wide = dt_ens_wide[dt_map]
  ##-----------------------------------------------
  
  ##------- Now fill in Observations --------
  setkey(dt_ens_wide, "Lon_Obs", "Lat_Obs")
  setkey(dt_obs,'Lon','Lat')
  dt_combine = dt_obs[dt_ens_wide]
  dt_combine[,i.Lon := NULL]
  dt_combine[,i.Lat := NULL]
  ##-----------------------------------------
  
  ##------ Form a convenient key -------
  lon_all = sort(dt_combine[,unique(Lon)])
  n_lon = length(lon_all)
  cutoff_lon = c(-Inf,head(lon_all,-1) + diff(lon_all)/2)
  f_lon = approxfun(cutoff_lon, 1:n_lon, method="constant", rule = 2)
  
  lat_all = sort(dt_combine[,unique(Lat)])
  n_lat = length(lat_all)
  cutoff_lat = c(-Inf, head(lat_all,-1) + diff(lat_all)/2)
  f_lat = approxfun(cutoff_lat, 0:(n_lat - 1) * n_lon, method="constant", rule = 2)
  dt_combine[,grid_id := f_lon(Lon) + f_lat(Lat)]
  ##--------------------------------------
  return(dt_combine)
}

#'  Load a pre-rendered combined wide dataset
#'
#' @description Given a pre-rendered dataset, this loads it from the correct location
#'
#' @param data_dir Root location of the data
#' @param vintage Indication of which vintage we are using, most likely 'mr'
#' @param bias Boolean. Should the bias corrected ensemble be loaded?
#' @param var Boolean. Should the data with estimated bias and variance be loaded? 
#' @param model String.  Currently taking values \code{NorESM} or \code{senorge} to indicate whether the NorESM or senorge data are of interest
#'
#' @return A data.table containing all data for the parameters specified, provided these data have been rendered.
#' @export
load_combined_wide = function(data_dir = "~/PostClimDataNoBackup/SFE/Derived/", 
                              vintage = "mr", 
                              bias = FALSE,
                              var = FALSE,
                              output_name = NULL
                              )
{
  if(is.null(output_name))
    {
    if(var)
    {
      file = paste0(data_dir,"/dt_combine_wide_bc_var.RData")
    }else if(bias)
      {
        file = paste0(data_dir,"/dt_combine_wide_bias.RData")
      } else {
        file = paste0(data_dir,"/dt_combine_",vintage,"_wide.RData")
      }
    }else{
      file = paste0(data_dir,"/",output_name)
    }
  load(file)
  return(dt)
}
