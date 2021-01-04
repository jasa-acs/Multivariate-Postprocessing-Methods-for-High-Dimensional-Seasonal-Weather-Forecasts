#' A utility function that makes and stores the grid mapping
#'
#' @description This is essentially a wrapper for \code{construct_grid_map}
#'
#' @param out.file Location to store the output, which is then saved at \code{/dt_map.RData}
#'
#' @export
make_grid_mapping = function(out.file = "./Data/PostClim/SFE/Derived")
{
  dt_map = construct_grid_map()
  save(dt_map, file = paste0(out.file,"dt_map.RData"))
}


#' Make the combined dataset for all years and store in a wide format
#'
#' @param y_start The first year for which we have obs/ens data
#' @param y_stop The final year for which we have obs/ens data
#' @param vintage Which vintage are we using, often "mr" for most recent
#' @param data_dir The root directory that stores our data
#' @param grid_mapping_loc The location where the precomputed grid alignment object is stored.
#' @param output_loc The folder where the output should be stored. If \code{NULL} then put in the same location as \code{data_dir}
#' @param output_name The name of the output file.  If \code{NULL} then \code{paste0("dt_combine_",vintage,"_wide.RData")} will be used
#' @param lat_box The upper and lower bounds for the latitude.  By default the entire globe is used.
#' @param lon_box The upper and lower bounds for the longitude.  By default the entire globe is used.
#'
#' @export
make_combined_wide_dataset = function(y_start = 1985,
                                      y_stop = 2010,
                                      vintage = "mr",
                                      data_dir = "~/PostClimDataNoBackup/",
                                      grid_mapping_loc = "~/PostClimDataNoBackup/SFE/Derived/",
                                      output_loc = "~/PostClimDataNoBackup/SFE/Derived/",
                                      output_name = NULL,
                                      lat_box = c(-Inf,Inf),
                                      lon_box = c(-Inf,Inf))
{

  ##--- Error check
  if(lat_box[1] >= lat_box[2])
  {
    stop(paste("In make_combined_wide_data lat_box[1] is",
               round(lat_box[1],3),
               "while lat_box[2] is ", round(lat_box[2],3),
               "the order should be reversed"))
  }

  if(lon_box[1] >= lon_box[2])
  {
    stop(paste("In make_combined_wide_data lon_box[1] is",
               round(lon_box[1],3),
               "while lon_box[2] is ", round(lon_box[2],3),
               "the order should be reversed"))
  }


  
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
    for(m in 1:12)
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
  ##------------------------

  ##--------- Combine -----
  dt = rbindlist(dt_combine_all)
  dt[, YM := year * 12 + month]
  setkey(dt, "YM", "Lon", "Lat")
  ##------------------------

  ##---- Restrict ---
  dt = dt[ (Lon >= lon_box[1]) & (Lon <= lon_box[2]) & (Lat >= lat_box[1]) & (Lat <= lat_box[2])]
  
  ##----- Should I save or should I go? -----
  if(is.null(data_dir))
  {
    return(dt)
  }else{
    if(is.null(output_name))
    {
      output_name = paste0("dt_combine_",vintage,"_wide.RData")
    }
    f_name = paste0(output_loc,"/",output_name)
    save(dt,file = f_name)
    return(1)
  }
  ##-------------------------------------------

  
}

#' @export
make_NorCPM_wide_again_precip_2mtemp = function(save_dir = "~/PostClimDataNoBackup/SFE/Derived/",
                                                data.names = c("dt_combine_mr_wide_new.RData","dt_prect_NorCPM_wide.RData","dt_2mtemp_NorCPM_wide.RData")){
  
  data = fread("~/PostClimDataNoBackup/SFE/NorCPM_2mTemp/wide_data")
  
  old_names = paste0("wide_data_2.",c("ym",
                                      "lon",
                                      "lat",
                                      paste0("sst",1:10),
                                      "sst_bar",
                                      "sst_sd",
                                      paste0("ens",1:9),
                                      "ens_bar",
                                      "ens_sd",
                                      paste0("ts",1:9),
                                      "ts_bar",
                                      "ts_sd",
                                      paste0("prect",1:9),
                                      "prect_bar",
                                      "prect_sd",
                                      "sst_noaa",
                                      "year",
                                      "month"))
  
  new_names = c("YM",
                "Lon",
                "Lat",
                paste0("SST",1:10),
                "SST_bar",
                "SST_sd",
                paste0("Ens",1:9),
                "Ens_bar",
                "Ens_sd",
                paste0("ts",1:9),
                "ts_bar",
                "ts_sd",
                paste0("prect",1:9),
                "prect_bar",
                "prect_sd",
                "sst_noaa",
                "year",
                "month")
  
  
  setnames(data,old_names,new_names)
  
  # --- change "NULL" into NA and convert strings into doubles ---
  
  is.na(data) <- data == "NULL"
  
  num_names = c(paste0("SST",1:10),
                "SST_bar",
                "SST_sd",
                paste0("Ens",1:9),
                "Ens_bar",
                "Ens_sd",
                paste0("ts",1:9),
                "ts_bar",
                "ts_sd",
                paste0("prect",1:9),
                "prect_bar",
                "prect_sd",
                "sst_noaa")
  
  data[,(num_names) := lapply(.SD,as.numeric),.SDcols = num_names]
  
  ##------ include grid_id -------
  lon_all = sort(data[,unique(Lon)])
  n_lon = length(lon_all)
  cutoff_lon = c(-Inf,head(lon_all,-1) + diff(lon_all)/2)
  f_lon = approxfun(cutoff_lon, 1:n_lon, method="constant", rule = 2)
  
  lat_all = sort(data[,unique(Lat)])
  n_lat = length(lat_all)
  cutoff_lat = c(-Inf, head(lat_all,-1) + diff(lat_all)/2)
  f_lat = approxfun(cutoff_lat, 0:(n_lat - 1) * n_lon, method="constant", rule = 2)
  data[,grid_id := f_lon(Lon) + f_lat(Lat)]
  ##--------------------------------------
  
  #----break into smaller chunks and save
  
  dt = data[,.SD,.SDcols = c("Lon","Lat", paste0("SST",1:10), "SST_bar","SST_sd", paste0("Ens",1:9),"Ens_bar","Ens_sd","grid_id","year","month","YM","sst_noaa")]
  dt[,,.SDcols =  c(paste0("SST",1:10), "SST_bar","SST_sd","Ens_bar","Ens_sd")]
  save(dt, file = paste0(save_dir,data.names[1]))
  
  dt_prect = data[,.SD,.SDcols = c("Lon","Lat",paste0("prect",1:9),"prect_bar","prect_sd","grid_id","year","month","YM")]
  save(dt_prect, file = paste0(save_dir,data.names[2]))
  
  dt_2mtemp = data[,.SD,.SDcols = c("Lon","Lat",paste0("ts",1:9),"ts_bar","ts_sd","grid_id","year","month","YM")]
  save(dt_2mtemp, file = paste0(save_dir,data.names[3]))
  
}

