
#' Transforming the datatable (which has been univariately post-processed) by replacing temperatures by temperature anamolies
#' 
#' @description replaces all temperature-valued columns (forecasts and observations) in dt by anomalies w.r.t. the mean observed temp.
#'              over all climatology years. For the climatology years it leaves the current year out (out of sample).
#' 
#' @param dt the data table.
#' @param clim_years The years of climatology.
#' 
#' @return the transformed data table dt, containing an extra column \code{clim}.
#' 
#' @examples \dontrun{}
#' 
#' @author Claudio Heinrich        
#' 
#' @export

dt_transform_center = function(dt,clim_years,
                               trafo_colnames = c("SST_bar","SST_hat","Ens_bar",paste0("Ens",1:9)))
{
  Y = unique(dt[,year])
  M = unique(dt[,month])
  
  # get climatology parallelized
  
  clim_by_month = function(m)
  {
    dt_temp = list()
    for(y in Y)
    { 
      temp = dt[month ==m & year !=y & year %in% clim_years,mean(SST_bar),by = grid_id][,V1]
      
      if(dim(dt[month ==m & year ==y,])[1] != 0)
      {
        dt_temp_2 = dt[month == m & year == y,][, clim := temp]
        dt_temp = rbindlist(list(dt_temp,dt_temp_2)  )
      }
      
    }
    return(dt_temp)
  }
  
  temp = parallel::mclapply(M,clim_by_month,mc.cores = length(M))
  dt = rbindlist(temp)[order(year,month,Lon,Lat)]
  
  trafo_colnames = trafo_colnames[which(trafo_colnames %in% colnames (dt))]

  dt[,eval(trafo_colnames) := .SD - clim,.SDcols = trafo_colnames]  

  return(dt)  
}




#' Transforming the datatable (which has been univariately post-processed) by replacing temperatures by standardized temperature anamolies
#' 
#' @description replaces all temperature-valued columns (forecasts and observations) in dt by standardized anomalies w.r.t. 
#'              the mean and sd of theobserved temp. over all climatology years. For the climatology years it leaves the current
#'              year out (out of sample). Locations with (almost) 0 sd in climatology get assigned NAs.
#' 
#' @param dt the data table.
#' @param clim_years The years of climatology.
#' @param crit_ts Locations with climatology standard deviation below this threshold are transformed into NAs
#' 
#' @return the transformed data table dt, containing extra columns \code{clim} and \code{clim_sd}.
#' 
#' @examples \dontrun{}
#' 
#' @author Claudio Heinrich        
#' 
#' @export

dt_transform_stan = function(dt,clim_years,crit_ts = 0.01)
{
  Y = unique(dt[,year])
  M = unique(dt[,month])
  
  # get climatology and clim_sd, parallelized
  
  clim_by_month = function(m)
  {
    dt_month = list()
    for(y in Y)
    { 
      temp_dt = dt[month ==m & year !=y & year %in% clim_years,.(mean(SST_bar),sd(SST_bar)),by = grid_id]
      
      if(dim(dt[month ==m & year ==y,])[1] != 0)
      {
        dt_ym = dt[month == m & year == y,][, c("clim","clim_sd") := temp_dt[,.(V1,V2)]]
        dt_month = rbindlist(list(dt_month,dt_ym)  )
      }
      
    }
    return(dt_month)
  }
  
  temp = parallel::mclapply(M,clim_by_month,mc.cores = length(M))
  dt = rbindlist(temp)[order(year,month,Lon,Lat)]
  
  trafo_colnames = c("SST_bar","SST_hat","Ens_bar",paste0("Ens",1:9))
  trafo_colnames = trafo_colnames[which(trafo_colnames %in% colnames (dt))]
  
  dt[,eval(trafo_colnames) := (.SD - clim)/clim_sd,.SDcols = trafo_colnames]  
  
  # set the rows with 0 climatology sd to NAs:
  
  crit_rows = dt[,which(clim_sd < crit_ts)]
  dt[crit_rows,eval(trafo_colnames) := NA]  
  dt[crit_rows,c("clim","clim_sd") := NA]  
  
  
  return(dt)  
}

