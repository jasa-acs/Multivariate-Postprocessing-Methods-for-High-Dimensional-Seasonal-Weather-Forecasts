
#######################################################################
##  Estimate a monthly Exponential semi-variogram function.
#######################################################################




#' preparation for geostationary forecasts
#' 
#' @description Estimates the variogram for the residuals assuming anexponential covariance model and saves the results.
#' 
#' 
#' @param dt the data table.
#' @param training_years,m Integer vectors containing the year(s) and month(s) of the training dataset.
#' @param ens_size Integer. Size of the NWP ensemble.
#' @param saveorgo Logical, whether we save or not. 
#' @param save_dir,file_name The directory to save in. The name of the file of the saved variogram for a given month is \code{file_name <month> .RData}.
#' @param nintv Integer. How many distance bins are considered for the empirical variogram.
#' @param truncate Logical. The empirical variogram oftentimes is far from the fitted variogram for the 10 percent largest distances considered. If truncate == TRUE, those are ignored leading to a visually much better fit of the variogram.
#' 
#' @return data table containing n columns with noise and n columns with forecasts.
#' 
#' 
#' @author Claudio Heinrich, Yuan Qifen        
#' 
#' @importFrom spacetime STFDF
#' @importFrom sp SpatialPoints spDists SpatialPointsDataFrame
#' @importFrom gstat variogramST fit.variogram 
#' 
#' @export


geostationary_training = function (dt ,
                                   training_years = 1985:2000,
                                   m = 1:12,
                                   save_dir,
                                   file_name = "variogram_exp_m",
                                   nintv = 75,
                                   truncate = TRUE,
                                   mc_cores = 1){

  
  
  Lon_min  = dt[,range(Lon)][1]
  Lon_max  = dt[,range(Lon)][2]
  Lat_min  = dt[,range(Lat)][1]
  Lat_max  = dt[,range(Lat)][2]
  
  
  #parallelize by month
  
  gt_bm = function(mon)
  {
      print(paste0("month = ",mon))
      
      dt_temp = dt[month == mon & year %in% training_years,]
      
      land_ids <- which(dt_temp[, is.na(Ens_bar) | is.na(SST_bar)])
      if(!identical(land_ids,integer(0)))
      {
        dt_temp = dt_temp[-land_ids,]
      }
      
      sp <- sp::SpatialPoints(cbind(x=dt_temp[YM == min(YM), Lon],
                                    y=dt_temp[YM == min(YM), Lat]), 
                              proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
      
      # function to convert YM into the right date format:
      
      time_convert = function(YM){
        M = YM %% 12
        M[M == 0] = 12
        Y = (YM - M)/ 12
        M[M < 10]  = paste0(0,M[M < 10])
        as.Date(paste0(Y,"-",M,"-15"))
      }
      
      time = as.POSIXct( time_convert(unique(dt_temp[,YM])), tz = "GMT")
      
      setkey(dt_temp,YM,Lon,Lat) # for creating STFDFs the data should be ordered such that the 'spatial index is moving fastest'
      data = dt_temp[,.(trc(Ens_bar+Bias_Est)-SST_bar)] 
      setnames(data,"Res")
      
      stfdf = spacetime::STFDF(sp, time, data)
      
      #### Calculate the empirical semi-variogram
      #### ----------------------------------------------
      
      ## calculate the distance matrix [km], full symmetric matrix
      Dist <- sp::spDists(sp, longlat = TRUE)
      
      ## set the intervals
      up_Dist <- Dist[upper.tri(Dist, diag = FALSE)] 
      sort_up <- sort(up_Dist)
      
      bound_id = seq(1,length(up_Dist),length.out = nintv + 1)
      boundaries <- sort_up[bound_id]
      
      empVgm <- gstat::variogramST( Res ~ 1, stfdf, tlags=0, boundaries = boundaries, assumeRegular = TRUE, na.omit = TRUE) 
    
      #### Fit to the Exponential semi-variogram function
      #### -------------------------------------------------- 
    
      
      ## set the cutoff - sometimes performance improves if up to 10% of the largest distances are ignored in the variogram fit
      if(truncate) 
      {
        cutoff_ind  = ceiling(0.9*nintv)
        cutoff = empVgm$dist[cutoff_ind]
      } else {
        cutoff = max(empVgm$dist)
        }
      
    
      ## prepare empirical variograms for fitting
      spEmpVgm <- empVgm[empVgm$dist<=cutoff,] 
      spEmpVgm <- spEmpVgm[spEmpVgm$timelag==0,]
      sSpEmpVgm <- spEmpVgm[spEmpVgm$np!=0,] 
      spEmpVgm <- sSpEmpVgm[,1:3] 
      class(spEmpVgm) <- c("gstatVariogram", "data.frame")
      spEmpVgm$dir.hor <- 0
      spEmpVgm$dir.ver <- 0
      
      ## Exponential semi-variogram function with nugget, see the 1st argument. 
      ## Fixing "psill", fit "nugget" and "range", the 2nd arg.
      
      Mod <- gstat::fit.variogram(spEmpVgm, gstat::vgm("Exp"), fit.sills = c(T,T))
      
      # get and save covariance matrix
      
      psill <- Mod$psill[2]
      range <- Mod$range[2]
      nugget <- max(Mod$psill[1],0)
      
      Sigma <- psill*exp(-Dist/range)
      sills <- diag(Sigma) + nugget
      diag(Sigma) <- sills
      
      ns <- length(sp)
      
      svd_Sigma = svd(Sigma)
      
      sqrt_Sigma = svd_Sigma$u %*% sqrt(diag(svd_Sigma$d)) 
      
      save(sp,stfdf,Mod,Dist,sqrt_Sigma,file = paste0(save_dir,file_name,mon,".RData"))
  }
  
  if(mc_cores == 1)
  {
    for(mon in m)
    {
      gt_bm(mon)
    }
  } else {
    parallel::mclapply(X = m,FUN = gt_bm,mc.cores = mc_cores)
  }
  
}





#' geostationary forecasts 
#' 
#' @description Generates forecasts using a geostationary model with exponential covariance function for post-processing
#'      
#' 
#' @param dt the data table.
#' @param Y,M Integer vectors containing year(s) and month(s).
#' @param n Integer. Size of the desired forecast ensemble. 
#' @param noise Logical. Should the return data contain columns with noise?
#' @param var_dir,var_file_names Where the fitted variograms are stored.
#' 
#' @return data table containing n columns with forecasts.
#' 
#' 
#' @author Claudio Heinrich        
#' 
#' @importFrom MASS mvrnorm
#' 
#' @export

forecast_GS = function(dt, Y, M,
                       n=10,
                       mc_cores = 12,
                       noise = FALSE,
                       var_dir, var_file_name = "variogram_exp")
{

  dt = dt[year %in% Y,][month %in% M,]
  
  #find land grid ids:
  
  land_ids <- which(dt[, is.na(Ens_bar) | is.na(SST_bar)])
  if(!identical(land_ids,integer(0)))
  {
    dt_water = dt[-land_ids,]
  } else {
    dt_water = dt
  }
  
  SD_cols = c("Lon","Lat","grid_id","month","year","YM",
              "SST_hat","SST_bar","Ens_bar","Bias_Est","var_bar","SD_hat")
  SD_cols = SD_cols[which(SD_cols %in% colnames(dt))]
  
  fc_water <- na.omit( dt_water[,.SD,.SDcols = SD_cols])
  
  # parallelize:
  forecast_by_month = function(m)
  {
    print( paste0("Month = ",m))
    
    load(file = paste0(var_dir,var_file_name,"_m",m,".RData"))  
    
    psill <- Mod$psill[2]
    range <- Mod$range[2]
    nugget <- max(Mod$psill[1],0)
    
    Sigma <- psill*exp(-Dist/range)
    diag(Sigma) <- diag(Sigma) + nugget
    
    
    ns <- length(sp)
    
    dt_month = fc_water[month == m,]
    
    for (y in Y)
      {print(c(y,m))
      # marginally correct Sigma:
      mcf = dt_month[year == y, SD_hat]/sqrt(Sigma[1])
      sqrt_Sigma_hat =   diag(mcf) %*% sqrt_Sigma
      
      # generate noise:
      no <- sqrt_Sigma_hat %*% matrix(data = rnorm(n * dim(sqrt_Sigma_hat)[1]),ncol = n)
      for (i in 1:n)
        {
        dt_month[year == y, paste0("no",i) := no[,i]]
        dt_month[year == y, paste0("fc",i) := trc(Ens_bar + Bias_Est + .SD), 
                 .SDcols = paste0("no",i)]
        if(!noise)
        {
          dt_month[,paste0("no",i) := NULL]
        }
        
        }
      }
    return(dt_month)
  }
  
  
  temp = parallel::mclapply(M,forecast_by_month,mc.cores = mc_cores)
  fc_water = rbindlist(temp)
    
  #-------- add land --------------
  
  if(!identical(land_ids,integer(0)))
  {
    fc_land = dt[land_ids,]
    fc_land[,  paste0("fc",1:n):= NA]
    if(noise)
    {
      fc_land[,  paste0("no",1:n):= NA]
    }
   
    
    fc = rbindlist(list(fc_water,fc_land), fill = TRUE)
  }
  
  # order:
  
  fc = fc[ order(year,month,Lon,Lat)]
  
  return(fc)
}  
      
 