
#' computes simple moving averages, resistant to missing values  
#'
#' @param l length of the averaging window.
#' @param vec,years vectors of the same length, vec[i] contains the value corresponding to year years[i]
#' @skip Integer. If skip = n > 0 the moving average skips the most recent n years. Useful if realizations in the most recent past are missing.
#'                                    
#' @return a vector of the same length as the input vectors. At location j it contains the average of the values contained in vec that fall into the period of the last l years (excluding the present). First entry is 0. NAs are ignored.
#'
#' @author Claudio Heinrich
#' @examples sim_mov_av(5, rnorm(10), c(1990,1993,1995:2002))
#' 
#' 
#' @export

sim_mov_av = function(l,vec, years, skip = 0, twosided = FALSE){
  
  all_years = min(years):max(years)
  
  #set NAs in vec to 0 in order to ignore them
  vec_na_rm = vec
  na_loc = which(is.na(vec))
  if(!identical(na_loc,integer(0)))
  {
    vec_na_rm[na_loc] = 0
  }
  
  sma = rep(0,length(vec))
  
  if(!twosided)
  {
    for (i in 1:length(vec)){
      y=years[i]
      weight_vec = rep(0,length(vec))
      rel_loc = (y > years) & ((y - years) <= l) & (!is.na(vec))
      skip_loc =  y-years <= skip
      weight_vec = rel_loc & !skip_loc
      if(TRUE %in% weight_vec) weight_vec = weight_vec/sum(weight_vec)
      sma[i] = sum(weight_vec * vec_na_rm)
    }
  } else if(twosided)
  {
    for (i in 1:length(vec)){
      y=years[i]
      weight_vec = rep(0,length(vec))
      rel_loc = abs(years - y) <= l & (!is.na(vec))
      skip_loc = abs(years - y) <= skip
      weight_vec = rel_loc & !skip_loc
      if(TRUE %in% weight_vec) weight_vec = weight_vec/sum(weight_vec)
      sma[i] = sum(weight_vec * vec_na_rm)
    }
  }
  
  return(sma)  
}

#' computes exponentially weighted moving averages, is resistant to missing values  
#'
#' @param a weight parameter.
#' @param vec,years vectors of the same length, vec[i] contains the value corresponding to year years[i]
#'                  
#' @return a vector of the same length as the input vectors. At location j it contains the average of the past entries of vec, weighted by exp(-a*d) where d is the distance to the current year. First entry is 0.
#' @skip Integer. If skip = n > 0 the moving average skips the most recent n years. Useful if realizations in the most recent past are missing.
#'
#' @author Claudio Heinrich
#' @examples exp_mov_av(.1, rnorm(10), c(1990,1993,1995:2002))
#' 
#' 
#' @export


exp_mov_av = function( a,vec, years, skip = 0,twosided = FALSE){
  
  lv = length(vec)
  
  all_years = min(years):max(years)
  
  #set NAs in vec to 0 in order to ignore them
  vec_na_rm = vec
  na_loc = which(is.na(vec))
  if(!identical(na_loc,integer(0)))
  {
    vec_na_rm[na_loc] = 0
  }
  
  exp_weights = exp(-a * (0:(length(all_years)-1)))
  
  ema = rep(0,lv)
  
  if(lv > 1)
  {
    if(!twosided)
    {
      for (i in (skip + 2):lv){
        weight_vec = rev(exp_weights[which(all_years %in% years[1:(i-1-skip)])])
        #normalize
        if(!identical(na_loc,integer(0)))
        {
          weight_vec = weight_vec/sum(weight_vec[-na_loc])
        } else {
          weight_vec = weight_vec/sum(weight_vec)
        }
        ema[i] = sum(weight_vec*vec_na_rm[1:(i-1-skip)])
      }
    } else if(twosided)
    {
      for (i in 1:lv){
        year_dist = abs(years[1 : lv] - years[i])
        weight_vec = rep(0,lv)
        weight_vec[!(year_dist <= skip)] = exp_weights[year_dist]
        #normalize
        if(!identical(na_loc,integer(0)))
        {
          weight_vec = weight_vec/sum(weight_vec[-na_loc])
        } else {
          weight_vec = weight_vec/sum(weight_vec)
        }
        ema[i] = sum(weight_vec*vec_na_rm)
      }
    }
  }
  
  return(ema)  
}



#' crps score for normal distribution, resistant to missing values and 0 standard deviation  
#'
#' @param y vector of observations.
#' @param mean,sd mean and sd of the forecasting distribution
#'                   
#' @return vector of the same length as y containing crps scores.
#'
#' @author Claudio Heinrich
#' 
#' @importFrom scoringRules crps
#' 
#' @export


crps_na_rm = function(y, mean, sd){
  
  na_loc = which( is.na(y) | is.na(mean) | is.na(sd) | sd == 0)
  
  if(identical(na_loc,integer(0))){
    x = scoringRules::crps(y, family = "normal", mean = mean, sd = sd)
  }else{
    x = rep(0,length(y))
    x[-na_loc] = scoringRules::crps(y[-na_loc], family = "normal", mean = mean[-na_loc], sd = sd[-na_loc])
    x[na_loc] = NA
  }
  return(x)
}

#' Computes scores for bias correction and variance estimation
#'
#' @param DT The data table
#' @param ens_size Size of the ensemble.
#' @param eval_years The years to compute the scores for
#' @param var Logical. If TRUE the CRPS is computed, else the MSE of the bias corrected forecast is computed
#'                   
#' @return A one-row data table containing the mean score.
#'
#' @author Claudio Heinrich
#' @examples \dontrun{
#' save_dir = "~/PostClimDataNoBackup/SFE/Derived/NAO/"
#' DT = load_combined_wide(data_dir = save_dir, output_name = "dt_combine_wide_bc_var.RData")
#' global_mean_scores(DT)}
#' 
#' 
#' @export


global_mean_scores_EnsBMA = function (DT, ens_size = 9, eval_years = 2001:2010, var = TRUE)
{
  if(var){
    # marginal distribution is a mixture distribution of normal distributions centered at the ensemble members + Bias_Est, all with variance SD_hat^2.
    
    na_ids = DT[year %in% eval_years,which(is.na(Ens_bar) | is.na(Bias_Est) | is.na(SD_hat))]
    DT_temp = DT[year %in% eval_years,][-na_ids,]
    
    obs = DT_temp[,SST_bar]
    means = as.matrix(DT_temp[,trc(.SD + Bias_Est),.SDcols = paste0("Ens",1:ens_size)])
    sds = rep(DT_temp[,SD_hat], ens_size)
    sds = matrix(sds, ncol = ens_size)
    # get CRPS for each location for Gaussian mixture distribution
    mar_CRPS = scoringRules::crps_mixnorm(obs,means,sds)
    glob_mean_sc = data.table("CRPS" = mean(mar_CRPS))
    
  } else {
    glob_mean_sc = DT[year %in% eval_years, 
                      .("MSE" = mean( (SST_bar - SST_hat)^2, na.rm=TRUE))]
  }
  
  return(glob_mean_sc)
}


global_mean_scores = function (DT, eval_years = 2001:2010, var = TRUE){
  
  if(var){
    glob_mean_sc = DT[year %in% eval_years, 
                    .("CRPS" = mean(crps_na_rm(SST_bar, mean = SST_hat,sd = SD_hat),na.rm = TRUE))]
  } else glob_mean_sc = DT[year %in% eval_years, 
                           .("MSE" = mean( (SST_bar - SST_hat)^2, na.rm=TRUE))]
  
  return(glob_mean_sc)
}

#' sets up the estimation of marginal variance by computing the ensemble variance by location, month and year.
#'
#' @param dt The data table containing the estimated bias.
#' @param ens_size integer. Size of the ensemble.
#' @param saveorgo logical. Do we save the data table?
#' @param save_dir,file_name character strings. Directory and file name for saving.
#' @param mean_est For internal use - don't change this.
#' 
#' @return the data table dt with a new column labelled 'var_bar'
#' 
#' @author Claudio Heinrich
#' 
#' @examples 
#' \dontrun{load_combined_wide()
#'          ens_sd_est(dt)}
#'          
#' @export

ens_sd_est = function(dt, 
                      ens_size = 9,
                      saveorgo = TRUE,
                      mean_est = "sv",
                      save_dir = "~/PostClimDataNoBackup/SFE/Derived/",
                      file_name = "dt_combine_wide_bias_sd.RData")
{# get the average variation of the bias corrected truncated ensemble members around the true value by year, month and grid_id:
  if(mean_est == "bcf"){
    get_variances = function(i)
    {
      var_setup_dt = as.data.table(dt[, (SST_hat - SST_bar)^2/ens_size,.SDcols = paste0("Ens",i)])
      #var_setup_dt = as.data.table(var_setup_dt)
      return(var_setup_dt)
    }
    
    var_list = parallel::mclapply(1:ens_size,get_variances,mc.cores = ens_size)
    var_list = as.data.table(var_list)
    setnames(var_list,paste0("temp",1:ens_size))
  
    dt[,"var_bar" := rowSums(var_list)]
      
  } else if (mean_est == "sm"){
    get_variances = function(i)
    {
    var_setup_dt = as.data.table(dt[, (.SD - Ens_bar)^2/(ens_size-1),.SDcols = paste0("Ens",i)])
    #var_setup_dt = as.data.table(var_setup_dt)
    return(var_setup_dt)
    }
    
    var_list = parallel::mclapply(1:ens_size,get_variances,mc.cores = ens_size)
    var_list = as.data.table(var_list)
    setnames(var_list,paste0("temp",1:ens_size))
  
    dt[,"var_bar" := rowSums(var_list)]
  } else if (mean_est == "sv") 
  {
    dt[,"var_bar" := (SST_bar-SST_hat)^2]
  }
  
  
return(dt)
}



#' Applies bias correction with a specified method to the data and saves or returns scores.
#'
#' @param dt The data table.
#' @param method Method of bias correction. Takes "sma" for simple moving average and "ema" for exponential moving average. 
#' @param par_1 Numeric. If method == "sma", par_1 is the (integer) length of the moving average window, if method == "ema", par_1 is the scale parameter for the exponential downscaling, typically in (0,1).
#' @param ... arguments passed on to the moving average function, e.g. two_sided for the training period.                   
#'                   
#' @return A data table containing columns .(year,month,Lon,Lat,Bias_est,SST_hat).
#'
#' @author Claudio Heinrich
#' 
#' @export

bias_correct = function(dt, method, par_1, ...)
{
  if(method == "sma"){
    b_hat =  dt[,.(year,Bias_est = sim_mov_av(l = par_1,  
                                              vec = SST_bar - Ens_bar, 
                                              years = year,
                                              ...)),
                by = .(Lon,Lat,month)]
  }
  if (method == "ema"){
    b_hat =  dt[,.(year,Bias_est = exp_mov_av(a = par_1,  
                                              vec = SST_bar - Ens_bar, 
                                              years = year,
                                              ...)),
                by = .(Lon,Lat,month)]
  }
  
  setkey(b_hat,year,month,Lon,Lat)
  
  # estimated temperature
  ret_val = data.table(b_hat, SST_hat = trc(dt[,Ens_bar] + b_hat[,Bias_est]))
  
  return(ret_val) 
  
}



#' computes the RMSE when instead of bias correction the linear model SST_hat = a + b Ens_bar is used (as in standard NGR), with data grouped by month.
#' 
#' @param DT the data table.
#' @param months The considered months.
#' @param training_years,validation_years Training years and validation years.
#' 
#'
#' @examples \dontrun{ DT = load_combined_wide()
#'                     bias_lr_bm(DT = DT)}
#'   
#' @author Claudio Heinrich
#' 
#' @export

bias_lr_bm = function(DT,
                      months = 1:12,
                      validation_years = 2001:2010)
{
  
  for(y in validation_years) # parallelizing runs into memory issues, and for some reason is not faster?
  {
    print(y)
    # grouped by month
    
    fits_by_month = lme4::lmList(formula = SST_bar ~ 1 + Ens_bar | month,
                                 data=DT[year < y,.(SST_bar,Ens_bar,month)])
    
    months = as.character(DT[year == y,month])
    DT[year == y, c("a_bm","b_bm"):= coef(fits_by_month)[months,]]
    DT[year == y, T_hat_lr_bm := a_bm + b_bm * Ens_bar]
  }
  return(DT)
}


#' computes climatology and sample climatology, always oos starting from the last year
#' 
#' @param DT the data table.

#' @author Claudio Heinrich
#' 
#' @export

compute_clim = function(DT)
{
  DT[,clim := (cumsum(SST_bar) - SST_bar)/(year - min(year) ),.(month,Lon,Lat)]
  DT[,fc_clim := (cumsum(Ens_bar) - Ens_bar)/(year - min(year)),.(month,Lon,Lat)]
  DT[,ano:= SST_bar - clim]
  DT[,fc_ano:= SST_hat - fc_clim]
  
  return(DT)
}


#' computes the MSE when instead of bias correction the locally adaptive model SST_ano = a + b Ens_bar_ano is used (as in standard NGR), with data grouped by month.
#' 
#' @param DT the data table.
#' @param months The considered months.
#' @param validation_years validation years.
#' 
#'
#' @examples \dontrun{ DT = load_combined_wide()
#'                     bias_lr_bm(DT = DT)}
#'   
#' @author Claudio Heinrich
#' 
#' @export

bias_lr_bm_la = function(DT,
                      months = 1:12,
                      validation_years = 2001:2010)
{
  
  for(y in validation_years) # parallelizing runs into memory issues, and for some reason is not faster?
  {
    print(y)
    # grouped by month
    
        fits_by_month = lme4::lmList(formula = ano ~ 1 + fc_ano | month,
                                 data=DT[year < y,.(ano,fc_ano,month)])
    
    months = as.character(DT[year == y,month])
    DT[year == y, c("a_bm_la","b_bm_la"):= coef(fits_by_month)[months,]]
    DT[year == y, T_hat_lr_bm_la := clim + a_bm_la + b_bm_la * fc_ano]
  }
  return(DT)
}

#' computes the RMSE when instead of bias correction the linear model SST_hat = a + b Ens_bar is used (as in standard NGR), with data grouped by location.
#' 
#' @param DT the data table.
#' @param months The considered months.
#' @param training_years,validation_years Training years and validation years.
#' 
#'
#' @examples \dontrun{ DT = load_combined_wide()
#'                     bias_lr_bl(DT = DT)}
#'   
#' @author Claudio Heinrich
#' 
#' @export


bias_lr_bl = function(DT,
                      months = 1:12,
                      validation_years = 2001:2010)
{
  
  # grouped by location
  for(y in validation_years)
  {
    print(y)
    
    fits_by_loc = suppressWarnings(lme4::lmList(formula = SST_bar ~ 1 + Ens_bar | grid_id,
                                                data=DT[year < y,.(SST_bar,Ens_bar,grid_id)]))
    
    grid_ids = as.character(DT[year == y,grid_id])
    DT[year == y,c("a_bl","b_bl"):= coef(fits_by_loc)[grid_ids,]]
    DT[year == y,T_hat_lr_bl := a_bl + b_bl * Ens_bar]
  }
  return(DT)
}

#' computes the RMSE when instead of bias correction the linear model SST_hat = a + b Ens_bar is used (as in standard NGR), with data grouped by location.
#' 
#' @param DT the data table.
#' @param months The considered months.
#' @param training_years,validation_years Training years and validation years.
#' 
#'
#' @examples \dontrun{ DT = load_combined_wide()
#'                     bias_lr_bl(DT = DT)}
#'   
#' @author Claudio Heinrich
#' 
#' @export


bias_lr_bl_la = function(DT,
                      months = 1:12,
                      validation_years = 2001:2010)
{
  
  # grouped by location
  for(y in validation_years)
  {
    print(y)
    
    fits_by_loc = suppressWarnings( lme4::lmList(formula =  ano ~ 1 + fc_ano | grid_id,
                               data=DT[year < y,.(ano,clim,fc_ano,grid_id)]))
    
    grid_ids = as.character(DT[year == y,grid_id])
    DT[year == y,c("a_bl_la","b_bl_la"):= coef(fits_by_loc)[grid_ids,]]
    DT[year == y,T_hat_lr_bl_la := clim + a_bl_la + b_bl_la * fc_ano]
  }
  return(DT)
}

#' computes the RMSE when instead of bias correction the linear model SST_hat = a + b Ens_bar is used (as in standard NGR), with data grouped by both month and location.
#' 
#' @param DT the data table.
#' @param months The considered months.
#' @param training_years,validation_years Training years and validation years.
#' 
#'
#' @examples \dontrun{ DT = load_combined_wide()
#'                     bias_lr_bm(DT = DT)}
#'   
#' @author Claudio Heinrich
#' 
#' @export

bias_lr_bb = function(DT,
                      months = 1:12,
                      validation_years = 2001:2010)
{
  
  # grouped by location
  for(y in validation_years)
  {
    for(m in months)
    {
      print(c(y,m))
      
      data = DT[month == m][year < y,.(SST_bar,Ens_bar,grid_id)]
      data[,grid_id := as.factor(grid_id)]
      fits_by_both = suppressWarnings(lme4::lmList(formula = SST_bar ~ 1 + Ens_bar | grid_id,
                                                   data=data))
      gids = as.character(DT[month == m][year == y,grid_id])
      DT[month == m & year == y,c("a_bb","b_bb"):= coef(fits_by_both)[gids,]]
      DT[month == m & year == y,T_hat_lr_bb := a_bb + b_bb * Ens_bar]
      
    }
    
  }
  return(DT)
}


#' computes the RMSE when instead of bias correction the linear model SST_hat = a + b Ens_bar is used (as in standard NGR), with data grouped by both month and location.
#' 
#' @param DT the data table.
#' @param months The considered months.
#' @param training_years,validation_years Training years and validation years.
#' 
#'
#' @examples \dontrun{ DT = load_combined_wide()
#'                     bias_lr_bm(DT = DT)}
#'   
#' @author Claudio Heinrich
#' 
#' @export

bias_lr_bb_la = function(DT,
                         months = 1:12,
                         validation_years = 2001:2010)
{
  
  # grouped by location
  for(y in validation_years)
  {
    for(m in months)
    {
      print(c(y,m))
      
      data = DT[month == m][year < y,.(clim,ano,fc_ano,grid_id)]
      data[,grid_id := as.factor(grid_id)]
      fits_by_both = suppressWarnings(lme4::lmList(formula = ano ~ 1 + fc_ano | grid_id,
                                  data=data))
      gids = as.character(DT[month == m][year == y,grid_id])
      DT[month == m & year == y,c("a_bb_la","b_bb_la"):= coef(fits_by_both)[gids,]]
      DT[month == m & year == y,T_hat_lr_bb_la := clim + a_bb_la + b_bb_la * fc_ano]
      
    }
    
  }
  return(DT)
}


#' parallelized version of bias_lr_bb
#' 
#' @param DT the data table.
#' @param months The considered months.
#' @param validation_years Validation years.
#' @param ... Arguments passed on to mclapply, most importantly mc.cores.
#' 
#'
#' @examples \dontrun{ DT = load_combined_wide()
#'                     bias_lr_bm(DT = DT)}
#'   
#' @author Claudio Heinrich
#' 
#' @export

bias_lr_bb_par = function(DT,
                          months = 1:12,
                          validation_years = 2001:2016,
                          mc_cores
                          )
{
  if('a' %in% colnames(DT))
  {
    DT[,c('a','b','T_hat_lr_both'):= NULL]
  }
  
  gids = DT[,unique(grid_id)]
  gids_char = as.character(gids)
  
  bias_by_month_par = function(m)
  {
    ret_dt = as.data.table(expand.grid(month = m, year = validation_years,grid_id = gids))
    for(y in validation_years)
    {
      data = DT[month == m][year < y,.(SST_bar,Ens_bar,grid_id)]
      data[,grid_id := as.factor(grid_id)]
      fits_by_both = lme4::lmList(formula = SST_bar ~ 1 + Ens_bar | grid_id,
                                  data=data)
      
      ret_dt[month == m & year == y,c("a","b"):= coef(fits_by_both)[gids_char,]]
      ret_dt[month == m & year == y,T_hat_lr_both := a + b * DT[month == m & year ==y,Ens_bar]]
    }
    return(ret_dt)
  }
  
  T_hat = rbindlist(parallel::mclapply(months,bias_by_month_par, mc.cores = mc_cores))
  
  setkey(T_hat,year,month,grid_id)
  DT = merge(DT,T_hat,all.x = TRUE)
  
  return(DT)
}





#' Estimates standard deviation with a specified method to the data and saves or returns scores.
#'
#' @param dt The data table.
#' @param method Method of variance estimation. Takes "sma" for simple moving average and "ema" for exponential moving average. 
#' @param par_1 Numeric. If method == "sma", par_1 is the (integer) length of the moving average window, if method == "ema", par_1 is the scale parameter for the exponential downscaling, typically in (0,1).
#' @param scores Logical. If true, the CRPS is returned.
#' @param eval_years Numerical vector. The years for evaluating the score.
#' @param saveorgo Logical. If TRUE, the data table with new column SD_hat is saved.
#' @param save_dir,file_name Directory and name for the saved file.
#'                   
#' @return The data table dt containing a new column SD_hat.
#'
#' @author Claudio Heinrich
#' @examples \dontrun{sd_est(saveorgo = FALSE)}
#' 
#' 
#' @export

sd_est = function(dt ,
                  method = "sma", # "sma" for 'simple moving average',
                        # "ema" for 'exponential moving average'
                  par_1 = 16,   # if method == sma, par_1 is the length of window for the sma
                        # if method == ema, par_1 is the ratio of the exp. mov. av.
                  scores = FALSE,
                  eval_years = 2001:2010,
                  saveorgo = TRUE,
                  save_dir = "~/PostClimDataNoBackup/SFE/Derived/",
                  file_name = "dt_combine_wide_bc_var.RData")
{

  if(method == "sma"){
    dt[year != min(year),"SD_hat" := sqrt(sim_mov_av(l = par_1, 
                                                     vec = var_bar, 
                                                     years = year)),
       by = .(grid_id, month)]
  }
  if (method == "ema"){
    dt = dt[year != min(year),"SD_hat" := sqrt(exp_mov_av(a = par_1,
                              vec = var_bar, 
                              years = year)),
            by = .(grid_id, month)]
  }
  
  if(saveorgo){
    save(dt, file = paste0(save_dir,file_name))
  }
  
  if(scores){
    mean_sc = global_mean_scores(dt, eval_years = eval_years, var = TRUE)
    return(mean_sc)
  } else return(dt)
  
}



#' Estimates standard deviation with a specified method to the data and saves or returns scores.
#'
#' @param dt The data table.
#' @param method Method of variance estimation. Takes "sma" for simple moving average and "ema" for exponential moving average. 
#' @param par_1 Numeric. If method == "sma", par_1 is the (integer) length of the moving average window, if method == "ema", par_1 is the scale parameter for the exponential downscaling, typically in (0,1).
#' @param scores Logical. If true, the CRPS is returned.
#' @param eval_years Numerical vector. The years for evaluating the score.
#' @param saveorgo Logical. If TRUE, the data table with new column SD_hat is saved.
#' @param save_dir,file_name Directory and name for the saved file.
#'                   
#' @return The data table dt containing a new column SD_hat.
#'
#' @author Claudio Heinrich
#' @examples \dontrun{sd_est(saveorgo = FALSE)}
#' 
#' 
#' @export

sd_est_2 = function(dt, method, par_1,...)
{
  
  if(method == "sma"){
    sd_hat = dt[,"SD_hat" := sim_mov_av(l = par_1, 
                                      vec = var_bar, 
                                      years = year,
                                      ...),
            by = .(Lon,Lat, month)]
  }
  if (method == "ema"){
    sd_hat = dt[,"SD_hat" := exp_mov_av(a = par_1,
                                 vec = var_bar,
                                 years = year,
                                 ...),
       by = .(Lon,Lat, month)]
  }
  
  sd_hat[,SD_hat := sqrt(SD_hat)]
  sd_hat = sd_hat[,.(year,month,Lon,Lat,SD_hat)]
  
  return(sd_hat)

}



#' Estimates the variance as has been suggested for NGR, grouped by month
#' 
#' @param dt the data table.
#' @param months The considered months.
#' @param validation_years Training years and validation years.
#' 
#' @return data table containing the RMSEs for the model above, where the coefficients are estimated in three different ways: grouped by month, location and by both.
#'
#'   
#' @author Claudio Heinrich
#' 
#' @export

var_est_NGR_bm = function(dt,
                          months = 1:12,
                          validation_years = 2001:2010)
{
  na_loc = which(dt[,is.na(SST_bar) | is.na(Ens_bar)])
  dt_new = dt[-na_loc,]
  
  #grouped by month:
  
  var_est_by_month = as.data.table(expand.grid(validation_years,months))
  setnames(var_est_by_month,c("year","month"))
  
  
  print("minimize CRPS for data grouped by month:")
  for(y in validation_years)
  {
    print(y)
    # grouped by month
    for(m in months)
    {
      temp = dt_new[year < y & month == m,][,Ens_var := Ens_sd^2]
      score_by_month = function(cd)
      {
        return(mean((temp[,var_bar]  - (cd[1]^2 + cd[2]^2 * temp[,Ens_var]))^2, na.rm = TRUE))
      } 
      opt_par = optim(par = c(0,1),fn = score_by_month)
      var_est_by_month[year == y & month == m, "c":= opt_par$par[1]]
      var_est_by_month[year == y & month == m, "d":= opt_par$par[2]]
    }
  }
  
  if("c" %in% colnames(dt))
  {
    dt[,c("c","d"):=NULL]
  }
  
  dt = merge(dt,var_est_by_month,by = c("year","month"),all.x = TRUE)
  dt[,SD_hat_lr_bm := sqrt(c^2 + d^2*Ens_sd^2)]
  
  return(dt)
}


#' Estimates the variance as has been suggested for NGR, grouped by location
#' 
#' @param dt the data table.
#' @param months The considered months.
#' @param validation_years Training years and validation years.
#' @param mc.cores Should we parallelize and with how many cores. Large datasets run into memory issues if mc.cores>1.
#' 
#' @return data table containing the RMSEs for the model above, where the coefficients are estimated in three different ways: grouped by month, location and by both.
#'
#'   
#' @author Claudio Heinrich
#' 
#' @export


var_est_NGR_bl = function(dt,
                          months = 1:12,
                          validation_years = 2001:2010,
                          mc.cores = 1)
{
  dt_new = dt[!(is.na(SST_bar) | is.na(Ens_bar)),]
  
  #grouped by location:
  
  print("minimize variance score for data grouped by location:")
  
  if(mc.cores == 1)
  {
    dt_gids = dt_new[month == min(month) & year == min(year),grid_id]
    
    var_est_bl = as.data.table(expand.grid(validation_years,dt_gids))
    setnames(var_est_bl,c('year','grid_id'))
    
    
    
    n_gid = length(dt_gids)
    brk = min(ceiling(n_gid/50),100)
    
    for(y in validation_years) 
    {
      print(y)
      temp_1 = dt_new[year < y, ]
      
      i = 0
      for(gid in dt_gids)
      {
        if(gid %% ceiling(n_gid/brk) == 0)
        {
          i=i+100/brk
          print(paste0(y,': ',floor(i),'%'))
        }
        temp = temp_1[ grid_id == gid,][,Ens_var := Ens_sd^2]
        score_by_gid = function(cd)
        {
          return(mean((temp[,var_bar]  - (cd[1]^2 + cd[2]^2 * temp[,Ens_var]))^2, na.rm = TRUE))
        } 
        opt_par = optim(par = c(0,1),fn = score_by_gid)
        var_est_bl[year == y & grid_id == gid, "c" := opt_par$par[1]]
        var_est_bl[year == y & grid_id == gid, "d" := opt_par$par[2]]
      }
      print(paste0(y,': 100%'))
    }
    
    
    if("c" %in% colnames(dt))
    {
      dt[,c("c","d"):=NULL]
    }
    
    dt = merge(dt,var_est_bl,by = c("year","grid_id"),all.x = TRUE)
    dt[,SD_hat_lr_bl := sqrt(c^2 + d^2*Ens_sd^2)]
    
  }
  
  # parallelize for reasonably small data tables
  
  if(mc.cores > 1)
  {
    dt_gids = dt_new[month == min(month) & year == min(year),grid_id]
    
    vebl_parallel = function(y)
    { return_DT = data.table(grid_id = dt_gids, year = y)
    
    temp_1 = dt_new[year < y, ]
    
    for(gid in dt_gids)
    {
      temp = temp_1[ grid_id == gid,][,Ens_var := Ens_sd^2]
      score_by_gid = function(cd)
      {
        return(mean((temp[,var_bar]  - (cd[1]^2 + cd[2]^2 * temp[,Ens_var]))^2, na.rm = TRUE))
      } 
      opt_par = optim(par = c(0,1),fn = score_by_gid)
      return_DT[grid_id == gid, "c" := opt_par$par[1]]
      return_DT[grid_id == gid, "d" := opt_par$par[2]]
    }
    return(return_DT)
    }
    
    var_est_by_loc = parallel::mclapply(validation_years,vebl_parallel,mc.cores =  mc.cores)
    var_est_by_loc = rbindlist(var_est_by_loc)
    
    if("c" %in% colnames(dt))
    {
      dt[,c("c","d"):=NULL]
    }
    
    dt = merge(dt,var_est_by_loc,by = c("year","grid_id"),all.x = TRUE)
    dt[,SD_hat_lr_bl := sqrt(c^2 + d^2*Ens_sd^2)]
  }
  
  
  return(dt) 
}

#' Estimates the variance as has been suggested for NGR, grouped by both month and location
#' 
#' @param dt the data table.
#' @param months The considered months.
#' @param validation_years Training years and validation years.
#' @param mc.cores Should we parallelize and with how many cores. Large datasets run into memory issues if mc.cores>1.
#' 
#' @return data table containing the RMSEs for the model above, where the coefficients are estimated in three different ways: grouped by month, location and by both.
#'
#'   
#' @author Claudio Heinrich
#' 
#' @export


var_est_NGR_bb = function(dt,
                          months = 1:12,
                          validation_years = 2001:2010,
                          mc.cores = 1)
{
  na_loc = which(dt[,is.na(SST_bar) | is.na(Ens_bar)])
  dt_new = dt[-na_loc,][,c:=NA][,d:=NA]
  
  if(mc.cores == 1)
  {
    dt_gids = dt_new[month == min(month) & year == min(year),grid_id]
    
    var_est_bb = as.data.table(expand.grid(validation_years,months,dt_gids))
    setnames(var_est_bb,c('year','month','grid_id'))
    
    # get number of gridpoints and how often to print progress
    n_gid = length(dt_gids)
    brk = min(ceiling(n_gid/50),100)
    
    for(y in validation_years) 
    {
      for(m in months)
      {
        temp_1 = dt_new[year < y & month == m,] 
        i = 0
        for(gid in dt_gids)
        {
          if(gid %% ceiling(n_gid/brk) == 0)
          {
            i=i+100/brk
            print(paste0(y,', ',m,': ',i,'%'))
          }
          temp = temp_1[grid_id == gid,][,Ens_var := Ens_sd^2]
          score_by_gid = function(cd)
          {
            return(mean((temp[,var_bar]  - (cd[1]^2 + cd[2]^2 * temp[,Ens_var]))^2, na.rm = TRUE))
          } 
          opt_par = optim(par = c(0,1),fn = score_by_gid)
          var_est_bb[year == y & month == m & grid_id == gid, "c" := opt_par$par[1]]
          var_est_bb[year == y & month == m & grid_id == gid, "d" := opt_par$par[2]]
        }
        print(paste0(y,', ',m,': 100%'))
      }
    }
    
    if("c" %in% colnames(dt))
    {
      dt[,c("c","d"):=NULL]
    }
    
    dt = merge(dt,var_est_bb,by = c("year",'month',"grid_id"),all.x = TRUE)
    dt[,SD_hat_lr_bb := sqrt(c^2 + d^2*Ens_sd^2)]
    
  }
  
  #parallelized version:
  if(mc.cores > 1)
  {
    
    dt_gids = dt_new[month == min(month) & year == min(year),grid_id]
    
    print("minimize variance score for data grouped by both:")
    
    res_dt = list()
    for (y in validation_years)
    {print(y)
      vebb_parallel = function(m)
      {
        temp = dt_new[year < y & month == m,]
        
        return_DT = data.table(grid_id = dt_gids)
        return_DT[,"year":=y][,"month":=m]
        for(gid in dt_gids)
        {
          temp_2 = temp[grid_id == gid,][,Ens_var := Ens_sd^2]
          score_by_both = function(cd)
          {
            return(mean((temp_2[,var_bar] - (cd[1]^2 + cd[2]^2 * temp_2[,Ens_var]))^2, na.rm = TRUE))
          } 
          opt_par = optim(par = c(0,1),fn = score_by_both)
          return_DT[ grid_id == gid, "c":= opt_par$par[1]]
          return_DT[ grid_id == gid, "d":= opt_par$par[2]]
        }  
        return(return_DT)
      }
      var_est_by_both = parallel::mclapply(months,vebb_parallel,mc.cores =  mc.cores)
      var_est_by_both = rbindlist(var_est_by_both)
      
      res_dt = rbindlist(list(res_dt,var_est_by_both[,year := y]))
    }
    
    if("c" %in% colnames(dt))
    {
      dt[,c("c","d"):=NULL]
    }
    
    dt  = merge(dt,res_dt,by = c("year","grid_id","month"),all.x = TRUE,)
    
    dt[,SD_hat_lr_bb := sqrt(c^2 + d^2*Ens_sd^2)]
    
  }
  
  return(dt)
}




