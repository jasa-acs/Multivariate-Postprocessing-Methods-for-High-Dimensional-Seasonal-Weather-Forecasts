########################################################################################
############### Creates a postprocessed forecast by Schaake Shuffle ####################
########################################################################################


#' Schaake shuffle forecasts 
#' 
#' @description Generates forecasts post-processed by Schaake shuffle
#' 
#' @param dt the data table containing bias and variance estimates.
#' @param obs_ens_years years to be used as historical archive of observations to determine dependence structure.
#' @param Y,M vectors contining the years and months for the forecast, default are all months in dt, and all years that are not contained in obs_ens_years.
#' @param method takes "R" or "Q", similar as ECC, see forecast_ECC.
#' 
#' @return data table containing \code{ens_size = length(obs_ens_years)} forecast columns, labelled "fc"1:ens_size, and the corresponding noise columns, defined as fc - SST_hat
#' 
#' @author Claudio Heinrich        
#' 
#' @importFrom matrixStats rowRanks
#' 
#' @export



forecast_Schaake = function(dt,obs_ens_years,
                            Y = NULL, M = NULL, 
                            method = "Q")
{
  if(length(intersect(Y,obs_ens_years)) > 0)
  {
    warning('forecast years and historical catalogue are overlapping!')
  }
  
  # reduce data 
  
  if(!is.null(M))
  {
    dt = dt[month %in% M,]
  }  
  
  # historical catalogue
  
  dt_hist = dt[year %in% obs_ens_years]
  
  if(is.null(Y))
  {
    Y = dt[,unique(year)]
    Y = Y[!(Y %in% obs_ens_years)]
  }
  
  
  dt = dt[year %in% Y,]
  
  # remove land:
  
  dt_temp = dt[!(is.na(Bias_Est) | is.na(Ens_bar) | is.na(SD_hat))]
  dt_hist = dt_hist[!(is.na(Bias_Est) | is.na(Ens_bar) | is.na(SD_hat))]
  
  na_rows = dt[,(is.na(Bias_Est) | is.na(Ens_bar) | is.na(SD_hat))]
  
  # do univariate postprocessing:  
  
  print('univariate postprocessing...')
  
  ens_size = length(obs_ens_years)
  
  length_norm = dt_temp[,.N]
  
  if(method == "R")
  {
    for(i in 1:ens_size)
    {
      dt_temp[ ,paste0("fc",i) := rnorm(n = length_norm, 
                                        mean = dt_temp[,SST_hat], 
                                        sd = dt_temp[,SD_hat])]
    }
  }
  
  if(method == "Q")
  {
    for(i in 1:ens_size)
    {
      dt_temp[ ,paste0("fc",i) := qnorm(i/(ens_size + 1), 
                                        mean = dt_temp[,SST_hat], 
                                        sd = dt_temp[,SD_hat])]
    }
  }
  
  
  # get rank order of the historical catalogue and reorder the post-processed ensemble to match the rank order statistic of the ensemble.
  # the following is slow, would be nice to find a faster option:
  
  
  #get historical catalogue:
  
  ind = 0 
  
  hist_obs = NULL
  
  for(yy in obs_ens_years)
  {
    ind = ind + 1
    hist_obs[[ind]] = dt_hist[year == yy, SST_bar]
  }
  
  # plug into forecast dt
  
  for(yy in dt_temp[,unique(year)])
  {
    dt_temp[year == yy, paste0('Ens',1:ens_size) := hist_obs]
  }
  
  
  print('computing ranks...')
  
  rks_ens = t(apply(dt_temp[,.SD,.SDcols = paste0("Ens",1:ens_size)], MARGIN = 1, rank, ties = 'random'))
  rks_fc = t(apply(dt_temp[,.SD,.SDcols = paste0("fc",1:ens_size)], MARGIN = 1, rank, ties = 'random'))
  
  #Here is a faster option that does almost the same, but unfortunately it cannot resolve ties at random:
  #rks_ens = dt_temp[,matrixStats::rowRanks(as.matrix(.SD)),.SDcols = paste0("Ens",1:ens_size)]
  #rks_fc = dt_temp[,matrixStats::rowRanks(as.matrix(.SD)),.SDcols = paste0("fc",1:9)]
  
  
  fcs = as.matrix(dt_temp[,.SD,.SDcols = paste0("fc",1:ens_size)])
  
  num_row = nrow(fcs)
  
  schaake = matrix(NA,nrow = num_row,ncol = ens_size)
  for(i in 1:num_row){
    if( i %% 100000 == 0)
    {
      print(paste0(i," / ",num_row))
    }
    if(!is.na(fcs[i,1]))
    { 
      schaake[i,] = fcs[i, match(rks_ens[i,],rks_fc[i,])]
    }
  }
  
  schaake = data.table(schaake)
  setnames(schaake,paste0("fc",1:ens_size))
  
  dt_temp[,paste0("fc",1:ens_size) := NULL]
  
  dt_temp = data.table(dt_temp,schaake)
  
  # put back missing locations
  
  if(sum(na_rows) > 0)
  {
    key_dt = key(dt)
    dt = rbindlist(list(dt_temp,dt[na_rows,]),fill = TRUE)
    setkeyv(dt,key_dt)
  } else {
    dt = dt_temp
  }
  
  
  return(dt)
}
