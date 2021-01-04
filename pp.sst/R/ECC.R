
############################################################################
############### Creates a postprocessed forecast by ECC ####################
############################################################################


#' ECC forecasts 
#' 
#' @description Generates forecasts post-processed by ensemble copula coupling (ECC)
#' 
#' @param dt the data table containing bias and variance estimates.
#' @param Y,M vectors contining the years and months for the forecast, default are all years and months in dt.
#' @param ens_size Integer. Size of the NWP ensemble.
#' @param method takes "R" or "Q". Method of ECC, see Schefzik et al. 2013.
#' 
#' @return data table containing \code{ens_size} forecast columns, labelled "fc"1:ens_size, and the corresponding noise columns, defined as fc - SST_hat
#' 
#' @examples \dontrun{}
#' 
#' @author Claudio Heinrich        
#' 
#' @importFrom matrixStats rowRanks
#' 
#' @export



forecast_ECC = function(dt, Y = NULL, M = NULL,
                        ens_size = 9,
                        method = "Q")
{
  # prepare data table:
  
  if(!is.null(Y))
  {
    dt = dt[year %in% Y,]
  }
  
  
  if(!is.null(M))
  {
    dt = dt[month %in% M,]
  }  
  
  
  dt_temp = dt[!(is.na(Bias_Est) | is.na(Ens_bar) | is.na(SD_hat))]
  
  # do univariate postprocessing:  
  
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
  
  
  # get rank order of the ensemble and reorder the post-processed ensemble to mach the rank order statistic of the ensemble.
  # the following is slow, would be nice to find a faster option:
  rks_ens = t(apply(dt_temp[,.SD,.SDcols = paste0("Ens",1:ens_size)], MARGIN = 1, rank, ties = 'random'))
  rks_fc = t(apply(dt_temp[,.SD,.SDcols = paste0("fc",1:ens_size)], MARGIN = 1, rank, ties = 'random'))
  
  #Here is a faster option that does almost the same, but unfortunately it cannot resolve ties at random:
  #rks_ens = dt_temp[,matrixStats::rowRanks(as.matrix(.SD)),.SDcols = paste0("Ens",1:ens_size)]
  #rks_fc = dt_temp[,matrixStats::rowRanks(as.matrix(.SD)),.SDcols = paste0("fc",1:9)]
  
  
  fcs = as.matrix(dt_temp[,.SD,.SDcols = paste0("fc",1:ens_size)])
  
  num_row = nrow(fcs)
  
  ecc_fcs = matrix(NA,nrow = num_row,ncol = ens_size)
  for(i in 1:num_row){
    if( i %% 100000 == 0)
    {
      print(paste0(i," / ",num_row))
    }
    if(!is.na(fcs[i,1]))
    { 
      ecc_fcs[i,] = fcs[i, match(rks_ens[i,],rks_fc[i,])]
    }
  }
  
  ecc_fcs = data.table(ecc_fcs)
  setnames(ecc_fcs,paste0("fc",1:ens_size))
  
  dt_temp[,paste0("fc",1:ens_size) := NULL]
  
  dt_temp = data.table(dt_temp,ecc_fcs)
  
  dt_temp[,lapply(X = .SD,FUN = trc),.SDcols = paste0("fc",1:ens_size)]
  
  # put back missing locations

  if(!(dt[is.na(Bias_Est) | is.na(Ens_bar) | is.na(SD_hat),.N]==0))
  {
    key_dt = key(dt)
    dt = rbindlist(list(dt_temp,dt[is.na(Bias_Est) | is.na(Ens_bar) | is.na(SD_hat),]),fill = TRUE)
    setkeyv(dt,key_dt)
  } else {
    dt = dt_temp
  }
  
  
  return(dt)
}





