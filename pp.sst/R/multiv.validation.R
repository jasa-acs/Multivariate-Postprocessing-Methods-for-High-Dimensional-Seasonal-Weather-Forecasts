#########################################
### tools for multivariate validation ###
#########################################


#' new, faster plotting of multivariate rank histograms
#' 
#' @description computes the multivariate rank histograms for the provided forecast data table. Plots the average rank and the band depth rankq.
#'  
#' @param dt_fc the forecast data table.
#' @param fc_ens_size Size of the forecast ensemble: dt_fc ought to have columns named fc1,...,fcn, where n = fc_ens_size.
#' @param mc_cores number of cores used for parallelization.
#' @param mn Title of the plot.
#' @param breaks Number of breaks in the rank histograms. Should at most be 1 + fc_ens_size (which is the size of the ensemble for ECC).
#' @param save_pdf,plot_dir,file_name Whether and where the rank histogram should be saved as pdf. \code{file_name} should not contain '.pdf'.
#' 
#' @return if save_pdf = FALSE the rank histogram is shown. If save_pdf = TRUE, the plots are saved and a data table with the ranks of the observations is returned.
#' 
#' @examples \dontrun{}
#' 
#' @author Claudio Heinrich        
#' 
#' @export

mv_rank_hist_new = function(dt_fc,
                            fc_ens_size,
                            mc_cores = 12,
                            mn = "",
                            breaks = 11,
                            save_pdf = FALSE, plot_dir = "", file_name = "")
{
  dt_fc = dt_fc[!is.na(SST_hat)]
  
  ym = unique(dt_fc[,YM])
  
  ym_ind = 1:length(ym)
  
  ranks_dt = data.table(YM = ym)
  
  ranks_matrix = matrix(ym,nrow = length(ym),ncol = 1+2*(fc_ens_size +1)) 
  #ncol: 1 col for YM, 2 methods of ranking, for each we get ranks for observation and each Monte Carlo sample
  
  drm = dim(ranks_matrix)
  
  YM_ind = 0
  
  ranks_ym = function(yearmonth)
  {print(floor(yearmonth/12))
    fc_obs_dt = dt_fc[YM == yearmonth,.SD,.SDcols = c("SST_bar",paste0("fc",1:fc_ens_size))]
    
    temp_avg = avg.rank(fc_obs_dt)[1]
    temp_bd = bd.rank(fc_obs_dt)[1]
    
    return_dt = data.table(YM = yearmonth,av.rk.obs = temp_avg, bd.rk.obs = temp_bd)
    return(return_dt)
  }
  
  ranks = rbindlist(parallel::mclapply(X = ym,
                                       FUN = ranks_ym,
                                       mc.cores = mc_cores,
                                       mc.silent = FALSE))
  
  if(save_pdf)
  {
    pdf(file=paste0(plot_dir,file_name,".pdf"),width=24,height=7,points=12)
  }
  
  par(mfrow=c(1,2),mex=0.5,oma = c(0,0,2.5,0),mar=c(2.5,2.5,2.5,2.5)+0.1,mgp=c(0.5,0,0))
  
  rhist_dt(ranks[,.(YM,av.rk.obs)], max_rk = fc_ens_size +1, breaks = breaks, hist_xlab = "average")
  rhist_dt(ranks[,.(YM,bd.rk.obs)], max_rk = fc_ens_size +1, breaks = breaks, hist_xlab = "band depth")
  
  title(mn,outer = TRUE)
  
  if(save_pdf)
  {
    dev.off()  
    return(ranks)}
}


##### multivariate rank histograms ######

#' making and plotting of multivariate rank histograms
#' 
#' @description computes the multivariate rank histograms for the provided forecast data table. Plots the average rank, the band depth rank, and the minimum spanning tree rank.
#'  
#' @param dt_fc the forecast data table.
#' @param fc_ens_size Size of the forecast ensemble: dt_fc ought to have columns named fc1,...,fcn, where n = fc_ens_size.
#' @param mn Title of the plot
#' @param breaks Number of breaks in the rank histograms. Should at most be 1 + fc_ens_size (which is the size of the ensemble for ECC).
#' @param save_pdf,plot_dir,file_name Whether and where the rank histogram should be saved as pdf. \code{file_name} should not contain '.pdf'.
#' 
#' @return none, if save_pdf = FALSE, the rank histogram is shown
#' 
#' @examples \dontrun{}
#' 
#' @author Claudio Heinrich        
#' 
#' @export

mv_rank_hist = function(dt_fc,
                        fc_ens_size,
                        mn = "",
                        breaks = 10,
                        save_pdf = FALSE, plot_dir = "", file_name = "")
{
  dt_fc = dt_fc[!is.na(SST_hat)]
  
  ym = unique(dt_fc[,YM])
  
  ranks_matrix = matrix(ym,nrow = length(ym),ncol = 1+3*(fc_ens_size +1)) 
  #ncol: 1 col for YM, 3 methods of ranking, for each we get ranks for observation and each Monte Carlo sample
  
  drm = dim(ranks_matrix)
  
  YM_ind = 0
  
  for(yearmonth in ym)
  {
    print(paste0("YM = ",yearmonth-min(ym)+1,"/",ym[length(ym)]-min(ym)+1))  
    
    YM_ind = YM_ind + 1
    fc_obs_mat = dt_fc[YM == yearmonth,.SD,.SDcols = c("SST_bar",paste0("fc",1:fc_ens_size))]
    
    # get ranks
    ranks_matrix[YM_ind,2:drm[2]] = c(mst.rank(as.matrix(fc_obs_mat)),
                                      avg.rank(as.matrix(fc_obs_mat)),
                                      bd.rank(as.matrix(fc_obs_mat)))
  }
  
  names_vec = c("YM","mst.rk.obs",paste0("mst.r.",1:fc_ens_size),
                "av.rk.obs",paste0("av.r.",1:fc_ens_size),
                "bd.rk.obs",paste0("bd.rk",1:fc_ens_size))
  
  ranks = data.table(ranks_matrix)
  
  setnames(ranks, names_vec)
  
  if(save_pdf)
  {
    pdf(file=paste0(plot_dir,file_name,".pdf"),width=8,height=2,points=12)
  }
  
  par(mfrow=c(1,3),mex=0.5,oma = c(0,0,2.5,0),mar=c(2.5,2.5,2.5,2.5)+0.1,mgp=c(0.5,0,0))
  
  rhist_dt(ranks[,.(YM,mst.rk.obs)], max_rk = fc_ens_size +1, breaks = breaks,  hist_xlab = "minimum spanning tree")
  rhist_dt(ranks[,.(YM,av.rk.obs)],max_rk = fc_ens_size +1, breaks = breaks, hist_xlab = "average")
  rhist_dt(ranks[,.(YM,bd.rk.obs)], max_rk = fc_ens_size +1, breaks = breaks, hist_xlab = "band depth")
  
  title(mn,outer = TRUE)
  
  if(save_pdf)
  {
    dev.off()  }
  
}

########################


### preranking functions from 'Assessing the Calibration of High-Dimensional Ensemble Forecasts Using Rank Histograms' ###

#' Minimum spanning tree ranks
#' @export

mst.rank <- function (x) {
  l.mst <- NULL
  for(f in 1:(dim(x)[2])) {
    euc.dist <- rdist::rdist(t(x[,-f]))
    l.mst <- c(l.mst,sum(vegan::spantree(euc.dist)$dist))
  }
  x.rank <- rank(l.mst,ties="random")
  return(x.rank)
}

########################


#' Multivariate ranks 
#' @export

mv.rank <- function(x)
{
  d <- dim(x)
  x.prerank <- numeric(d[2])
  for(i in 1:d[2]) {
    x.prerank[i] <- sum(apply(x<=x[,i],2,all))
  }
  x.rank <- rank(x.prerank,ties="random")
  return(x.rank)
}

########################


#' Average ranks
#' @export

avg.rank <- function(x)  {
  x.ranks <- apply(x,1,rank,ties = 'random')
  x.preranks <- apply(x.ranks,1,mean)
  x.rank <- rank(x.preranks,ties="random")
  return(x.rank)
}

########################


#' Band depth ranks
#' @export

bd.rank <- function(x)
{
  d <- dim(x)
  x.prerank <- array(NA,dim=d)
  for(i in 1:d[1]) {
    tmp.ranks <- rank(x[i,],ties = 'random')
    x.prerank[i,] <- (d[2] - tmp.ranks) * (tmp.ranks - 1)
  }
  x.rank <- apply(x.prerank,2,mean) + d[2] - 1
  x.rank <- rank(x.rank,ties="random")
  return(x.rank)
} 

#######################


#' rank histograms for data tables
#' 
#' @param B A data table, containing the key variable (most commonly YM) as first column and the ranks of the observations as second
#' @param breaks Number of breaks desired for the rank histogram.
#' @param hist_xlab,hist_ylab Labels for axis.
#' @param hist_ylim Range of y-axis.
#' 
#' @author Claudio Heinrich
#'  
#' @export

rhist_dt <- function(B, max_rk, breaks = 11, hist_xlab=NULL, hist_ylab=NULL, hist_ylim=NULL,mn = '')
{
  hist(as.vector(B[[2]]),breaks=seq(0, max_rk, length.out = breaks), 
       main=mn,
       xlab=hist_xlab, ylab=hist_ylab,
       axes=FALSE, col="gray80", border="gray60", ylim=hist_ylim)
  abline(a=length(B[[1]])/(breaks-1), b=0, lty=2, col="gray30")
}


#' histogram transformation 
#' 
#' when working with rank histograms problems arise when the number of bins does not divide the number of possible ranks
#' This function solves this problem by randomizing the data such that it is uniform on [0,1] if it was uniform before
#' 
#' @param vec A vector containing ranks. 
#'
#' @return A transformed vector containing data on [0,1]
#'   
#' @author Claudio Heinrich
#'  
#' @export

hist_transform <- function(vec)
{
  randomize = runif(n = length(vec))
  return( (vec - 1 + randomize))
}



############################################
####### 'Empirical' variogram scores #######
############################################


#' computing variogram score for a given month and year
#' 
#' @param m,y Month and year.
#' @param dt_fc The forecast data table.
#' @param n number of forecasts available for moment estimation.
#' @param p Power used in the variogram score.
#' 
#' @return a data table with one row and three columns: year, month and variogram score.
#' 
#' @author Claudio Heinrich
#'  
#' @export

var_sc_emp = function(m,y,dt_fc,n,p)
{
  # setting up:
  
  dt_fc = dt_fc[year %in% y & month %in% m,]
  
  land_ids <- which(dt_fc[, is.na(Ens_bar) | is.na(SST_bar)])
  if(!identical(land_ids,integer(0)))
  {
    dt_fc = dt_fc[-land_ids,]
  }
  
  # build data table that contains pairs of coordinates for all coordinates contained in dt:
  
  dt_coor_1 = dt_fc[,.(Lon,Lat,grid_id)]
  setnames(dt_coor_1,c("Lon1","Lat1","grid_id1"))
  # add dummy key, do outer full merge with a duplicate, and remove key again:
  dt_coor_1[,"key" := 1]
  dt_coor_2 = copy(dt_coor_1) 
  setnames(dt_coor_2,c("Lon2","Lat2","grid_id2","key"))
  var_sc_prereq = merge(dt_coor_1,dt_coor_2, by = "key",allow.cartesian = TRUE)
  var_sc_prereq[, "key" := NULL]
  
  # get pairs of indices of locations
  
  n_loc = dt_coor_1[,.N]
  id_1 = as.vector(mapply(rep,1:n_loc,times = n_loc))
  id_2 = rep(1:n_loc, times = n_loc)
  
  # get weights for variogram score
  
  var_sc_prereq[,"weights" := 1/.N]
  
  
  ##--- For each pair of coordinates (i,j), estimate the pth absolute moment of X_i-X_j 
  
  
  mom_est = 0
  
  for(i in 1:n)
  {
    mom_est = mom_est + (1/n) * abs(dt_fc[,.SD,.SDcols = paste0("fc",i)][id_1] - dt_fc[,.SD,.SDcols = paste0("fc",i)][id_2])^p
  }
  
  var_sc_prereq[,mom_est := mom_est]
  
  #get observation
  var_sc_prereq[,mom_obs := abs(dt_fc[,SST_bar][id_1] - dt_fc[,SST_bar][id_2])^p]
  
  vs = var_sc_prereq[,sum(weights * (mom_obs - mom_est)^2)] 
  
  return(data.table(year = y, month = m, vs = vs)) 
}

########################


#' computing empirical variogram scores
#' 
#' @description Computes the variogram scores for each month and year for a given forecast data table. Parallelized along months.
#' 
#' @param dt_fc The forecast data table.
#' @param years,ms The years and months for which to compute the variogram score.
#' @param n number of forecasts available for moment estimation.
#' @param p Power used in the variogram score.
#' @param save_dir,file_name If we should save and where. \code{file_name} ought not contain '.RData'.
#' @param mc.cores How many cores do you want to use?
#' @return a data table with three columns: year, month and vs.
#' 
#' @author Claudio Heinrich
#'  
#' @export

var_sc_par = function(dt_fc, years, ms, n , p = 0.5,
                      save_dir = NULL, file_name = "vs",
                      mc_cores = NULL)
{
    if(is.null(mc_cores)) mc_cores = length(ms)
    
    vs_all = list()
    for(i in 1:length(years))
    {
        year_i = years[i]
        print(year_i)
        ## to avoid errors in mclapply:
        if(mc_cores == 1){
            l = list()
            for(j in 1:length(ms)){
                l[[j]] = var_sc_emp(m = ms[j],
                                    y = year_i,
                                    dt_fc = dt_fc,
                                    n = n,
                                    p = p)
            }
        }else{
            l = parallel::mclapply(X = ms,
                                   FUN = var_sc_emp, 
                                   y = year_i,
                                   dt_fc = dt_fc,
                                   n = n,
                                 p = p,
                                   mc.cores = mc_cores)
        }
        vs_all[[i]] = rbindlist(l)
    }
    vs = rbindlist(vs_all)
    if(!is.null(save_dir))
    {
        save(vs,file = paste0(save_dir,file_name,".RData"))
    }
  
    return(vs)
}
