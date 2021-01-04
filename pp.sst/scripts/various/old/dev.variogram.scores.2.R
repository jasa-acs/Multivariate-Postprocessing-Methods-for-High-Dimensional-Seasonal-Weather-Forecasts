
######################################################
###### Functions for computing variogram scores ######
######################################################


#' Computing variogram scores
#' 
#' Computes the variogram for a normal forecasting distribution
#'
#' @param dt a data table containing the columns .(sq_mean_diff, dvar, sq_obs_diff), containing the squared differences of forecast means, the variances for differences of the forecast distribution,
#' as well as the square of pairwise differences of observations.
#' @param p power for the variogram score.
#' 
#' @return a double var_sc
#' 
#' @author Claudio Heinrich
#' 
#' @export


variogram_score_nrm_p2 = function(dt,weights = TRUE)
  {
  
  # get the p-th moment of the forecasting distribution, which is non-central Gaussian
  if(!weights)
  {
    var_sc = sum(dt[,(sq_mean_diff + dvar - sq_obs_diff)^2])/dt[,.N]  
  }
  if(weights)
  {
    var_sc = sum(dt[,weights * (sq_mean_diff + dvar - sq_obs_diff)^2])
  }
  
  return(var_sc)
}


get_moment_nrm = function(dt,p = 0.5)
{
  fc_mean = dt[,fc_mean]
  fc_var = dt[,fc_var]
  
  # get the p-th abs. moment of the forecasting distribution, which is non-central Gaussian
  f = 2^(p/2) * gamma((p+1)/2)/sqrt(pi)
  p_mom := f * fc_var^(p/2) * hypergeo::genhypergeo(-p/2,1/2,-1/2*(fc_mean/sqrt(fc_var))^2)
  return(p_mom)
}

#' Computing variogram scores
#'
#' @param dt a data table containing the columns .(year, grid_id1,grid_id2,fc_var, obs_diff) 
#' @param p power for the variogram score.
#' @param eval_years integer vector containing the years you want to use for validation, need to be contained in dt[,year]
#' 
#' @return a double var_sc
#' 
#' @author Claudio Heinrich
#' 
#' @export


variogram_score_old = function(dt,
                           p = 0.5,
                           eval_years=2001:2010){
  
  
  # get the p-th moment of a standard normal distribution
  p_mom_sn = 2^(p/2)*gamma((p+1)/2)/sqrt(pi) 
  
  var_sc = sum(((dt[,fc_var])^(p/2) * p_mom_sn - 
                  (dt[,obs_diff])^p)^2)
  
  return(var_sc)
}


#####################################################
####################### PCA: ########################
#####################################################


var_sc_PCA_standardized = function(m, y, dt,
                                   PCA = NULL, PCA_DT = NULL,
                                   dvec = c(1:10,12,14,16,18,20,25,30,50,70),
                                   marginal_correction = TRUE,
                                   cov_dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                                   ens_size = 9,
                                   save_dir = "~/PostClimDataNoBackup/SFE/Derived/PCA/",
                                   file_name = "var_sc_by_PC_stan",
                                   weighted = FALSE,
                                   weight_fct = NULL)
{
  if(weighted)
  {
    # the weights are a function of the SST difference averaged over training_years
    weight_dt = dt[year %in% training_years & month == m,][,mean(SST_bar), by = grid_id]
    setnames(weight_dt,c("grid_id","mean_SST"))
  }
  
  # setting up:
  
  dt = dt[year %in% y & month %in% m,]
  
  land_ids <- which(dt[, is.na(Ens_bar) | is.na(SST_bar)])
  if(!identical(land_ids,integer(0)))
  {
    dt = dt[-land_ids,]
    if(weighted)
    {
      weight_dt = weight_dt[-land_ids,]
    }
  }
  
  
  trc_level = 0.01 
  dt[clim_sigma < trc_level, clim_sigma := trc_level]
  
  clim = dt[,clim]
  clim_sd = dt[,clim_sigma]
  
  
  #get covariance matrix
  
  if(is.null(PCA))
  {
    load(file = paste0(cov_dir,"CovRes_mon",m,".RData"))
    
    PCA = irlba::irlba(res_cov, nv = max(dvec))
  }
  
  
  if(is.null(PCA_DT))
  {
    PCA_DT = dt[,.(Lon,Lat,grid_id)]
    
    for(d in  min(dvec):max(dvec)){
      PCA_DT [,paste0("PC",d) := PCA$d[d]*PCA$u[,d]]
    } 
    
    # also get marginal variances
    variances = list()
    d = min(dvec) 
    vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
    variances[[d]] = vec^2
    
    if(length(dvec)>1)
    {
      for(d in (min(dvec)+1) : max(dvec)){
        vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
        variances[[d]] = variances[[d-1]] + vec^2
      }  
    }
    
    names(variances) = paste0("var",min(dvec):max(dvec))
    PCA_DT = data.table(PCA_DT,rbindlist(list(variances)))  
  }
  
  
  # marginal SD correction factor:
  SD = dt[,SD_hat]
  
  crit_ind  = which(SD < 0.00001 | PCA_DT[,var1] < 1e-15)
  
  PCA_DT[,paste0("marSDcf",dvec) := SD/(clim_sd * sqrt(.SD)),.SDcols = paste0("var",dvec)]
  
  # at critical indices we just use the SD of the PCA:
  PCA_DT[crit_ind, paste0("marSDcf",dvec) := 1]
  
  
  
  # build data table that contains pairs of coordinates for all coordinates contained in dt:
  
  dt_coor_1 = PCA_DT[,.(Lon,Lat,grid_id)]
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
  
  #get weights for variogram score
  
  if(weighted)
  {
    if(is.null(weight_fct))
    {
      weight_fct = function(x)
      { y = rep(1, length(x))
      y[abs(x)>1] = (1/x[abs(x)>1]^2)
      return(y)
      }
    }
    
    var_sc_prereq[,"weights" := weight_fct(weight_dt[,mean_SST][id_1] - weight_dt[,mean_SST][id_2])]
    #normalizing:
    var_sc_prereq[,"weights" := weights/sum(weights)]
  } else {
    var_sc_prereq[,"weights" := 1/.N]
  }
  
  
  
  ##--- For each pair of coordinates (i,j) compute the variance of X_i-X_j 
  ##--- for the predictive distribution with d principal components
  
  print("getting variances of differences:")
  
  diff_var = list()
  
  if(marginal_correction)
  {
    list_page = 1
    for(d in dvec)
    {
      print(paste0("d = ",d))
      mar_cor = PCA_DT[,eval(parse(text = paste0("marSDcf",d)))]
      temp = mar_cor * PCA_DT[,.SD,.SDcols = paste0("PC",1:d)] 
      diff_var[[list_page]] =  rowSums((temp[id_1]-temp[id_2])^2)
      list_page = list_page + 1
    }  
    names(diff_var) = paste0("dvar",dvec)
    
  } else if(! marginal_correction)
  {
    d = min(dvec) 
    vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]/clim_sd
    diff_var[[1]] = (vec[id_1]-vec[id_2])^2
    
    list_page = 2
    for(d in (min(dvec)+1):max(dvec)){
      vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]/clim_sd
      diff_var[[list_page]] = diff_var[[list_page-1]] + (vec[id_1]-vec[id_2])^2
      list_page = list_page +1
    }
    names(diff_var) = paste0("dvar",min(dvec):max(dvec))
  }
  
  diff_var = rbindlist(list(diff_var))
  
  var_sc_prereq = data.table(var_sc_prereq,diff_var)
  
  #complement var_sc_prereq by the squared differences of the mean vectors and squared differences of observation:
  
  var_sc_prereq[,sq_m_diff := (dt[,(SST_hat-clim)/clim_sd][id_1]-dt[,(SST_hat-clim)/clim_sd][id_2])^2]
  
  var_sc_prereq[,sq_obs_diff := (dt[,(SST_bar-clim)/clim_sd][id_1]-dt[,(SST_bar-clim)/clim_sd][id_2])^2]
  
  
  # get variogram scores
  print("done - compute variogram scores:")
  scores = list()
  
  #var_sc_by_PC = function(d){
  list_page = 1
  
  for(d in dvec){
    print(paste0("d = ",d))
    return_data =  data.table(year = y, d = d)
    
    dt_temp = var_sc_prereq[,.SD,.SDcols = c("sq_m_diff",paste0("dvar",d),"sq_obs_diff","weights")]
    setnames(dt_temp,c("sq_mean_diff", "dvar", "sq_obs_diff","weights"))
    return_data[,sc := variogram_score_nrm_p2(dt_temp)]
    
    
    scores[[list_page]] = return_data[,month := m]  
    list_page = list_page +1
  }
  
  scores = rbindlist(scores)
  
  if(!marginal_correction)
  {
    file_name = paste0(file_name,"_nmc")
  }
  
  save(scores, file = paste0(save_dir,file_name,"_m",m,"_y",y,".RData"))
}



#' Computing variogram scores for PCA post-processing 
#'
#' @description  This function computes the variogram scores for
#' the PCA post-processing method for a given month and year for a range of considered numbers of PCs.
#'
#' @param m,y The month and year.
#' @param dt The data table.
#' @param p The power for the variogram score.
#' @param PCA,PCA_DT Optional. These objects depend only on the month, not the year. Therefore, to speed things up, they can be computed prior to running the function for a range of years.
#' @param dvec Integer vector. Contains the numbers of principal components we want to test for.
#' @param marginal_correction Logical. If TRUE, PCA with corrected marginal variance is considered (see paper).
#' @param cov_dir String. The covariance directory for the PCA.
#' @param ens_size The number of ensemble members.
#' @param save_dir The directory to save in.
#' @param file_name The file is saved as \code{file_name(_nmc)_m#_y##.RData}, where # is the month and ## the year and _nmc is added if no marginal correction is chosen, i.e. if marginal_correction == F.
#' 
#' @examples \dontrun{
#' save_dir = "~/PostClimDataNoBackup/SFE/Derived/NAO/"
#' DT = load_combined_wide(data_dir = save_dir, output_name = "dt_combine_wide_bc_var.RData")
#' var_sc_PCA(m = 1,y = 2000, dt = DT)}
#'
#' @author Claudio Heinrich
#' 
#' @export


var_sc_PCA_new = function(m, y, dt, weight_mat,
                          PCA = NULL, PCA_DT = NULL, 
                          dvec = c(1:10,12,14,16,18,20,25,30,50,70),
                          marginal_correction = TRUE,
                          cov_dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                          ens_size = 9,
                          save_dir = "~/PostClimDataNoBackup/SFE/Derived/PCA/",
                          file_name = "var_sc_by_PC",
                          weighted = FALSE,
                          weight_fct = NULL)
{
  # setting up:
  
  dt = dt[year %in% y & month %in% m,]
  
  land_ids <- which(dt[, is.na(Ens_bar) | is.na(SST_bar)])
  if(!identical(land_ids,integer(0)))
  {
    dt = dt[-land_ids,]
  }
  
  # get covariance matrix
  
  if(is.null(PCA))
  {
    load(file = paste0(cov_dir,"CovRes_mon",m,".RData"))
    
    PCA = irlba::irlba(res_cov, nv = max(dvec))
    
  }
  
  
  if(is.null(PCA_DT))
  {
    PCA_DT = dt[,.(Lon,Lat,grid_id)]
    
    for(d in  min(dvec):max(dvec)){
      PCA_DT [,paste0("PC",d) := PCA$d[d]*PCA$u[,d]]
    } 
  
    # also get marginal variances
    variances = list()
    d = min(dvec) 
    vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
    variances[[d]] = vec^2
    
    if(length(dvec)>1)
    {
      for(d in (min(dvec)+1) : max(dvec)){
        vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
        variances[[d]] = variances[[d-1]] + vec^2
      }  
    }
    
    names(variances) = paste0("var",min(dvec):max(dvec))
    PCA_DT = data.table(PCA_DT,rbindlist(list(variances)))  
  }
  
  
  # marginal SD correction factor:
  SD = dt[,SD_hat]
  
  crit_ind  = which(SD < 0.00001 | PCA_DT[,var1] < 1e-20)
  
  PCA_DT[,paste0("marSDcf",dvec) := SD/sqrt(.SD),.SDcols = paste0("var",dvec)]
  
  # at critical indices we just use the SD of the PCA:
  PCA_DT[crit_ind, paste0("marSDcf",dvec) := 1]

  
  
  # build data table that contains pairs of coordinates for all coordinates contained in dt:
  
  dt_coor_1 = PCA_DT[,.(Lon,Lat,grid_id)]
  setnames(dt_coor_1,c("Lon1","Lat1","grid_id1"))
  # add dummy key, do outer full merge with a duplicate, and remove key again:
  dt_coor_1[,"key" := 1]
  dt_coor_2 = copy(dt_coor_1) 
  setnames(dt_coor_2,c("Lon2","Lat2","grid_id2","key"))
  var_sc_prereq = merge(dt_coor_1,dt_coor_2, by = "key",allow.cartesian = TRUE)
  var_sc_prereq[, "key" := NULL]
  
  # get pairs of indices of locations
  
  n_loc_pair = dt_coor_1[,.N]
  id_1 = as.vector(mapply(rep,1:n_loc,times = n_loc_pair))
  id_2 = rep(1:n_loc_pair, times = n_loc_pair)
  
  
  #finding variance and covariance locations in Sigma:
  n_loc = PCA_DT[,.N]
  var_id_1 = id_1 + (id_1-1)*n_loc
  var_id_2 = id_2 + (id_2-1)*n_loc
  cov_id = id_1 + (id_2 - 1)*n_loc
  
  
  #get weights for variogram score
  
  if(weighted)
  {
    if(is.null(weight_fct))
    {
      weight_fct = function(x)
      { y = rep(1, length(x))
      y[abs(x)>1] = (1/x[abs(x)>1]^2)
      return(y)
      }
    }
    
    var_sc_prereq[,"weights" := weight_fct(weight_dt[,mean_SST][id_1] - weight_dt[,mean_SST][id_2])]
    #normalizing:
    var_sc_prereq[,"weights" := weights/sum(weights)]
  } else {
    var_sc_prereq[,"weights" := 1/.N]
  }
  
  
  
  ##--- For each pair of coordinates (i,j) compute the variance of X_i-X_j 
  ##--- for the predictive distribution with d principal components
  
    print("getting variances of differences:")
  
  diff_var = list()
  
  if(marginal_correction)
  {
    list_page = 1
    for(d in dvec)
    {
      print(paste0("d = ",d))
      mar_cor = PCA_DT[,eval(parse(text = paste0("marSDcf",d)))]
      temp = as.matrix(mar_cor * PCA_DT[,.SD,.SDcols = paste0("PC",1:d)]) 
      Sigma_temp = weight_mat * (temp %*% t(temp))
      diff_var[[list_page]] =  Sigma_temp[var_id_1] + Sigma_temp[var_id_2] - 2*Sigma_temp[cov_id]
      list_page = list_page + 1
    }  
    names(diff_var) = paste0("dvar",dvec)
    
  } else if(!marginal_correction)  {
    list_page = 1
    for(d in dvec)
    {
      print(paste0("d = ",d))
      temp = as.matrix( PCA_DT[,.SD,.SDcols = paste0("PC",1:d)]) 
      Sigma_temp = weight_mat * (temp %*% t(temp))
      diff_var[[list_page]] =  Sigma_temp[var_id_1] + Sigma_temp[var_id_2] - 2*Sigma_temp[cov_id]
      list_page = list_page + 1
    }  
    names(diff_var) = paste0("dvar",dvec)
    
  }
  
  diff_var = rbindlist(list(diff_var))
  
  var_sc_prereq = data.table(var_sc_prereq,diff_var)
  
  #complement var_sc_prereq by the squared differences of the mean vectors and squared differences of observation:
  
  var_sc_prereq[,sq_m_diff := (dt[,SST_hat][id_1]-dt[,SST_hat][id_2])^2]
  
  var_sc_prereq[,sq_obs_diff := (dt[,SST_bar][id_1]-dt[,SST_bar][id_2])^2]
  
  
  # get variogram scores
  print("done - compute variogram scores:")
  scores = list()
  
  #var_sc_by_PC = function(d){
  list_page = 1
  
  for(d in dvec){
    print(paste0("d = ",d))
    return_data =  data.table(year = y, d = d)
    
    dt_temp = var_sc_prereq[,.SD,.SDcols = c("sq_m_diff",paste0("dvar",d),"sq_obs_diff","weights")]
    setnames(dt_temp,c("sq_mean_diff", "dvar", "sq_obs_diff","weights"))
    return_data[,sc := variogram_score_nrm_p2(dt_temp)]
    
    
    scores[[list_page]] = return_data[,month := m]  
    list_page = list_page +1
  }
  
  scores = rbindlist(scores)
  
  if(!marginal_correction)
  {
    file_name = paste0(file_name,"_nmc")
  }
  
  save(scores, file = paste0(save_dir,file_name,"_m",m,"_y",y,".RData"))
}
  


####################


