
# fooling around with weight functions:

PCs = c(0:10,15,20,30,40)


weight_fct = function(x)
{y = rep(0,length(x))
  y[abs(x)<1 & abs(x)>0.00001]=1/abs(x[abs(x)<1 & abs(x)>0.00001])
  y[abs(x)<0.00001] = 1/0.00001
return(y)
}


months = 1
validation_years = 2001:2003


for(m in months)
{ 
  # setup: get principal components and marginal variances for the given month m:
  
  print(paste0("m = ",m))
  load(file = paste0(PCA_dir,"CovRes_mon",m,".RData"))
  PCA = irlba::irlba(res_cov, nv = max(PCs))
  
  land_ids = which(DT[year == min(year) & month == min(month),is.na(Ens_bar) | is.na(SST_bar)])
  
  PCA_DT = DT[year == min(year) & month == min(month),][-land_ids,.(Lon,Lat,grid_id)]
  
  for(d in  min(PCs):max(PCs))
  {
    if(d ==0) PCA_DT [,paste0("PC",d) := 0]
    if(d>0)
    {
      PCA_DT [,paste0("PC",d) := PCA$d[d]*PCA$u[,d]]
    }
  } 
  
  # also get marginal variances
  variances = list()
  d = min(PCs[PCs>0]) 
  vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
  variances[[d]] = vec^2
  
  if(length(PCs)>1)
  {
    for(d in (min(PCs[PCs>0])+1):max(PCs)){
      vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
      variances[[d]] = variances[[d-1]] + vec^2
    }
  }
  names(variances) = paste0("var",min(PCs[PCs>0]):max(PCs))
  PCA_DT = data.table(PCA_DT,rbindlist(list(variances)))    
  PCA_DT[,var0:=0]
  
  # now move to getting the variogram scores:
  
  # without marginal correction:
  
  dummy_fct = function(y)
  {
    var_sc_PCA(m, y, DT, PCA = PCA, PCA_DT = PCA_DT, dvec = PCs, ens_size = ens_size,
               marginal_correction = FALSE, 
               cov_dir = PCA_dir, save_dir = PCA_dir,weight_fct = weight_fct)
  }
  parallel::mclapply(validation_years,dummy_fct,mc.cores = length(validation_years))
  
  # with marginal correction:
  
  dummy_fct = function(y)
  {
    print(c(paste0("y = ",y),paste0("m = ",m)))
    var_sc_PCA(m, y, DT, PCA = PCA, PCA_DT = PCA_DT, dvec = PCs,ens_size = ens_size,
               cov_dir = PCA_dir,save_dir = PCA_dir,weight_fct = weight_fct)
  }
  parallel::mclapply(validation_years,dummy_fct,mc.cores = length(validation_years))
}



###### combine: ######

scores_nmc = list()
k=0
for(m in months){
  for(y in validation_years)
  {k=k+1
  load(file = paste0(PCA_dir,"var_sc_by_PC_nmc_m",m,"_y",y,".RData"))
  scores_nmc[[k]] = scores
  }
}
scores_nmc = rbindlist(scores_nmc)

scores_mc = list()
k=0
for(m in months){
  for(y in validation_years)
  {k=k+1
  load(file = paste0(PCA_dir,"var_sc_by_PC_m",m,"_y",y,".RData"))
  scores_mc[[k]] = scores
  }
}
scores_mc = rbindlist(scores_mc)


#############################################

mean_sc = scores_mc[,mean_sc := mean(sc), by = d][,.(d,mean_sc)]
mean_sc = unique(mean_sc)
mean_sc = mean_sc[d>0,]


mean_sc_nmc = scores_nmc[,mean_sc := mean(sc), by = d][,.(d,mean_sc)]
mean_sc_nmc = unique(mean_sc_nmc)

no_pp_score = mean_sc_nmc[d == 0, mean_sc]

mean_sc_nmc = mean_sc_nmc[d>0,]

print(paste0("Minimal variogram score is achieved for ",mean_sc[mean_sc == min(mean_sc),d][1],
             " principal components with a score of ",mean_sc[,min(mean_sc)],
             ". Without marginal correction it is achieved for ",mean_sc_nmc[mean_sc == min(mean_sc),d][1],
             " principal components with a score of ",mean_sc_nmc[,min(mean_sc)],"."))


rr = range(c(mean_sc[[2]],mean_sc_nmc[[2]],no_pp_score))


#### with marginal correction ####

plot(x = mean_sc[[1]],
     y = mean_sc[[2]],
     ylim = rr,
     type = "b",
     col = "blue",
     main = paste0("mean variogram scores for ",name_abbr),
     xlab = "number of principal components",
     ylab = "mean score")

#---- add minima: -----
abline(h = no_pp_score, lty = "dashed", col = adjustcolor("blue",alpha = .5))

opt_num_PCs = mean_sc[,which.min(mean_sc)]

points(x = mean_sc[[1]][opt_num_PCs],
       y = mean_sc[[2]][opt_num_PCs],
       col = "blue",
       bg = "blue",
       pch = 21)

#### without marginal correction ####

lines(x = mean_sc_nmc[[1]],
      y = mean_sc_nmc[[2]],
      type = "b",
      col = "darkgreen")

#---- add minima: -----

opt_num_PCs = mean_sc_nmc[,which.min(mean_sc)]

points(x = mean_sc_nmc[[1]][opt_num_PCs],
       y = mean_sc_nmc[[2]][opt_num_PCs],
       col = "darkgreen",
       bg = "darkgreen",
       pch = 21)

