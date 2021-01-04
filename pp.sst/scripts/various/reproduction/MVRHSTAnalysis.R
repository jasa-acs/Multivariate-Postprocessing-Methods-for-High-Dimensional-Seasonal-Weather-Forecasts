# figure out whats wrong with the route scoring


rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "NAO/lv/2" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))

time_s54 = proc.time()

mc_cores = 5

y = 2003

par(mfrow = c(3,4),oma = c(0,0,1,0))

for(m in 1:12){
  test = matrixStats::rowRanks(as.matrix(GS_fc_route[year == y &month == m & !is.na(SST_hat),.SD,.SDcols = c('SST_bar',paste0('fc',1:500))]))
  hist(test[,1],main = c(m,y))
}
title('GS',outer = TRUE)


for(m in 1:12){
  test2 = matrixStats::rowRanks(as.matrix(PCA_mc_fc_route[ year == y & month == m & !is.na(SST_hat),.SD,.SDcols = c('SST_bar',paste0('fc',1:500))]))
  hist(test2[,1],main = c(m,y))
}
title('PCA_mc',outer = TRUE)

for(m in 1:12){
  test3 = matrixStats::rowRanks(as.matrix(PCA_ac_fc_route[year == y & month == m & !is.na(SST_hat),.SD,.SDcols = c('SST_bar',paste0('fc',1:500))]))
  hist(test3[,1],main = c(m,y))
}
title('PCA_ac',outer = TRUE)

for(m in 1:12)
{
  test4 = matrixStats::rowRanks(as.matrix(ECC_fc_route[year == y & month == m & !is.na(SST_hat),.SD,.SDcols = c('SST_bar',paste0('fc',1:9))]))
  hist(test4[,1],main = c(m,y))
}
title('ECC',outer = TRUE)



#all together


par(mfrow = c(2,2),oma = c(0,0,1,0))

test = matrixStats::rowRanks(as.matrix(GS_fc_route[ !is.na(SST_hat),.SD,.SDcols = c('SST_bar',paste0('fc',1:500))]))
hist(test[,1],main = '')
title('GS',outer = FALSE)


  test2 = matrixStats::rowRanks(as.matrix(PCA_mc_fc_route[ !is.na(SST_hat),.SD,.SDcols = c('SST_bar',paste0('fc',1:500))]))
  hist(test2[,1],main = '')
title('PCA_mc',outer = FALSE)


  test3 = matrixStats::rowRanks(as.matrix(PCA_ac_fc_route[ !is.na(SST_hat),.SD,.SDcols = c('SST_bar',paste0('fc',1:500))]))
  hist(test3[,1],main = '')
  title('PCA_ac',outer = FALSE)

  test4 = matrixStats::rowRanks(as.matrix(ECC_fc_route[ !is.na(SST_hat),.SD,.SDcols = c('SST_bar',paste0('fc',1:9))]))
  hist(test4[,1],main = '')
  title('ECC',outer = FALSE)

  
 ### # case study for band depth rank###
  
  N = 100 # size of ensemble
  
  d = 100 # dimension
  
  rep = 500 # number of fc - obs - pairs
  
  par(mfrow = c(1,2))
  
  # case 1: perfect correlation in the forecast
  
  offset = -0.25
  
  offset_sd = 0.1
  
  sd_var = 1
  
  
  band_depth_rk_obs = c()
  for(i in 1:rep)
  {
    if(i %% 50 == 0) print(i)
    # observation: rep repetitions of a standard normal, perfectly correlated along dimension
    
    obs = rep(rnorm(n=1, mean = 0, sd = 1), d)
    
    # fc model 1 
    
    fc = rep(rnorm(n=N, mean = offset, sd = sd_var + offset_sd),each = d)
    
    ###
    data =matrix(c(obs,fc),nrow = d)  
    
    band_depth_rk_obs = c(band_depth_rk_obs,bd.rank(data)[1])
    
  }
  
  
  hist(band_depth_rk_obs, main = 'right correlation')
  
  # model 2: 0 correlation along d
  
  band_depth_rk_obs = c()
  for(i in 1:rep)
  {
    if(i %% 50 == 0) print(i)
    # observation: rep repetitions of a standard normal, perfectly correlated along dimension
    
    obs = rep(rnorm(n=1, mean = 0, sd = 1 + offset_sd), d)
    
    # fc model 1 
    
    fc = rnorm(n=N*d, mean = offset, sd = sd_var)
    
    ###
    data =matrix(c(obs,fc),nrow = d)  
    
    band_depth_rk_obs = c(band_depth_rk_obs,bd.rank(data)[1])
    
  }
  
  
  hist(band_depth_rk_obs, main = ' too weak correlation')
  

  
  
  # check how this ranks ties
  
  data = matrix(1,nrow = d,ncol = N)
  a = bd.rank(data)
