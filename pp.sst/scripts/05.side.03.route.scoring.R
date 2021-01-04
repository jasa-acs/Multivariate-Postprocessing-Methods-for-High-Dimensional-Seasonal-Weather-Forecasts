########################################################################################################

#############  side script 5.3 - score various functionals that depend on areas of SST  ################

########################################################################################################

# We score max and min SST over a predefined areal. 
# The example contained in the paper is minimum SST along the shipping route Norfolk to Bordeaux.
# Take care that the considered area is actually contained in DT for your current run!

#### setting up ######

rm(list = ls())

options(max.print = 1e3)

library(pp.sst)

name_abbr = "Full" 

save_dir = file.path('~','SST','Derived', name_abbr)

load(file = paste0(save_dir,"setup.RData"))

time_s53 = proc.time()

####################################
####### area specification #########
####################################

# specify the area to score over. Version 1: specify a route:

# Fix two points p1 and p2. These are then connected by a route (as the crow flies, not caring about land in between).
# A grid point in dt is then considered for scoring if it is not on land and is nearest neighbour to a point on the route. 

Bordeaux = c(-0.57,44.8)
Norfolk = c(-76.3,36.9)

p1 = data.table(Lon = Bordeaux[1], Lat = Bordeaux[2], Loc = 'Bordeaux') 
p2 = data.table(Lon = Norfolk[1], Lat = Norfolk[2], Loc = "Norfolk") 

route_name = "Bordeaux to Norfolk"
file_name = "scores_Bd_to_Nf"


# p1 = data.table(Lon = 5.32, Lat = 60.4, Loc = "Bergen") 
# p2 = data.table(Lon = -21.83, Lat = 64.13, Loc = "Reykyavik") 
# 
# route_name = "Bergen to Reykyavik"
# file_name = "scores_Bergen_to_Reykyavik"
# 
# p1 = data.table(Lon = -46.5, Lat = -23.8, Loc = "Sao Paolo") 
# p2 = data.table(Lon = 18.6, Lat = -34.5, Loc = "Capetown") 
# 
# route_name = "Sao Paolo to Capetown"
# file_name = "scores_SP_to_Capetown"


# get grid ids_along this route: fix n and use gcIntermediate to find n coordinates on the route from p1 to p2 as the crow flies:

n = 1000
route = geosphere::gcIntermediate(p1[,.(Lon,Lat)],p2[,.(Lon,Lat)],n = n)

# find the grid_ids in DT closest to the coordinates on the route:

grid_id_dt = unique(DT[,.(Lon,Lat,grid_id)])

point_match = NULL
for(j in 1:dim(route)[1])
{
  a = geosphere::distHaversine(as.vector(route[j,]),as.matrix(grid_id_dt[,.(Lon,Lat)]))
  point_match[j] = which.min(a)
}

dt_route = unique(grid_id_dt[point_match,])
route_ids = dt_route[,grid_id]

###################################

# specify the area to score over. Version 2: specify a window:

# option two: Lon/Lat window

# Lon_min = -25
# Lon_max = -15
# Lat_min = 60
# Lat_max = 64
# 
# grid_id_dt = unique(DT[Lon >= Lon_min & Lon <= Lon_max,][Lat >= Lat_min & Lat <= Lat_max,.(Lon,Lat,grid_id)])
# route_ids = grid_id_dt[,grid_id]
# 
# route_name = "South of Island"
# file_name = "scores_south_of_Island"


######################################

# get forecasts

load(paste0(PCA_dir,"fc_mc.RData"))
load(paste0(PCA_dir,"fc_ac.RData"))
load(paste0(GS_dir,"fc.RData"))
load(paste0(ECC_dir,"fc.RData"))
load(paste0(Schaake_dir,"fc.RData"))

PCA_ac_fc = PCA_fc_ac
PCA_mc_fc = PCA_fc_mc

########################################

mod_vec = c('PCA_mc','PCA_ac','GS','ECC','Schaake')

mod_vec_all = c('PFC',mod_vec)

# reduce data to route

DT_route = DT[grid_id %in% route_ids & year %in% validation_years,]

for(mod in mod_vec){
  assign(paste0(mod,'_fc_route'),get(paste0(mod,'_fc'))[grid_id %in% route_ids,])
}


rm(PCA_mc_fc,PCA_fc_mc,PCA_ac_fc,PCA_fc_ac,GS_fc,ECC_fc,Schaake_fc)

########################################

# put considered functionals in a list

fun_nlist = list('max','min','mean')
fun_list = lapply(fun_nlist,get)


######################################################
###### Get MSEs for the  maximum along route ########
######################################################

# initialize data tables:

validation_dt = as.data.table(expand.grid(months,validation_years))
setnames(validation_dt,c("month","year"))

scores_dt = as.data.table(expand.grid(  model = mod_vec_all, fun = fun_nlist,MSE = 0,CRPS = 0))

# get observation

for(name in fun_nlist)
{
  fun = get(name)
  
  DT_route[,paste0('SST_',name) := fun(SST_bar, na.rm = TRUE), by = .(month,year)]
  
  temp = DT_route[grid_id == min(grid_id) ,eval(parse(text = paste0('SST_',name)))]
  validation_dt[,paste0('SST_',name) := temp]
  
  # get value of point forecast
  
  # the weird column name is for unification with other models
  DT_route[,paste0('PFC_',name,'_fc') := fun(SST_hat,na.rm = TRUE),by = .(month,year)]
  
  temp = DT_route[grid_id == min(grid_id) ,eval(parse(text = paste0('PFC_',name,'_fc')))]
  validation_dt[,paste0('PFC_',name,'_fc') := temp]
  
  # get value of multivariate forecasts for the different models:
  
  for(mod in mod_vec)
  {
    print(c(name,mod))
    
    # how large is the ensemble
    es_mod = fc_ens_size 
    if( mod == 'ECC')
    {
      es_mod = 9
    }
    if( mod == 'Schaake')
    {
      es_mod = length(training_years)
    }
    
      
    
    fc_temp = get(paste0(mod,'_fc_route'))
    
    # get the forecasted values:
    for(i in 1:es_mod)
    {
      fc_temp[,paste0(mod,'_',name,'_',i) := lapply(X = .SD,FUN = fun,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
    }
    
    # take mean over all forecasts and put this into validation_dt
    
    fc_temp = fc_temp[grid_id == min(grid_id) ,.SD,.SDcols = paste0(mod,'_',name,'_',1:es_mod)]
    fc_temp[,paste0(mod,'_',name,'_fc') := rowMeans(.SD),.SDcols = paste0(mod,'_',name,'_',1:es_mod)]
    
    temp = fc_temp[,.SD,.SDcols = paste0(mod,'_',name,'_fc')]
    
    validation_dt = data.table(validation_dt,temp)
    
    ### get MSE ###
    
    obs = validation_dt[,eval(parse(text = paste0('SST_',name)))]
    fcs = validation_dt[,eval(parse(text = paste0(mod,'_',name,'_fc')))]
    
    validation_dt[,paste0('MSE_',mod,'_',name,'_fc'):=(obs-fcs)^2]
    
    mse = mean((obs - fcs)^2)
    
    scores_dt[model == mod & fun == name,MSE:=mse]
    
    ### get CRPS ###
    
    data = as.matrix(fc_temp[,.SD,.SDcols = paste0(mod,'_',name,'_',1:es_mod)])
    
    crps = scoringRules::crps_sample(y = obs, 
                                     dat = data)
    
    validation_dt[,paste0('CRPS_',mod,'_',name,'_fc'):=crps]
    
    scores_dt[model == mod & fun == name, CRPS := mean(crps)]
  }

  # scores for point forecasts:
  mod = 'PFC'
  
  obs = validation_dt[,eval(parse(text = paste0('SST_',name)))]
  fcs = validation_dt[,eval(parse(text = paste0(mod,'_',name,'_fc')))]
  
  mse = mean((obs - fcs)^2)
  crps = mean(abs(obs-fcs))
  
  scores_dt[model == mod & fun == name,MSE:=mse]
  scores_dt[model == mod & fun == name,CRPS:=crps]
}
  

save(scores_dt,file = paste0(save_dir,'route_scores.RData'))

### permutation tests ###

funn = 'min'
score = 'MSE'

N = 20000

#pdf(paste0(plot_dir,'perm_test_routescores.pdf'),width = 14,height = 14)

par( oma=c(0,6,6,0), mfrow=c(length(mod_vec),length(mod_vec)), mar=c(3,3,2,2)+0.1 )

for(mod1 in mod_vec)
{
  
  for(mod2 in mod_vec)
  {
    col_name_1 = paste0(score,'_',mod1,'_',funn,'_fc')
    col_name_2 = paste0(score,'_',mod2,'_',funn,'_fc')
    
        # permutation test for mod1 ~ mod2
    pt_vs = permutation_test_difference(validation_dt[,get(col_name_1)],validation_dt[,get(col_name_2)], N=N)
    
    x_lim_max = 1.1*max(abs(c(pt_vs$d_bar,pt_vs$D)))
    
    hist(pt_vs$D, xlim = c(-x_lim_max,x_lim_max),main = paste0(mod1,' vs. ',mod2),breaks = 20)
    
    qq_1 = quantile(pt_vs$D,c(0.05))
    
    abline(v = pt_vs$d_bar,col = 'red')
    
    abline(v = qq_1,lty = 2)
    
    qq_2 = quantile(pt_vs$D,c(0.01))
    
    abline(v = qq_2,lty = 3)
  }
}


### get p-value for all tests ###

p_values_route_scoring = as.data.table(expand.grid(MOD1 = mod_vec,MOD2 = mod_vec,SCORE = c('MSE','CRPS'),FUN = c('max','min')))

p_values_route_scoring[,p_value:=NA_real_]

for(funn in c('max','min') )
{
  for(score in c('CRPS','MSE'))
  {
    print(c(funn,score))
    
    for(mod1 in mod_vec)
    {
      
      for(mod2 in mod_vec)
      {
        col_name_1 = paste0(score,'_',mod1,'_',funn,'_fc')
        col_name_2 = paste0(score,'_',mod2,'_',funn,'_fc')
        
        # permutation test for mod1 ~ mod2
        pt_vs = permutation_test_difference(validation_dt[,get(col_name_1)],validation_dt[,get(col_name_2)], N=N)
        
        p_val = pt_vs$p_val
        
        p_values_route_scoring[MOD1 == mod1 & MOD2 == mod2 & SCORE == score & FUN == funn, p_value := p_val]
        
      }
    }
    
  }
}

save(p_values_route_scoring,file = paste0(save_dir,'p_values.RData'))

#dev.off()



save.image(file = paste0(save_dir,'setup.RData'))

  


# ##################################
# ############# plots ##############
# ##################################
# 
# # max
# pdf(paste0(plot_dir,"RMSE_max_along_route.pdf"))
# 
# rr = range(RMSEs[,.SD,.SDcols = c("PCA_max","GS_max","ECC_max")])
# 
# plot(RMSEs[,PCA_max],type = "b", col = "blue",
#      xlab = "route", ylab = "RMSE", main = "RMSE of max. temp. on a route",
#      ylim = rr)
# lines(RMSEs[,ECC_max],type = "b", col = "darkgreen")
# lines(RMSEs[,GS_max],type = "b", col = "darkred")
# 
# legend(x = "topright", legend = c("PCA","GS","ECC"),col = c("blue","darkred","darkgreen"),lty = c(1,1,1))
# dev.off()
# 
# # min
# pdf(paste0(plot_dir,"RMSE_min_along_route.pdf"))
# 
# rr = range(RMSEs[,.SD,.SDcols = c("PCA_min","GS_min","ECC_min")])
# 
# plot(RMSEs[,PCA_min],type = "b", col = "blue",
#      xlab = "route", ylab = "RMSE", main = "RMSE of min. temp. on a route",
#      ylim = rr)
# lines(RMSEs[,ECC_min],type = "b", col = "darkgreen")
# lines(RMSEs[,GS_min],type = "b", col = "darkred")
# 
# legend(x = "topright", legend = c("PCA","GS","ECC"),col = c("blue","darkred","darkgreen"),lty = c(1,1,1))
# dev.off()
# 
# 
# # mean
# pdf(paste0(plot_dir,"RMSE_mean_along_route.pdf"))
# 
# rr = range(RMSEs[,.SD,.SDcols = c("PCA_mean","GS_mean","ECC_mean")])
# 
# plot(RMSEs[,PCA_mean],type = "b", col = "blue",
#      xlab = "route", ylab = "RMSE", main = "RMSE of mean temp. on a route",
#      ylim = rr)
# lines(RMSEs[,ECC_mean],type = "b", col = "darkgreen")
# lines(RMSEs[,GS_mean],type = "b", col = "darkred")
# 
# legend(x = "topright", legend = c("PCA","GS","ECC"),col = c("blue","darkred","darkgreen"),lty = c(1,1,1))
# dev.off()
