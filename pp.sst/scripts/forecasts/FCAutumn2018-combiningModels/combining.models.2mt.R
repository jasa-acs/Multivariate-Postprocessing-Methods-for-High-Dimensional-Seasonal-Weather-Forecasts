
# surface temperature forecast averaged over east and west Norway


# setup

rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)


# get NorCPM temperature forecast and the regions of Norway

data_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/ts/")

load(file = paste0(data_dir,"dt_combine_mr_wide.RData"))
load(file = paste0(data_dir,"../../GCFS1/DT_2mt.RData"))
load(file = paste0(data_dir,"../Norway_regions.RData"))

# reduce data to ost- and vestlandet

id_ol = Norway_regions[ind_nor == 2,grid_id]
id_vl = Norway_regions[ind_nor == 1,grid_id]

#select relevant months

mon_sel = 8:10

#### gettin YM and grid_ids for DT ####

DT[,YM := 12 * year + month]

DT = DT[order(year,month,Lon,Lat)]
dt = dt[order(year,month,Lon,Lat)]

gids = rep(dt[year == 1985 & month == 1, grid_id],length(unique(DT[,YM])))
DT[,grid_id:= gids]

dt[,YM := 12 * year + month]

################

dt_small = dt[grid_id %in% c(id_ol,id_vl) & YM %in% DT[,unique(YM)] ,]

dt_sm_GCFS = DT[grid_id %in% c(id_ol,id_vl) & YM %in% dt_small[,unique(YM)],]

st_mean = merge(unique(dt_small[,.(year,month)])[,region := "ol"],unique(dt_small[,.(year,month)])[,region := "vl"],all = TRUE)
st_mean[, ym := 12*year+month]

st_mean = st_mean[month %in% mon_sel]

# get means NorCPM
st_ol = dt_small[grid_id %in% id_ol & month %in% mon_sel, mean(ts_bar), by = .(year,month)][,3]
st_vl = dt_small[grid_id %in% id_vl & month %in% mon_sel, mean(ts_bar), by = .(year,month)][,3]

st_mean[region == "ol",st_NCPM := st_ol]
st_mean[region == "vl",st_NCPM := st_vl]

# get means GCFS
st_ol = dt_sm_GCFS[grid_id %in% id_ol & month %in% mon_sel, mean(Ens_bar), by = .(year,month)][,3]
st_vl = dt_sm_GCFS[grid_id %in% id_vl & month %in% mon_sel, mean(Ens_bar), by = .(year,month)][,3]

st_mean[region == 'ol',st_GCFS := st_ol]
st_mean[region == 'vl',st_GCFS := st_vl]

###############

# getting observations:

vl_obs_dt = as.data.table(read.csv(file = paste0("~/PostClimDataNoBackup/SFE/ST_mean_Nor/mean_st_vl.csv")))
ol_obs_dt = as.data.table(read.csv(file = paste0("~/PostClimDataNoBackup/SFE/ST_mean_Nor/mean_st_ol.csv")))

# replace months by number, add ym

month_names = vl_obs_dt[1:12,month]
m = match(vl_obs_dt[,month],month_names)

vl_obs_dt[,month:=m]
ol_obs_dt[,month:=m]

vl_obs_dt[, ym := 12*year+month]
ol_obs_dt[, ym := 12*year+month]

vl_obs_dt = vl_obs_dt[ym %in% st_mean[,ym],]
ol_obs_dt = ol_obs_dt[ym %in% st_mean[,ym],]

# plugging into st_mean

vl_obs = c(vl_obs_dt[,st],rep(NA,st_mean[region == "vl",.N]-vl_obs_dt[,.N]))
ol_obs = c(ol_obs_dt[,st],rep(NA,st_mean[region == "ol",.N]-ol_obs_dt[,.N]))

st_mean[region == "vl",t2m := vl_obs]
st_mean[region == 'ol',t2m := ol_obs]


############################# 

###### post-processing ######

#############################

training_years = 1985:2009
validation_years = 2010:2017

# bias correction:

st_mean[,b_hat_NCPM := sim_mov_av(l = 50,vec = st_NCPM - t2m,years = year),by = .(month,region)]
st_mean[,b_hat_GCFS := sim_mov_av(l = 50,vec = st_GCFS - t2m,years = year),by = .(month,region)]

st_mean[,st_hat_NCPM := st_NCPM - b_hat_NCPM]
st_mean[,st_hat_GCFS := st_GCFS - b_hat_GCFS]

# climatology

st_mean[year %in% training_years, clim := mean(t2m), by = .(month,region) ]

st_mean[year %in% validation_years, clim := rep(st_mean[year == min(year),clim],length(validation_years))]

st_mean[year == 2018 , clim := st_mean[year == min(year) & month %in% mon_sel, clim]]

st_mean = st_mean[order(year,month,region)]




###########################
##### mixture models ######
###########################

n = 100
lambda_vec = seq(0,1,length.out = 100)

RMSE_vec = c()

for(lambda_1 in lambda_vec)
{
  lambda_2 = 1-lambda_1
  
  # get a data table with mixed forecasts:
  mixed_fc = lambda_1*st_mean[,st_hat_NCPM]+lambda_2 * st_mean[,st_hat_GCFS]
  
  DT_m = st_mean[,.(year,month,region,ym,t2m)]
  DT_m[,st_hat := mixed_fc]
  
  RMSE_mixed = sqrt(DT_m[year %in% validation_years,mean( (st_hat - t2m)^2)])
  RMSE_vec = c(RMSE_vec,RMSE_mixed)
}

plot(lambda_vec,RMSE_vec,type = "l")

lambda_1 = lambda_vec[which(RMSE_vec == min(RMSE_vec))]
lambda_2 = 1-lambda_1

st_mean[,st_hat_mm_is := lambda_1 * st_hat_NCPM + lambda_2 * st_hat_GCFS]



#tryout out of sample for comparison

validation_years_2 = 1995:2010


# estimate lambda by month

RMSE_mm = data.table(month =integer(),
                     lambda_NCPM = numeric(),
                     lambda_GCFS = numeric(),
                     lambda_clim = numeric(),
                     RMSE = numeric())


for(mon in 8:11)
{ print(mon)
  lambda_est = function(i)
  {
    lambda_1 = lambda_vec[i]
    RMSE_mm = data.table(month = integer(),
                         lambda_NCPM = numeric(),
                         lambda_GCFS = numeric(),
                         lambda_clim = numeric(),
                         RMSE = numeric())
    for(j in 0:(100-i))
    {
      
      lambda_2 = lambda_vec[j+i] - lambda_1
      lambda_3 = 1 - lambda_1 - lambda_2
      
      # get a data table with mixed forecasts:
      mixed_fc = lambda_1*st_mean[,st_hat_NCPM]+lambda_2 * st_mean[,st_hat_GCFS] + lambda_3 * st_mean[,clim]
      
      DT_m = st_mean[,.(year,month,region,ym,t2m)]
      DT_m[,st_hat := mixed_fc]
      
      RMSE = sqrt(DT_m[year %in% validation_years_2 & month == mon , mean( (st_hat - t2m)^2)])
      RMSE_mm = rbindlist(list(RMSE_mm,
                               data.table(month = mon,
                                          lambda_NCPM = lambda_1,
                                          lambda_GCFS = lambda_2,
                                          lambda_clim = lambda_3,
                                          RMSE = RMSE)))
    }
    return(RMSE_mm)
  }
  
  RMSE_mm = rbindlist(list(RMSE_mm,rbindlist(parallel::mclapply(1:100,lambda_est,mc.cores = 10))))
}


RMSE_mm_min = RMSE_mm[0,]
for(m in months)
{
  RMSE_mm_min = rbindlist(list(RMSE_mm_min,RMSE_mm[month == m,][RMSE == min(RMSE),]))
  lambda_NCPM = RMSE_mm[month == m,][RMSE == min(RMSE),lambda_NCPM]
  lambda_GCFS = RMSE_mm[month == m,][RMSE == min(RMSE),lambda_GCFS]
  lambda_clim = RMSE_mm[month == m,][RMSE == min(RMSE),lambda_clim]
  st_mean[month == m, st_hat_mm_oos := lambda_NCPM * st_hat_NCPM + lambda_GCFS * st_hat_GCFS + lambda_clim * clim]
}




#### MSEs, RMSEs, and boxplots  ###

MSEs = st_mean[year %in% validation_years,.(year,month,region)]


for(r in c("ol","vl"))
{
  MSEs[region == r,NCPM := st_mean[year %in% validation_years & region == r, (st_hat_NCPM - t2m)^2]]
  MSEs[region == r,GCFS := st_mean[year %in% validation_years & region == r, (st_hat_GCFS - t2m)^2]]
  MSEs[region == r,clim := st_mean[year %in% validation_years & region == r, (clim - t2m)^2]]
  MSEs[region == r,mm_is := st_mean[year %in% validation_years & region == r, (st_hat_mm_is - t2m)^2]]
  MSEs[region == r,mm_oos := st_mean[year %in% validation_years & region == r, (st_hat_mm_oos - t2m)^2]]
  
}



#boxplots
plot_dir = paste0("./figures/Aut_2018_small/")

pdf(paste0(plot_dir,"bp_ol.pdf"))
boxplot(MSEs[region == "ol",.(sqrt(NCPM),sqrt(GCFS),sqrt(mm_oos),sqrt(clim))],names = c("NCPM","GCFS1","mm","clim"),main = "squared error ostlandet")
dev.off()

pdf(paste0(plot_dir,"bp_vl.pdf"))
boxplot(MSEs[region == "vl",.(sqrt(NCPM),sqrt(GCFS),sqrt(mm_oos),sqrt(clim))],names = c("NCPM","GCFS1","mm","clim"),main = "squared error vestlandet")
dev.off()
#RMSEs 

models = c("clim","NCPM","GCFS","mm_is","mm_oos")
regions = c("ol","vl")

RMSEs = as.data.table(expand.grid(models,regions))
setnames(RMSEs,c("model","region"))
RMSEs = RMSEs[order(model,region)]

for(m in models)
{
  for(r in regions)
  {
    RMSEs[region == r & model == m, RMSE := sqrt(MSEs[region == r,lapply(.SD,mean) ,.SDcols = m])]                         
  }
}



#RMSEs by months

models = c("clim","NCPM","GCFS","mm_is","mm_oos")
regions = c("ol","vl")
months = mon_sel

RMSEs_by_month = as.data.table(expand.grid(models,regions,months))
setnames(RMSEs_by_month,c("model","region","month"))
RMSEs_by_month = RMSEs_by_month[order(model,region)]

for(m in models)
{
  for(r in regions)
  {
    for(mon in months)
    {
      RMSEs_by_month[region == r & model == m & month == mon, RMSE := sqrt(MSEs[region == r & month == mon,lapply(.SD,mean) ,.SDcols = m])]                         
    }
    
  }
}


# plot
rr = range(RMSEs_by_month[region == "ol",RMSE])
plot(RMSEs_by_month[model == "clim" & region == "ol",.(month,RMSE)],type = "b",ylim = rr,xlab = 'RMSE',main = 'RMSE by month for ostlandet')
lines(RMSEs_by_month[model == "NCPM" & region == "ol",.(month,RMSE)],type = "b",col = 'blue')
lines(RMSEs_by_month[model == "GCFS" & region == "ol",.(month,RMSE)],type = "b",col = 'green')
lines(RMSEs_by_month[model == "mm_oos" & region == "ol",.(month,RMSE)],type = "b",col = 'red')


rr = range(RMSEs_by_month[region == "vl",RMSE])
plot(RMSEs_by_month[model == "clim" & region == "vl",.(month,RMSE)],type = "b",ylim = rr,xlab = 'RMSE',main = 'RMSE by month for ostlandet')
lines(RMSEs_by_month[model == "NCPM" & region == "vl",.(month,RMSE)],type = "b",col = 'blue')
lines(RMSEs_by_month[model == "GCFS" & region == "vl",.(month,RMSE)],type = "b",col = 'green')
lines(RMSEs_by_month[model == "mm_oos" & region == "vl",.(month,RMSE)],type = "b",col = 'red')

# RMSE combined

RMSE_combined = RMSEs[, .(sqrt(mean(RMSE^2))),by = model]
setnames(RMSE_combined,"V1","RMSE")



# get forecasts

forecasts = st_mean[year == 2018 & month %in% 8:10,.(month,region,st_hat_NCPM,st_hat_GCFS,st_hat_mm_oos,clim)]
forecasts[,ano_NCPM := st_hat_NCPM - clim][,ano_GCFS := st_hat_GCFS - clim][,ano_mm := st_hat_mm_oos - clim]

###################################



