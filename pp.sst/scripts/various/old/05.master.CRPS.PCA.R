
#######################################################################################

########  master script part 5 - CRPS for PCA without marginal correction  ############

#######################################################################################

# This script compares the marginal performance of the PCA method with and without marginal correction
#
# 
# Files generated:
#   
#
# Requires previous run of 03.master.var.est 
# with the same value of name_abbr as below.

##### setting up ######

rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "NAO_2" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))
DT = load_combined_wide(data_dir = save_dir, output_name = paste0("dt_combine_wide_bc_var.RData"))

###### getting marginal CRPS averaged over all ensemble members  ######

PCA_dir = paste0(save_dir,"PCA_new/")

months = 1:12

PCs = 1:50



setnames(mean_crps_by_pc,c("month","PCs"))


dummy_fct = function(m)
{ 
  mean_crps_by_pc = data.table("PCs" = PCs)
  #get covariance matrix
  load(file = paste0(PCA_dir,"CovRes_mon",m,".RData"))
  
  PCA = irlba::irlba(res_cov, nv = max(PCs))
  
  land_ids <- which(DT[year == min(year) & month == min(month), is.na(Ens_bar) | is.na(SST_bar)])  
  
  PCA_DT = DT[year == min(year) & month == min(month),][-land_ids,.(Lon,Lat,grid_id)]
  
  for(d in  PCs){
    PCA_DT [,paste0("PC",d) := PCA$d[d]*PCA$u[,d]]/sqrt(ens_size)
  } 
  
  # also get marginal variances
  variances = list()
  d = min(PCs) 
  vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
  variances[[d]] = vec^2
  
  for(d in PCs[2:length(PCs)]){
    vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
    variances[[d]] = variances[[d-1]] + vec^2
  }
  names(variances) = paste0("var",PCs)
  PCA_DT = data.table(PCA_DT,rbindlist(list(variances)))  
  
  # get marginal CRPS by number of principal component:
  
  crps_by_pc = matrix(0,nrow = PCA_DT[,.N],ncol = length(PCs))
  for(y in validation_years)
  {
    print(c(paste0("month = ",m),paste0("year = ",y)))
    obs = DT[year == y & month == m,][-land_ids,SST_bar]
    for(k in 1:ens_size){
        est = DT[year == y & month == m,][-land_ids,trc(eval(parse(text = paste0("Ens",k))) + Bias_Est)]
        for(d in PCs)
        {
          sd = PCA_DT[,sqrt(eval(parse(text = paste0("var",d)))) ]
          crps_by_pc[,d] = crps_by_pc[,d] + crps.na.rm(obs,est,sd)/(length(validation_years)*ens_size)
        }
    }
  }
  
  mean_crps_by_pc[,"CRPS" := colMeans(crps_by_pc,na.rm = TRUE)]
  return(mean_crps_by_pc)
  
}

mcrps = parallel::mclapply(1:12,dummy_fct,mc.cores = 12)  

mcrps = rbindlist(mcrps)[,mCRPS := mean(CRPS), by = PCs]
mcrps[PCs == 1,]

mcrps[,CRPS:=NULL]
mcrps= mcrps[unique(PCs),.(PCs,mCRPS)]


#### plotting ####

# get CRPS for marginally post-processed PCA:

load(paste0(save_dir,"/scores.bc.sd.sma.Rdata"))
load(paste0(save_dir,"/scores.bc.sd.ema.Rdata"))

CRPS_min = min(sc_sma_var[,CRPS],sc_ema_var[,CRPS])

y_range = range(c(CRPS_min,mcrps[,mCRPS]))

pdf(paste0(plot_dir,"/mean_CRPS_PCA_moem.pdf"))
plot(x = PCs,
     y = mcrps[,unique(mCRPS)],
     ylim = y_range,
     type = "l",
     col = "blue",
     main = paste0("CRPS by PC without variance correction, MOEM"),
     xlab = "PCs",
     ylab = "CRPS"
)

# highlight minimum and add minimum reference line 
abline(h = CRPS_min, lty = "dashed", col = adjustcolor("blue",alpha = .5))

dev.off()


############ Get Scatter plot for marginal correction factors for a certain number of principal components #######

m = 1


load(file = paste0(PCA_dir,"CovRes_mon",m,".RData"))

PCA = irlba::irlba(res_cov, nv = max(PCs))

land_ids <- which(DT[year == min(year) & month == min(month), is.na(Ens_bar) | is.na(SST_bar)])  

PCA_DT = DT[year == min(year) & month == min(month),][-land_ids,.(Lon,Lat,grid_id)]

for(d in  PCs){
  PCA_DT [,paste0("PC",d) := PCA$d[d]*PCA$u[,d]]
} 

# also get marginal variances
variances = list()
d = min(PCs) 
vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
variances[[d]] = vec^2

for(d in PCs[2:length(PCs)]){
  vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
  variances[[d]] = variances[[d-1]] + vec^2
}
names(variances) = paste0("var",PCs)
PCA_DT = data.table(PCA_DT,rbindlist(list(variances)))  


for(PC in c(1,5,10,15,20,25,30,40,50))
{
SD_PCA = sqrt(PCA_DT[,eval(parse(text = paste0("var",PC)))])
print(SD_PCA[1:5])
y = 2001

SD_hat = DT[year == y & month == m,][-land_ids,SD_hat]

plot(SD_hat/sqrt(PCA_DT[,eval(parse(text = paste0("var",PC)))]),
     ylim = c(0,2),
     ylab = SD_hat/SD_PCA,
     xlab = location)

}
