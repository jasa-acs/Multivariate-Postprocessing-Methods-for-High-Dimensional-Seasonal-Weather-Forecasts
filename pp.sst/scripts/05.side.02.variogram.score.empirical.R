
###################################################################

#############  side script 5.2 - variogram scores  ################

###################################################################

# This script compares the multivariate post-processing methods based on their variogram score.
# 

#### setting up ######

rm(list = ls())

options(max.print = 1e3)

library(pp.sst)

name_abbr = "Full" 

save_dir = file.path('~','SST','Derived', name_abbr)

load(file = paste0(save_dir,"setup.RData"))

time_s52 = proc.time()

##################

# choose power for variogram score:

pp = 0.5

#######################
######### PCA #########
#######################

###### mc #####

load(paste0(PCA_dir,"fc_mc.RData"))


vs_PCA_mc = var_sc_par(dt_fc = PCA_fc_mc,  p = pp,
                    years = validation_years, ms = months, n = fc_ens_size,
                    save_dir = NULL,
                    mc_cores = mc_cores)


save.image(file = paste0(save_dir,"setup.RData"))

rm(PCA_fc_mc)
gc()

#####  ac ######

load(paste0(PCA_dir,"fc_ac.RData"))


vs_PCA_ac = var_sc_par(dt_fc = PCA_fc_ac,  p = pp,
                       years = validation_years, ms = months, n = fc_ens_size,
                       save_dir = NULL,
                       mc_cores = mc_cores)




save.image(file = paste0(save_dir,"setup.RData"))

rm(PCA_fc_ac)
gc()

###### GS ######

load(paste0(GS_dir,"fc.RData"))

vs_GS = var_sc_par(dt_fc = GS_fc,  p = pp,
                   years = validation_years, ms = months, n = fc_ens_size,
                   save_dir = NULL,
                   mc_cores = mc_cores)

rm(GS_fc)

save.image(file = paste0(save_dir,"setup.RData"))

###### ECC ######

load(paste0(ECC_dir,"fc.RData"))

vs_ECC = var_sc_par(dt_fc = ECC_fc,  p = pp,
                    years = validation_years, ms = months, n = ens_size, 
                    save_dir = NULL,
                    mc_cores = mc_cores)
  
rm(ECC_fc)

save.image(file = paste0(save_dir,"setup.RData"))


###### Schaake ######

load(paste0(Schaake_dir,"fc.RData"))

vs_Schaake = var_sc_par(dt_fc = Schaake_fc,  p = pp,
                    years = validation_years, ms = months, n = length(training_years), 
                    save_dir = NULL,
                    mc_cores = mc_cores)

rm(Schaake_fc)

###########

# combine, boxplot and save

vs_dt = data.table( PCA_mc = vs_PCA_mc[,mean(vs)],
                    PCA_ac = vs_PCA_ac[,mean(vs)],
                    GS = vs_GS[,mean(vs)], 
                    ECC = vs_ECC[,mean(vs)],
                    Schaake = vs_Schaake[,mean(vs)])

pdf(paste0(plot_dir,"boxplot_vs.pdf"))
boxplot( list (vs_PCA_mc[,vs],vs_PCA_ac[,vs],vs_GS[,vs], vs_ECC[,vs],vs_Schaake[,vs]),names = c("PCA mc","PCA ac","GS","ECC",'Schaake'), main = "variogram scores")
dev.off()

save(vs_dt,vs_PCA_ac,vs_PCA_mc,vs_GS,vs_ECC,vs_Schaake, file = paste0(save_dir,"vs.RData"))



### permutation tests ###

mod_vec = c('PCA_mc','PCA_ac','GS','ECC','Schaake')
n_mod = length(mod_vec)

N = 10000

p_vals_pt_dt = data.table(expand.grid(model1 = mod_vec,model2 = mod_vec))

pdf(paste0(plot_dir,'perm_test_vs.pdf'),width = 14,height = 14)

par( oma=c(0,6,6,0), mfrow=c(n_mod,n_mod), mar=c(3,3,2,2)+0.1 )

for(mod1 in mod_vec)
{
  
  for(mod2 in mod_vec)
  {
    vs_1 = get(paste0('vs_',mod1))
    vs_2 = get(paste0('vs_',mod2))
    
    perm_test_dt = merge(vs_1,vs_2,by=c('year','month'))
    setnames(perm_test_dt,c('vs.x','vs.y'), c(mod1,mod2))
    
    # permutation test for mod1 ~ mod2
    pt_vs = permutation_test_difference(perm_test_dt[,get(mod1)],perm_test_dt[,get(mod2)], N=N)
    
    p_vals_pt_dt[model1 == mod1 & model2 == mod2, p_val := pt_vs$p_val]  
    
    x_lim_max = 1.1*max(abs(c(pt_vs$d_bar,pt_vs$D)))
    
    hist(pt_vs$D, xlim = c(-x_lim_max,x_lim_max),main = paste0(mod1,' vs. ',mod2),breaks = 20)
    
    qq = quantile(pt_vs$D,c(0.025,0.975))
    
    abline(v = pt_vs$d_bar,col = 'red')
    
    abline(v = qq,lty = 2)
    
    qq = quantile(pt_vs$D,c(0.1,0.9))
    
  }
}

dev.off()


### time and save ###

time_s52 = proc.time() - time_s52

save.image(file = paste0(save_dir,"setup.RData"))

