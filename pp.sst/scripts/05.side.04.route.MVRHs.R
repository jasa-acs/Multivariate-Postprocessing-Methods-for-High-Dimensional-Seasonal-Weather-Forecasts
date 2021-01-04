###################################################################################################

#############  side script 5.4 - multivariate rank histograms for SST along route  ################

###################################################################################################

# See case study in the paper. Requires previous run of side script 05.side.03.route.scoring.R




#### setting up ######

rm(list = ls())

options(max.print = 1e3)

library(pp.sst)

name_abbr = "Full" 

save_dir = file.path('~','SST','Derived', name_abbr)

load(file = paste0(save_dir,"setup.RData"))


mc_cores = 5

##################

brks = 6

for(mod in mod_vec)
{
  print(mod)
  ens_size = fc_ens_size
  if(mod == 'ECC') ens_size = 9
  if(mod == 'Schaake') ens_size = length(training_years)

  temp = mv_rank_hist_new(get(paste0(mod,'_fc_route')), fc_ens_size = ens_size,
                          mc_cores = mc_cores,
                          breaks = brks, mn = paste0(mod," rank histograms"),
                          save_pdf = TRUE, plot_dir = plot_dir, file_name = paste0("rank_histo_",mod,"_route"))

  assign(paste0('rks_',mod,'_route'),temp)
}
  


### multivariate rank histograms for route ###

# range of probability the probability plots
rr = c(0,100)

titles = c(latex2exp::TeX('$\\widehat{\\Sigma}^{mc}$'),
           latex2exp::TeX('$\\widehat{\\Sigma}^{ac}$'),
           'GS',
           'ECC',
           'Schaake'
           )

### average rank histograms ###

pdf(paste0(plot_dir,'avg_rhs_route.pdf'),width = 35,height = 8)

  par(oma = c(1,1,1,1), mfrow=c(1,5), mar=c(2,1,2,1) )
  par('cex' = 3,'cex.axis' = 0.75)
  
  ind = 1
  for(mod in mod_vec)
  {
  
    ens_size = fc_ens_size
    if(mod == 'ECC') ens_size = 9
    if(mod == 'Schaake') ens_size = length(training_years)
    
    
    temp = get(paste0('rks_',mod,'_route'))
    
    # for the Schaake shuffle the number of bins does not divide the number of ranks:
    if(mod == 'Schaake')
    {
      temp = data.table(YM = rks_Schaake_route[,YM],av.rk.obs = hist_transform(rks_Schaake_route[,av.rk.obs]))
    }
    
    rhist_dt(temp[,.(YM,av.rk.obs)], 
             max_rk = ens_size +1,
             breaks = brks, 
             hist_xlab = "",
             hist_ylim = rr)
    title(main = titles[ind],line = -1)
    
    ind = ind + 1
  }
  
  par('cex' = 3)
  title(main = 'Average rank',outer = TRUE,line = -1)

dev.off()


### band depth rank histograms ###

pdf(paste0(plot_dir,'bd_rhs_route.pdf'),width = 35,height = 8)

  par(oma = c(1,1,1,1), mfrow=c(1,5), mar=c(2,1,2,1) )
  par('cex' = 3,'cex.axis' = 0.75)
  
  ind = 1
  for(mod in mod_vec)
  {
    
    ens_size = fc_ens_size
    if(mod == 'ECC') ens_size = 9
    if(mod == 'Schaake') ens_size = length(training_years)
    
    
    temp = get(paste0('rks_',mod,'_route'))
    
    # for the Schaake shuffle the number of bins does not divide the number of ranks:
    if(mod == 'Schaake')
    {
      temp = data.table(YM = rks_Schaake_route[,YM],bd.rk.obs = hist_transform(rks_Schaake_route[,bd.rk.obs]))
    }
    
    rhist_dt(temp[,.(YM,bd.rk.obs)], 
             max_rk = ens_size +1,
             breaks = brks, 
             hist_xlab = "",
             hist_ylim = rr)
    title(main = titles[ind],line = -1)
    
    ind = ind + 1
  }
  
  par('cex' = 3)
  title(main = 'Band depth rank',outer = TRUE,line = -1)

dev.off()


##########

save.image(file = paste0(save_dir,"setup.RData"))

