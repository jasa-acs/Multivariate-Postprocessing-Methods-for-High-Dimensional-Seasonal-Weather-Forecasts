#############################################################################

#############  side script 3.5 - univariate rank histograms  ################

#############################################################################

# This script generates univariate rank histogram plots (not shown in the paper).

##### setting up ######

rm(list = ls())

options(max.print = 1e3)

library(pp.sst)

name_abbr = "Full" 

save_dir = file.path('~','SST','Derived', name_abbr)

load(file = paste0(save_dir,"setup.RData"))


###########

na_loc = which(DT[,is.na(SST_bar) | is.na(Ens_bar) ])  # land locations

# number of simulations and of breaks in the histogram: choose n = 10*k-1 and nbreaks = 10*k
n=99
nbreaks = 21


means = DT[-na_loc,][year %in% validation_years & month %in% months,SST_hat]
sds = DT[-na_loc,][year %in% validation_years & month %in% months,SD_hat]
obs = DT[-na_loc,][year %in% validation_years & month %in% months,SST_bar]

fcs = trc(rnorm(n*length(means),mean = means,sd = sds))

rh_mat = matrix(c(obs,fcs),ncol = n+1)
  
rank_mat = t(apply(rh_mat, MARGIN = 1, rank, ties = 'random'))
  
pdf(file=paste0(plot_dir,"rkh.pdf"))
  hist(rank_mat[,1],
       breaks=seq(0, n+1, length.out = nbreaks),
       main = 'univariate rank histogram',
       xlab = "", ylab = "", axes=FALSE, col="gray80", border="gray60")
  abline(a=dim(rank_mat)[1]/(nbreaks-1), b=0, lty=2, col="gray30")
dev.off()

### time and save ###

time_s35 = proc.time() - time_s35

save.image(file = paste0(save_dir,"setup.RData"))