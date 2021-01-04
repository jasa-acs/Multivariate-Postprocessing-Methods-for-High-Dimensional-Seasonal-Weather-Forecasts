
# plot forecasts for 2018 09-12, NorCPM, GCFS1, and mixture


##### setting up ######

rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)


name_abbr = "Aut_2018_small/mm" 
save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

plot_dir = "./figures/Aut_2018_small/"

load(file = paste0(save_dir,"setup.RData"))


#transform to climatology
clim_years = 1985:2014

DT_clim = dt_transform_center(DT,clim_years,
                              trafo_colnames = c("SST_bar","SST_hat","SST_NCPM","SST_GCFS1"))


rr = c(-5,5)
y=2018
cex = 1.5

for(m in 8:10)
{
plot_diagnostic(DT_clim[year == y & month == m, .(Lon,Lat,SST_NCPM)],
                mn = paste0("NorCPM prediction, ",m,"/",y),xlab = "",ylab = "",cex = cex,
                rr = rr,
                save_pdf = TRUE, save_dir = plot_dir,file_name = paste0("sst_fc_NCPM_m",m,"_y",y)
                )
plot_diagnostic(DT_clim[year == y & month == m, .(Lon,Lat,SST_GCFS1)],
                mn = paste0("GCFS1 prediction, ",m,"/",y),xlab = "",ylab = "",cex = cex,
                rr = rr,
                save_pdf = TRUE, save_dir = plot_dir,file_name = paste0("sst_fc_GCFS1_m",m,"_y",y)
                )
plot_diagnostic(DT_clim[year == y & month == m, .(Lon,Lat,SST_hat)],
                mn = paste0("mixed prediction, ",m,"/",y),xlab = "",ylab = "",cex = cex,
                rr = rr,
                save_pdf = TRUE, save_dir = plot_dir,file_name = paste0("sst_fc_mixed_m",m,"_y",y)
)

}

