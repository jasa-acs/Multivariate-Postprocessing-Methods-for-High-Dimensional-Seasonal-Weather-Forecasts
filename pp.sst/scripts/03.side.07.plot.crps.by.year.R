# This script generates plots of CRPS for variance correction by moving averages.
# The CRPS are averaged over all months and all years in the validation period, and all available locations
# Requires the oos. version of script 03.master to be run previously.


##### setting up ######

rm(list = ls())

options(max.print = 1e3)

library(pp.sst)

name_abbr = "Full" 

save_dir = file.path('~','SST','Derived', name_abbr)

load(file = paste0(save_dir,"setup.RData"))

time_s37 = proc.time()


#### plotting ####

# get data for last year
y = max(validation_years)


row_sma = msc_sd_sma[year == y,-1,with = FALSE]
row_ema = msc_sd_ema[year == y,-1,with = FALSE]

y_range = range(list(row_sma[,-'min_l',with = FALSE],row_ema[,-'min_a',with = FALSE]))  

## plot for sma ##

pdf(paste0(plot_dir,"mean_scores_sd_sma.pdf"))
plot(x = win_length, 
     y = as.vector(row_sma[,-c('min_crps','min_l'),with = FALSE]),
     ylim = y_range,
     type = "b",
     col = "blue",
     main = paste0("CRPS for variance estimation by SMA"),
     xlab = "window length",
     ylab = "CRPS"
)

# highlight minimum and add minimum reference line 
abline(h = row_sma[,min_crps], lty = "dashed", col = adjustcolor("blue",alpha = .5))

points(x = row_sma[,min_l],
       y = row_sma[,min_crps],
       col = "blue",
       bg = "blue",
       pch = 21)
dev.off()

## plot for ema ##

pdf(paste0(plot_dir,"mean_scores_sd_ema.pdf"))
plot(x = par_vec, 
     y = row_ema[,-c('min_crps','min_a'),with = FALSE],
     ylim = y_range,
     type = "b",
     col = "blue",
     main = paste0("CRPS for variance estimation by EMA"),
     xlab = "scale parameter",
     ylab = "CRPS"
)

# highlight minimum and add minimum reference line 
abline(h = row_ema[,min_crps], lty = "dashed", col = adjustcolor("blue",alpha = .5))

points(x = row_ema[,min_a],
       y = row_ema[,min_crps],
       col = "blue",
       bg = "blue",
       pch = 21)
dev.off()

#### time, update script counter, save ####

time_s37 = proc.time() - time_s37

save.image(file = paste0(save_dir,"setup.RData"))

