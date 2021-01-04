
# This script generates plots of RMSEs for bias correction by moving averages.
# The RMSEs are taken over all months and all years in the validation period, and all available locations
# Requires the oos. version of script 02.master


##### setting up ######

rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "test" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))

# start timer

time_s21 = proc.time()

#### plotting ####

# get data for last year
y = max(validation_years)


row_sma = msc_sma[year == y,-1,with = FALSE]
row_ema = msc_ema[year == y,-1,with = FALSE]

y_range = range(list(row_sma[,-'min_l',with = FALSE],row_ema[,-'min_a',with = FALSE]))  

y_range = sqrt(y_range)

## plot for sma ##

pdf(paste0(plot_dir,"mean_scores_sma.pdf"))
plot(x = win_length, 
     y = as.vector(sqrt(row_sma[,-c('min_MSE','min_l'),with = FALSE])),
     ylim = y_range,
     type = "b",
     col = "blue",
     main = paste0("RMSE for bias correction by SMA"),
     xlab = "window length",
     ylab = "RMSE"
)

# highlight minimum and add minimum reference line 
abline(h = row_sma[,sqrt(min_MSE)], lty = "dashed", col = adjustcolor("blue",alpha = .5))

points(x = row_sma[,min_l],
       y = row_sma[,sqrt(min_MSE)],
       col = "blue",
       bg = "blue",
       pch = 21)
dev.off()

## plot for ema ##

pdf(paste0(plot_dir,"mean_scores_ema.pdf"))
plot(x = par_vec, 
     y = sqrt(row_ema[,-c('min_MSE','min_a'),with = FALSE]),
     ylim = y_range,
     type = "b",
     col = "blue",
     main = paste0("RMSE for bias correction by EMA"),
     xlab = "scale parameter",
     ylab = "RMSE"
)

# highlight minimum and add minimum reference line 
abline(h = row_ema[,sqrt(min_MSE)], lty = "dashed", col = adjustcolor("blue",alpha = .5))

points(x = row_ema[,min_a],
       y = row_ema[,sqrt(min_MSE)],
       col = "blue",
       bg = "blue",
       pch = 21)
dev.off()

#### time, update script counter, save ####

time_s21 = proc.time() - time_s21

save.image(file = paste0(save_dir,"setup.RData"))

