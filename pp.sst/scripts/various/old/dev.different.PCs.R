##### setting up ######

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "Full" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))


ex_plot_dir = paste0(plot_dir,"Examples/")
dir.create(ex_plot_dir, showWarnings = FALSE)


### example plots ###

ex_month = 9

cex = 1.1

# get climatology

training_years = 1985:2000

ex_years = 2001


# SST residuals

r_plus = max(abs(range(DT[year %in% ex_years & month %in% ex_month,SST_bar - SST_hat],na.rm = TRUE)))

rr = c(-r_plus,r_plus)

for(m in ex_month)
{
  for(y in ex_years)
  {
    mn = paste0("Forecast residual, ",m," / ",y)

    for(k in 1:repetitions)
    {
      plot_diagnostic(DT[year == y & month == m,.(Lon,Lat,SST_bar - SST_hat)],
                      mn = mn,
                      ylab = "", xlab = "",
                      rr = rr,
                      cex = cex,
                      save_pdf = TRUE,
                      save_dir = ex_plot_dir,
                      file_name = paste0("res",k,"_m",m,"_y",y))  
    }
  }
}


# get principal components

PCs = 10

DT[,SD_hat := SD_hat_sv]

##### setting up ######

PCA_dir = paste0(save_dir,"PCA/")
dir.create(PCA_dir, showWarnings = FALSE)

versions = c("aggr_by_season","sum_of_squares","wrt_ens_mean")

for(ver in versions)
{
  
  for_res_cov(Y = training_years,
              dt = DT, 
              save_dir = PCA_dir,
              ens_size = ens_size,
              version = ver)
  
  Sys.sleep(10)
  
  prin_comp_dt = get_PCs(dt = DT, y = ex_years, m = ex_month, PCA_depth = PCs, cov_dir = PCA_dir)
  
  
  
  #without marginal correction
  
  for(y in ex_years){
    rr_PC_raw = range(prin_comp_dt[year == y,.SD,.SDcols = c(paste0("PC",1:PCs))],na.rm = TRUE)
    rr_PC_raw = c(-max(abs(rr_PC_raw)), max(abs(rr_PC_raw)))
    
    for(d in 1:PCs)
    {
      plot_diagnostic(prin_comp_dt[year == y,.SD,.SDcols = c("Lon","Lat",paste0("PC",d))],
                      rr = rr_PC_raw,
                      mn = paste0("PC ",d,", nmc, ",ver,", 0",ex_month," / ",y),
                      xlab = "", ylab = "",
                      cex = cex,
                      save_pdf = TRUE,
                      save_dir = ex_plot_dir,
                      file_name = paste0("PC",d,"_y",y,"_raw_",ver))
      
    }
  }
  
  #with marginal correction:
  
  for(y in ex_years){
    rr_PC_raw = range(prin_comp_dt[year == y,.SD,.SDcols = c(paste0("PC_marcor_",1:PCs))],na.rm = TRUE)
    rr_PC_raw = c(-max(abs(rr_PC_raw)), max(abs(rr_PC_raw)))
    
    for(d in 1:PCs)
    {
      plot_diagnostic(prin_comp_dt[year == y,.SD,.SDcols = c("Lon","Lat",paste0("PC_marcor_",d))],
                      rr = rr_PC_raw,
                      mn = paste0("PC ",d,", ",ver,", 0",ex_month," / ",y),
                      save_pdf = TRUE,
                      xlab = "", ylab = "",
                      cex = cex,
                      save_dir = ex_plot_dir,
                      file_name = paste0("PC",d,"_y",y,"_mc_",ver))
      
    }
  }
  
  
}


### how are they so similar??? ###

PCs = 14

ver = "sum_of_squares"

for_res_cov(Y = training_years,
            dt = DT, 
            save_dir = PCA_dir,
            ens_size = ens_size,
            version = ver)

prin_comp_ssq_dt = get_PCs(dt = DT, y = ex_years, m = ex_month, PCA_depth = PCs, cov_dir = PCA_dir)


ver = "wrt_ens_mean"

for_res_cov(Y = training_years,
            dt = DT, 
            save_dir = PCA_dir,
            ens_size = ens_size,
            version = ver)

prin_comp_sc_dt = get_PCs(dt = DT, y = ex_years, m = ex_month, PCA_depth = PCs, cov_dir = PCA_dir)

# reorienting

cols = c(paste0("PC",1:PCs),paste0("PC_marcor_",1:PCs))
for(i in 1:length(cols))
{
  cur_col = cols[i]
  print(cur_col)
  dist_unflipped = sqrt(sum((prin_comp_sc_dt[,eval(parse(text = cur_col))] - prin_comp_ssq_dt[,eval(parse(text = cur_col))])^2,na.rm = TRUE))
  dist_flipped = sqrt(sum((prin_comp_sc_dt[,eval(parse(text = cur_col))]   + prin_comp_ssq_dt[,eval(parse(text = cur_col))])^2,na.rm = TRUE))
  if(dist_flipped < dist_unflipped)
  {
    vec = prin_comp_sc_dt[,eval(parse(text = cur_col))]
    prin_comp_sc_dt[, (cur_col) := -vec]
  }
}



err_nmc = c()
for(i in 1:PCs){
  err_nmc = c(err_nmc,sqrt(sum((prin_comp_sc_dt[,eval(parse(text = paste0("PC",i)))]-prin_comp_ssq_dt[,eval(parse(text = paste0("PC",i)))])^2,na.rm = TRUE)))
  }
err_nmc

err = c()
for(i in 1:PCs){
  err = c(err,sqrt(sum((prin_comp_sc_dt[,eval(parse(text = paste0("PC_marcor_",i)))]-prin_comp_ssq_dt[,eval(parse(text = paste0("PC_marcor_",i)))])^2,na.rm = TRUE)))
}
err

comparison_dt = prin_comp_ssq_dt[,.SD,.SDcols = c(paste0("PC",1:PCs),paste0("PC_marcor_",1:PCs))] - prin_comp_sc_dt[,.SD,.SDcols = c(paste0("PC",1:PCs),paste0("PC_marcor_",1:PCs))]

comparison_dt = data.table(prin_comp_ssq_dt[,1:6],comparison_dt)

rr = c(-max(abs(range(comparison_dt[,.SD,.SDcols = paste0("PC",1:PCs)], na.rm = TRUE))),max(abs(range(comparison_dt[,.SD,.SDcols = paste0("PC",1:PCs)], na.rm = TRUE))))
for(i in 1:PCs)
{
  plot_diagnostic(comparison_dt[year == 2002 & month == 9,.(Lon,Lat,eval(parse(text = paste0("PC",i))))],
                  rr = rr,
                  set_white = 0,
                  mn = paste0("difference of PC ",i),
                  save_pdf = TRUE,
                  save_dir = ex_plot_dir,
                  file_name = paste0("diff_PC",i))  
}

