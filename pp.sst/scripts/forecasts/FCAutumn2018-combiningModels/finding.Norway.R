

# get the grid_ids in Norway and separate into three different regions:

rm(list = ls())

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/")

setwd(save_dir)

load(file = paste0("dt_combine_mr_wide.RData"))

DT = dt[year == 1985 & month == 1 & Lon < 32 & Lon > 0 & Lat <75 & Lat >55,]

land_ids = DT[is.na(Ens_bar),grid_id]

library(maps)

DT[grid_id %in% land_ids, country:=  map.where(database = "world", DT[grid_id %in% land_ids,Lon],DT[grid_id %in% land_ids,Lat])]

DT[,ind_nor := as.numeric(country == "Norway")]


# gives a decent approximation of Norway but needs manual correction:

cor_loc = matrix(c(rep(70.5,3),22.5,23.5,25.5),ncol = 2)
cor_loc = rbind(cor_loc,  c(69.5,18.5))
cor_loc = rbind(cor_loc,  matrix(c(rep(68.5,3),15.5,16.5,18.5),ncol = 2))
cor_loc = rbind(cor_loc,  matrix(c(rep(67.5,3),14.5,15.5,16.5),ncol = 2))
cor_loc = rbind(cor_loc,  c(66.5,15.5))
cor_loc = rbind(cor_loc,  c(65.5,12.5))
cor_loc = rbind(cor_loc,  c(64.5,10.5))
cor_loc = rbind(cor_loc,  c(63.5,10.5))
cor_loc = rbind(cor_loc,  c(62.5,6.5))
cor_loc = rbind(cor_loc,  c(59.5,5.5))
cor_loc = rbind(cor_loc,  c(59.5,10.5))
#cor_loc = rbind(cor_loc,  c(58.5,5.5))

cor_loc = as.data.table( cor_loc)
setnames(cor_loc,c("Lat","Lon"))
cor_loc[,"cor":=TRUE]

DT = merge(DT,cor_loc, by = c("Lon","Lat"),all.x = TRUE)
DT[cor == TRUE,ind_nor :=1]


DT[ind_nor == 0, ind_nor := NA]

#separate into three different regions: Norway west, east and north:

DT[Lat < 63.5 & Lon <= 7.5 & ind_nor == 1, ind_nor := 2]
DT[Lat >= 63.5 & ind_nor == 1, ind_nor := 3]

plot_diagnostic(DT[,.(Lon,Lat,ind_nor)])

# manually correct a few gridpoints:

DT[Lat == 62.5 & Lon == 8.5 , ind_nor := 3]
DT[Lat == 62.5 & Lon == 9.5 , ind_nor := 3]
DT[Lat == 58.5 & Lon == 7.5 , ind_nor := 1]
DT[Lat == 60.5 & Lon == 8.5 , ind_nor := 2]
DT[Lat == 61.5 & Lon == 8.5 , ind_nor := 2]

Norway_regions = DT[.(Lon,Lat,grid_id,ind_nor)]

plot_diagnostic(Norway_regions[,.(Lon,Lat,ind_nor)])

save(Norway_regions,file = paste0(save_dir,'Norway_regions.RData'))


