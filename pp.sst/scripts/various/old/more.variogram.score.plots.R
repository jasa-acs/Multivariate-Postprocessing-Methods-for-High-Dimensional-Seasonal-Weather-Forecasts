

sc_1 = copy(scores_SE)
setnames(sc_1,"sc","sc_SE")
sc_1[,mean_sc := NULL]
setkeyv(sc_1,c("month","year"))

sc_2 = scores_mc[d == 30,]
setnames(sc_2,"sc","sc_PCA")
sc_2[,c("d","mean_sc"):=NULL]
setkeyv(sc_2,c("month","year"))

sc_3 = copy(scores_geostat_mc)
setnames(sc_3,"sc","sc_gs")
sc_3[,mean_sc := NULL]
setkeyv(sc_3,c("month","year"))

sc_4 = copy(sc_ECC)
setnames(sc_4,"sc","sc_ECC")
setkeyv(sc_4,c("month","year"))



all_scores = merge(sc_1,merge(sc_2,merge(sc_3,sc_4)))


rr = range(all_scores[,.SD,.SDcols = c("sc_SE","sc_PCA","sc_gs","sc_ECC")])

y=2001

for(y in 2001:2010)
{
 pdf(file = paste0(plot_dir,"vs_bm_y",y,".pdf"))
  plot(all_scores[year == y,.(month,sc_gs)], 
       col = "blue", type = "b",
       main = paste0("Variogram scores by month for ",y),
       xlab = "month", ylab = "score",
       ylim = rr )
  
  
  lines(all_scores[year == y,.(month,sc_PCA)], 
        col = "darkgreen", type = "b")
  

  lines(all_scores[year == y,.(month,sc_SE)], 
        col = "darkred", type = "b")

  lines(all_scores[year == y,.(month,sc_ECC)], 
        col = "orchid4", type = "b")
  
    
  abline(h = rowMeans(all_scores[,lapply(.SD,mean),.SDcols = c("sc_PCA","sc_gs","sc_SE")]),col = "grey", lty = "dashed")
  
  legend(x = "topleft",legend = c("geostat","PCA","EV shrinking","ECC","mean"), col = c("blue","darkgreen","darkred","orchid4","grey"),lty = c(1,1,1,1,2))
  dev.off()
}


#by year

mon_names = c("Jan","Feb","Mar","Apr","May","June","July","Aug","Sep","Oct","Nov","Dec")

for(m in 1:12)
{
 pdf(file = paste0(plot_dir,"vs_by_m",m,".pdf"))
  plot(all_scores[month == m,.(year,sc_gs)], 
       col = "blue", type = "b",
       main = paste0("Variogram scores by year for ",mon_names[m]),
       xlab = "year", ylab = "score",
       ylim = rr )
  
  
  lines(all_scores[month == m,.(year,sc_PCA)], 
        col = "darkgreen", type = "b")
  
  
  
  
  lines(all_scores[month == m,.(year,sc_SE)], 
        col = "darkred", type = "b")
  
  lines(all_scores[month == m,.(year,sc_ECC)], 
        col = "orchid4", type = "b")
  
  abline(h = rowMeans(all_scores[,lapply(.SD,mean),.SDcols = c("sc_PCA","sc_gs","sc_SE")]),col = "grey", lty = "dashed")
  

  legend(x = "topleft",legend = c("geostat","PCA","EV shrinking","ECC","mean"), col = c("blue","darkgreen","darkred","orchid4","grey"),lty = c(1,1,1,1,2))
  dev.off()

}


pdf(file = paste0(plot_dir,"boxplot.pdf"))
boxplot(all_scores[,.(sc_gs,sc_PCA,sc_SE,sc_ECC)],use.cols = TRUE, main  = "variogram scores for NAO")
dev.off()


