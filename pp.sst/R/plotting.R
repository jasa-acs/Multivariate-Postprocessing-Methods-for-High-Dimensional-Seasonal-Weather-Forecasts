


#' Quick diagnostic plotting function
#' 
#' @description Takes a data table that contains cols Lon, Lat and var. If it has multiple years and months, the minimum is taken
#' 
#' @param dt The data table.
#' @param var Character string. The name of the column containing the values for the plot. Default is third column, for subset data tables of the form .(Lon,Lat,var).
#' @param mn Title of the plot.
#' @param col_scheme Either "bwr" for blue - white - red, "wr" for white - red, or "wb" for white - blue. Specifies the color scheme of the plot. 
#' @param set_white Forces the blue-white-red color scheme to center white at the set value if specified.
#' @param xlab,ylab Labeling.
#' @param brks Vector containing breaks for the temperature scale. If NULL, 10 equally spaced values are picked
#' @param save_pdf,save_dir,file_name Whether, where and under which name the plot should be saved.
#' @param stretch_par Numeric. Only used when save_pdf == TRUE. Stretches the pdf output. Default is NULL, where it is stretched to #Lons/#Lats.
#'
#' @return none
#'  
#' @export
#' 
#' @author Claudio Heinrich
#' @examples \dontrun{
#' dt = load_combined_wide()
#' plot_diagnostic(dt)
#' }
#' 
#' @importFrom fields designer.colors image.plot
#' @importFrom maps map

plot_diagnostic = function( dt, var = colnames(dt)[3], mn = var,
                            rr = NULL, 
                            col_scheme = "bwr", set_white = NULL,
                            xlab = "", ylab = "",
                            brks = NULL,
                            save_pdf = FALSE, save_dir = "./figures/", file_name = "diag_plot",stretch_par = NULL)
{
 
  # prepare data table
  
  if("year" %in% colnames(dt))
  {
    if("month" %in% colnames(dt))
    {
      dt = dt[year == min(year) & month == min(month),.SD,.SDcols = c('Lon','Lat',var)][order(Lat,Lon)]
    } else {
      dt = dt[month == min(month),.SD,.SDcols = c('Lon','Lat',var)][order(Lat,Lon)]
    }
  } else {
    dt = dt[,.SD,.SDcols = c('Lon','Lat',var)][order(Lat,Lon)]
  }
  
  
  #--- get longitudes and latitudes
  
  Lons = sort(unique(dt[,Lon]))
  Lats = sort(unique(dt[,Lat]))
  
  n_lon = length(Lons)
  n_lat = length(Lats)
  
  
  if(is.null(rr))  rr = range(dt[,3],na.rm=TRUE)
  if(!is.null(rr)){
    dt[dt[[3]]<min(rr),3] = min(rr)
    dt[dt[[3]]> max(rr),3] = max(rr)
  }
  
  A = matrix(dt[[3]],  n_lon, n_lat)
  
  # --- scaling and colors ---
  
  brk = seq(rr[1],rr[2],length = 500)
  if(is.null(brks))
  {
    brk.ind = round(seq(1,length(brk),length = 10))
    brk.lab = round(brk[brk.ind],2)
    brk.at = brk[brk.ind] 
  } else {
    brk.lab = brks
    brk.at = brks
  }
  
  if(col_scheme == "bwr"){
    if(is.null(set_white)){
    color <- fields::designer.colors(n=length(brk)-1, col = c("darkblue","white","darkred"))
    }else{
       zero.ind = min(which(brk > set_white))/length(brk)
       color <- fields::designer.colors(n=length(brk)-1, col = c("darkblue","white","darkred"), x = c(0,zero.ind,1))
    }
  }
  if(col_scheme == "wr"){
    color <- fields::designer.colors(n=length(brk)-1, col = c("white","darkred"))
  }
  if(col_scheme == "wb"){
    color <- fields::designer.colors(n=length(brk)-1, col = c("white","blue"))
  }
    
  # set up for coloring NAs gray
  
  newz.na <- rr[2]+(rr[2]-rr[1])/length(color) # new z for NA
  A[which(is.na(A))] <- newz.na 
  rr[2] <- newz.na      #extend the range to include the new value 
  color <- c(color, 'gray') # extend the color palette by 'gray'
  brk = c(brk,rr[2]) # extend the vector of breaks
  
  #--- plotting ---
  
  if(save_pdf) 
  {
    if (is.null(stretch_par)) stretch_par = n_lat/n_lon
    
    par_0 = par() # allow to set par manually before calling the function
    
    pdf(paste0(save_dir,file_name,".pdf"),width = 7,height = stretch_par * 7)
    
    # somehow this doesn't pass on font sizes?
    par('cex' = par_0$cex)
    par('cex.lab' = par_0$cex.lab)
    par('cex.axis' = par_0$cex.axis)
    par('mfrow' = par_0$mfrow)
  }
  
  par(mar = c(2,2,2,2))
  
  fields::image.plot(Lons,Lats,A,
             zlim=rr, main = mn,
             xlim = range(Lons), xlab=xlab,
             ylim = range(Lats), ylab=ylab,
             breaks=brk,
             col=color,
             axis.args=list(at = brk.at,
                            label = brk.lab))
   # add world map
  maps::map("world", add = TRUE)
  
  if(save_pdf) dev.off()
  
}



#' smooth plotting function
#' 
#' @description Takes a data table of the form .(Lon,Lat,value) or .(Lat,Lon,value) and plots values on the globe after applying a kernel smoothing
#' 
#' @param dt The data table.
#' @param var Character string. The name of the column containing the values for the plot. Default is the name of the third column of dt, to be used for subset data tables of the form .(Lon,Lat,var).
#' @param mn,rr Title and range of the plot.
#' @param theta parameter for the Gaussian smoothing kernel
#' @param pixels Resolution of the plot
#' @param col_scheme Either "bwr" for blue - white - red, "wr" for white - red, or "wb" for white - blue. Specifies the color scheme of the plot. 
#' @param set_white Forces the blue-white-red color scheme to center white at the set value if specified.
#' @param xlab,ylab Labeling.
#' @param brks vector of breaks for the temperature scale, if NULL, ten equally spaced numbers are picked.
#' @param save_pdf,save_dir,file_name Whether, where and under which name the plot should be saved.
#' @param stretch_par Numeric. Only used when save_pdf == TRUE. Stretches the pdf output. Default is NULL, where it is stretched to #Lons/#Lats.
#'
#' @return none
#'  
#' @export
#' 
#' @author Claudio Heinrich
#' @examples \dontrun{
#' dt = load_combined_wide()
#' plot_smooth(dt)
#' }
#' 
#' @importFrom fields as.image designer.colors image.plot image.smooth
#' @importFrom maps map
#' @import sp maptools


plot_smooth = function( dt, var = colnames(dt)[3], mn = var, rr = NULL,...,
                        theta = 0.5, pixels = 256,
                        col_scheme = "bwr", set_white = NULL,
                        xlab = "", ylab = "",
                        brks = NULL,
                        save_pdf = FALSE, save_dir = "./figures/", file_name = "diag_plot", stretch_par = NULL)
{
  # prepare data table
  
  if("year" %in% colnames(dt))
  {
    if("month" %in% colnames(dt))
    {
      dt = dt[year == min(year) & month == min(month),.SD,.SDcols = c('Lon','Lat',var)][order(Lat,Lon)]
    } else {
      dt = dt[month == min(month),.SD,.SDcols = c('Lon','Lat',var)][order(Lat,Lon)]
    }
  } else {
    dt = dt[,.SD,.SDcols = c('Lon','Lat',var)][order(Lat,Lon)]
  }
  
  
  #--- create image ---
  
  x = dt[,.(Lon,Lat)]
  setnames(x,c("Lon","Lat"), c("lon","lat"))
  
  Lons = unique(dt[,Lon])
  Lats = unique(dt[,Lat])
  
  n_lon = length(Lons)
  n_lat = length(Lats)
  
  A = matrix(dt[[3]],  n_lon, n_lat)
  
  im_0 = fields::image.smooth(fields::as.image(A,x = x,nx = pixels,ny = pixels),theta = theta)
  
  ## Find the points that fall over land
  
  if(!exists("wrld_simpl")) data(wrld_simpl, package = 'maptools') 
  
  all_loc = expand.grid(lat = im_0$x,lon = im_0$y)
  pts <- sp::SpatialPoints(all_loc, proj4string=sp::CRS(proj4string(wrld_simpl)))
  ii <- !is.na(over(pts, wrld_simpl)$FIPS)
  im_0$z[ii] = NA
  
  # --- fix range of plot and fill in values for points out of range ---
  
  if(is.null(rr))  rr = range(im_0$z,na.rm=TRUE)
  if(!is.null(rr)){
    im_0$z[im_0$z< min(rr)] = min(rr)
    im_0$z[im_0$z> max(rr)] = max(rr)
  }
  
  # --- scaling and colors ---
  
  brk = seq(rr[1],rr[2],length = 500)
  if(is.null(brks))
  {
    brk.ind = round(seq(1,length(brk),length = 10))
    brk.lab = round(brk[brk.ind],2)
    brk.at = brk[brk.ind] 
  } else {
    brk.lab = brks
    brk.at = brks
  }
  
  if(col_scheme == "bwr"){
    if(is.null(set_white)){
      color <- fields::designer.colors(n=length(brk)-1, col = c("darkblue","white","darkred"))
    }else{
      zero.ind = min(which(brk > set_white))/length(brk)
      color <- fields::designer.colors(n=length(brk)-1, col = c("darkblue","white","darkred"), x = c(0,zero.ind,1))
    }
  }
  if(col_scheme == "wr"){
    color <- fields::designer.colors(n=length(brk)-1, col = c("white","darkred"))
  }
  if(col_scheme == "wb"){
    color <- fields::designer.colors(n=length(brk)-1, col = c("white","blue"))
  }
  
  # color NAs grey
  
  newz.na <- rr[2]+(rr[2]-rr[1])/length(color) # new z for NA
  im_0$z[which(is.na(im_0$z))] <- newz.na 
  rr[2] <- newz.na # extend the range to include the new value 
  color <- c(color, 'gray') # extend the color range by gray
  brk = c(brk,rr[2]) # extend the vector of breaks
  
  #--- plotting ---
  
  if(save_pdf) 
    {
    if (is.null(stretch_par)) stretch_par = n_lat/n_lon
    
    par_0 = par() # allow to set par manually before calling the function
    
    pdf(paste0(save_dir,file_name,".pdf"),width = 7,height = stretch_par * 7)
    
    suppressWarnings(par(par_0))
    # somehow this doesn't pass on font sizes?
    par('cex' = par_0$cex)
    par('cex.lab' = par_0$cex.lab)
    par('cex.axis' = par_0$cex.axis)
    par('mfrow' = par_0$mfrow)
    }
  
  par(mar = c(2,2,2,2))
  
  fields::image.plot(im_0,
                     zlim=rr, main = mn,...,
                     xlim = range(Lons), xlab=xlab,
                     ylim = range(Lats), ylab=ylab,
                     breaks=brk,
                     col=color,
                     axis.args=list(at = brk.at,
                                    label = brk.lab))
  
  # add world map
  
  maps::map("world", add = TRUE)
  
  if(save_pdf) dev.off()
  
}


