## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	projfolder <- '../CEmodels_proj' # holds model projections (outside Git)
	modfolder <- '../CEModels' # holds the models (outside Git)
	climgridfolder <- '../data/'
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	projfolder <- 'CEmodels_proj'
	modfolder <- 'CEmodels'
	climgridfolder <- 'data/'
	}
# could add code for Lauren's working directory here

#######################
## Process NatCap data ##
#######################
require(rgdal)

# read in data
wind_AKw <- readGDAL('cmsp_data/Wind/output/npv_US_millions_Alaska_west1050.tif') # read in SpatialGridDataFrame
	# image(wind_AKw)
wind_AKe <- readGDAL('./cmsp_data/Wind/output/npv_US_millions_Alaska_east1050.tif') # read in SpatialGridDataFrame
wind_WC <- readGDAL('./cmsp_data/Wind/output/npv_US_millions_WestCoast1050.tif')
wind_GM <- readGDAL('./cmsp_data/Wind/output/npv_US_millions_GoMex1050.tif')
wind_NE <- readGDAL('./cmsp_data/Wind/output/npv_US_millions_Northeast1050.tif')
	# image(wind_NE)

load(paste(climgridfolder, 'climGrid_rcp85.proj2.RData', sep='')) # loads clim. has the grid cells we want to summarize to

# combine wind output into a list
wind <- list(wind_AKw, wind_AKe, wind_WC, wind_GM, wind_NE)

# project to LL. also converts to SpatialPointsDataFrame
wind.t <- wind
for(i in 1:length(wind)){
	wind.t[[i]] <- spTransform(wind[[i]], CRS('+proj=longlat +data=WGS84')) # produces warnings about coercing to points. this is ok
#		colfun <- colorRamp(c('white', 'blue'))
#		sc <- (wind.t[[i]]$band1- min(wind.t[[i]]$band1))/(max(wind.t[[i]]$band1) - min(wind.t[[i]]$band1))
#		plot(wind.t[[i]], col=rgb(colfun(sc), maxColorValue=255), pch=16, cex=0.05) # plots
#		plot(wind.t[[i]], col=c('blue', 'red')[1+(wind.t[[i]]$band1> -600)], pch=16, cex=0.05) # plots of <> a threshold
}
	
# Add a grid indicator
gridsize=0.25 # size of grid of the climate data, in degrees
for(i in 1:length(wind.t)){
	wind.t[[i]]$lat <- floor(coordinates(wind.t[[i]])[,2]/gridsize)*gridsize + gridsize/2 # round to nearest grid center
	wind.t[[i]]$lon <- floor(coordinates(wind.t[[i]])[,1]/gridsize)*gridsize + gridsize/2
}

# Summarize by grid cell
wind.sum <- wind.t
for(i in 1:length(wind.t)){
	wind.sum[[i]] <- aggregate(list(wind_npv = wind.t[[i]]$band1), by=list(lat = wind.t[[i]]$lat, lon = wind.t[[i]]$lon), FUN=mean)
}
	# plot to make sure it worked
	i <- 5 # pick the region to plot
	par(mfrow=c(1,2))
	plot(wind.t[[i]], col=c('blue', 'red')[1+(wind.t[[i]]$band1> 0)], pch=16, cex=0.05) # plots of <> a threshold
	plot(wind.sum[[i]]$lon, wind.sum[[i]]$lat, col=c('blue', 'red')[1+(wind.sum[[i]]$wind_npv > 0)], pch=16, cex=0.2)	

# concatenate the regions together
wind.out <- wind.sum[[1]]
if(length(wind.sum)>1){
	for(i in 2:length(wind.sum)){
		wind.out <- rbind(wind.out, wind.sum[[i]])
	}
}

# convert to positive longitude, to match climgrid
wind.out$lon[wind.out$lon<0] <- wind.out$lon[wind.out$lon<0] + 360
	range(wind.out$lon)
	
	# whole plot
	plot(wind.out$lon, wind.out$lat, col=c('blue', 'red')[1+(wind.out$wind_npv > 0)], pch=16, cex=0.3)	

# mark which grids are in clim
wind.out$keep <- FALSE
wind.out$keep[paste(wind.out$lat, wind.out$lon) %in% paste(clim$lat, clim$lon)] <- TRUE

	# compared untrimmed and trimmed
	par(mfrow=c(1,2))
	plot(wind.out$lon, wind.out$lat, col=c('blue', 'red')[1+(wind.out$wind_npv > 0)], pch=16, cex=0.3)	
	plot(wind.out$lon[wind.out$keep], wind.out$lat[wind.out$keep], col=c('blue', 'red')[1+(wind.out$wind_npv[wind.out$keep] > 0)], pch=16, cex=0.3)	

# remove grids not in clim
wind.out <- wind.out[wind.out$keep, c('lat', 'lon', 'wind_npv')]

# write out
write.csv(wind.out, 'cmsp_data/wind_npv.csv')