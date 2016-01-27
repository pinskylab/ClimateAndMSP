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
wind_AKw <- readGDAL('./cmsp_data/wind_output/npv_US_millions_Alaska_west1050.tif') # read in SpatialGridDataFrame
	# image(wind_AKw)
wind_AKe <- readGDAL('./cmsp_data/wind_output/npv_US_millions_Alaska_east1050.tif') # read in SpatialGridDataFrame
wind_notAK <- readGDAL('./cmsp_data/wind_output/npv_US_millions_notAlaska_1050.tif')

load(paste(climgridfolder, 'climGrid_rcp85.proj2.RData', sep=''))
#load(paste(climgridfolder, 'climGrid.proj2_2015-02-10.RData', sep='')) # on macbook, temporary hack

# combine into a list
wind <- list(wind_AKw, wind_AKe, wind_notAK)

# project to LL. also converts to SpatialPointsDataFrame
wind.t <- wind
for(i in 1:length(wind)){
	wind.t[[i]] <- spTransform(wind[[i]], CRS('+proj=longlat +data=WGS84'))
#		colfun <- colorRamp(c('white', 'blue'))
#		sc <- (wind.t[[i]]$band1- min(wind.t[[i]]$band1))/(max(wind.t[[i]]$band1) - min(wind.t[[i]]$band1))
#		plot(wind.t[[i]], col=rgb(colfun(sc), maxColorValue=255), pch=16, cex=0.05) # plots
#		plot(wind.t[[i]], col=c('blue', 'red')[1+(wind.t[[i]]$band1> -600)], pch=16, cex=0.05) # plots of <> a threshold
}
	
# Add a grid indicator
gridsize=0.25 # size of grid of the climate data, in degrees
for(i in 1:length(wind.t)){
	wind.t[[i]]$latgrid = floor(coordinates(wind.t[[i]])[,2]/gridsize)*gridsize + gridsize/2 # round to nearest grid center
	wind.t[[i]]$longrid = floor(coordinates(wind.t[[i]])[,1]/gridsize)*gridsize + gridsize/2
}

# Summarize by grid cell
wind.sum <- wind.t
for(i in 1:length(wind.t)){
	wind.sum[[i]] <- aggregate(list(wind_npv = wind.t[[i]]$band1), by=list(latgrid = wind.t[[i]]$latgrid, longrid = wind.t[[i]]$longrid), FUN=mean)
}
	# plot to make sure it worked
	i <- 1 # pick the region
	par(mfrow=c(1,2))
	plot(wind.t[[i]], col=c('blue', 'red')[1+(wind.t[[i]]$band1> -600)], pch=16, cex=0.05) # plots of <> a threshold
	plot(wind.sum[[i]]$lon, wind.sum[[i]]$lat, col=c('blue', 'red')[1+(wind.sum[[i]]$wind_npv > -600)], pch=16, cex=0.2)	

# concatenate the regions together
wind.out <- wind.sum[[1]]
if(length(wind.sum)>1){
	for(i in 2:length(wind.sum)){
		wind.out <- rbind(wind.out, wind.sum[[i]])
	}
}
# write out
write.csv(wind.out, 'cmsp_data/wind_output/wind_npv.sum')