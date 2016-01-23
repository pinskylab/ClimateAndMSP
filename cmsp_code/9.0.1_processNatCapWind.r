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
wind_ne <- readGDAL('./cmsp_data/wind_output/npv_US_millions_Northeast.tif') # read in SpatialGridDataFrame
	# image(wind_ne)
#wind_wc <- readGDAL('./cmsp_data/wind_output/npv_US_millions_WestCoast.tif')
#wind_gx <- readGDAL('./cmsp_data/wind_output/npv_US_millions_GoMex.tif')
#wind_ak <- readGDAL('./cmsp_data/wind_output/npv_US_millions_Alaska.tif')

load(paste(climgridfolder, 'climGrid_rcp85.proj2.RData', sep=''))

# project to LL. also converts to SpatialPointsDataFrame
wind_ne.t <- spTransform(wind_ne, CRS('+proj=longlat +data=WGS84'))
	colfun <- colorRamp(c('white', 'blue'))
	sc <- (wind_ne.t$band1- min(wind_ne.t$band1))/(max(wind_ne.t$band1) - min(wind_ne.t$band1))
	plot(wind_ne.t, col=rgb(colfun(sc), maxColorValue=255), pch=16, cex=0.05) # plots
	plot(wind_ne.t, col=c('blue', 'red')[1+(wind_ne.t$band1> -600)], pch=16, cex=0.05) # plots of <> a threshold
	
# Add a grid indicator
gridsize=0.25 # size of grid of the climate data, in degrees
wind_ne.t$latgrid = floor(coordinates(wind_ne.t)[,2]/gridsize)*gridsize + gridsize/2 # round to nearest grid center
wind_ne.t$longrid = floor(coordinates(wind_ne.t)[,1]/gridsize)*gridsize + gridsize/2

# Summarize by grid cell
wind_ne.sum <- aggregate(list(wind_npv = wind_ne.t$band1), by=list(latgrid = wind_ne.t$latgrid, longrid = wind_ne.t$longrid), FUN=mean)

	# plot to make sure it worked
	par(mfrow=c(1,2))
	plot(wind_ne.t, col=c('blue', 'red')[1+(wind_ne.t$band1> -600)], pch=16, cex=0.05) # plots of <> a threshold
	plot(wind_ne.sum$lon, wind_ne.sum$lat, col=c('blue', 'red')[1+(wind_ne.sum$wind_npv > -600)], pch=16, cex=0.2)	
