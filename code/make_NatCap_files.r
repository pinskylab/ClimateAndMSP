# Make a shapefile of our Areas of Interest, for use in the NatCap's InVest software
# Also make a csv of lat/lon for landing points and grid connection points


#######################
## Load libraries
#######################
#require(maptools)
require(RColorBrewer)
#require(rgdal)
#require(PBSmapping)
#require(rgeos) # for gBuffer
require(sf)

##############################################################
## Read in and process data to create Areas of Interest (AOI)
##############################################################
#climgrid <- read.csv('data/climGrid.csv', row.names=1)
climgrid <- readRDS('temp/SPsf2.rds') # the projection grid
#land <- readOGR(dsn='cmsp_data/', layer='global_polygon') # global land layer from NatCap
#	plot(land) # takes too long

# label by area of interest regions
climgrid$AOI <- NA
climgrid$AOI[climgrid$lon < -100 | climgrid$lon > 0] <- 'west'
climgrid$AOI[climgrid$lon >= -100 & climgrid$lon < 0] <- 'east'
sum(is.na(climgrid$AOI)) # 0
#    plot(climgrid['AOI'], lwd=0.001, axes=TRUE) # Aleutians are on the far right
#    plot(climgrid['AOI'], lwd=0.5, xlim=c(-70, -69), ylim=c(41,42), axes=TRUE) # zoom in (Cape Cod)

# project to North American Albers
climgridp <- st_transform(climgrid, crs='+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs') # North America Albers Equal Area Conic from https://epsg.io/102008
#    plot(climgridp['AOI'], lwd=0.001, axes=TRUE) # Aleutians are on the far right

# buffer slightly (1000 m) to aid merging the grid cells together
climgridbuf <- st_buffer(climgridp, 1000)

# merge all into multipart polygons based on AOI
climpoly <- aggregate(climgridbuf[,c('geometry')], by = list(AOI = climgridbuf$AOI), do_union = TRUE, FUN=function(x, ...) return(1))	
#	 plot(climpoly['AOI'], lwd=0.01, axes=TRUE)
#    plot(climpoly['AOI'], lwd=1, xlim=c(1.9e6, 2.2e6), ylim=c(4e5,6e5), axes=TRUE) # zoom in (Cape Cod)
	
# buffer whole analysis area
climpolyb <- st_buffer(climpoly, dist = 100000) # 100km
#	plot(climpolyb, lwd=0.2, axes=TRUE)
	

# split Alaska along 179 and 181deg longitude (3 pieces)
# maskPS <- as.PolySet(data.frame(PID=rep(c(1:3),rep(4,3)), POS=rep(1:4, 3), X=rep(c(150,179.99,180.01,230), c(2,4,4,2)), Y=rep(c(45,65,45,65,45,65,45), c(1,2,2,2,2,2,1))), projection='LL')
# mask <- PolySet2SpatialPolygons(maskPS)
# mask.t <- spTransform(mask, proj4string(gridSPD))
# 	#plot(mask.t, col=1:length(mask.t))
# 	plot(gridSP.b, lwd=0.2, axes=TRUE, col=1:4)
# 	plot(mask.t, add=TRUE)
# gridSPAK <- gIntersection(gridSPD[grep('Alaska', gridSPD$region),], mask.t, byid=TRUE)
# 	plot(gridSPAK, col=1:length(gridSPAK))
# 
# 	# convert Alaska to spatialpolygonsdataframe for writing
# 	pid <- sapply(slot(gridSPAK, "polygons"), function(x) slot(x, "ID"))
# 	regs <- c('Alaska_west', 'Alaska_middle', 'Alaska_east')
# 	gridSPDAK <- SpatialPolygonsDataFrame(gridSPAK, data.frame(region=regs, row.names=pid))

# write out each region as a separate shapefile
# will throw an error if file already exists
for(i in 1:nrow(climpolyb)){
	st_write(climpolyb[i,], paste0('temp/AOI_', climpolyb$AOI[i], '.shp')) # write the .prj file as part of the shapefile. Overwrite
}

# old writing out, usefull if need to split AK in pieces
# notAK <- grep('Alaska', gridSPD$region, invert=TRUE)
# for(i in notAK){ # for all except Alaska
# 	writeOGR(gridSPD[i,], "./cmsp_data/", paste('AOI', gridSPD$region[i], sep='_'), driver="ESRI Shapefile", overwrite_layer=TRUE) # write the .prj file as part of the shapefile
# }
# 
# for(i in 1:length(gridSPDAK)){ # for Alaska
# 	writeOGR(gridSPDAK[i,], "./cmsp_data/", paste('AOI', gridSPDAK$region[i], sep='_'), driver="ESRI Shapefile", overwrite_layer=TRUE) # write the .prj file as part of the shapefile
# }


#######################################
## Make land and grid connection points
#######################################
# set parameters
crs <- CRS('+proj=lcc +lat_1=32 +lat_2=44 +lat_0=40 +lon_0=-96 +datum=WGS84')
crslatlong = CRS("+init=epsg:4326")

# read in files
coastlines <- readOGR(dsn=natcapfolder, layer='NAmainland_lines') # global coast layer from NatCap
#	plot(coastlines)
towns <- readOGR(dsn=paste(natcapfolder, 'ne_10m_populated_places', sep='/'), layer='ne_10m_populated_places') # populated places layer from Natural Earth data
#	plot(towns, pch=16, cex=0.2)
	# plot(towns, pch=16, cex=0.2, add=TRUE) # to add to coastlines plot

# calculate coastline length
ln <- SpatialLinesLengths(coastlines, longlat=TRUE) # returns answer in km for each line
	ln
	
# trim towns to those >1000 people
towns1000 <- towns[which(towns@data$POP_MAX >= 1000),]
	length(towns) # 7343
	length(towns1000) # 6933
	
# convert shps to planar coordinates for spatial sampling
coastlines.p <- spTransform(coastlines, CRSobj=crs)
	plot(coastlines.p, lwd=0.2, axes=TRUE)
towns.p <- spTransform(towns1000, CRSobj=crs)
#	plot(towns.p, pch=16) # odd plot: some towns seem to have extreme coordinats
	plot(towns.p, pch=16, cex=0.5, add=TRUE, col='blue') # but plots on top of NA well

# trim town to those <20km from the NA coast
nearcoast <- gWithinDistance(coastlines.p, towns.p, byid=TRUE, dist=50*1000) # slow (a few min). units in meters (50km). returns a matrix with columns corresponding to each coastline (includes a few major islands)
nearcoastany <- rowSums(nearcoast) > 0 # sum across rows: we don't care which coast a town is close to
	sum(nearcoast) # 370
	sum(nearcoastany) # 347. shows that some towns were close to multiple coastlines
	points(towns.p[which(nearcoastany),], col='red', pch=16, cex=1) # adds to the plot before

towns.nearcoast <- towns.p[which(nearcoastany),]


# add points along the coastline. 1 every km
landpts <- spsample(coastlines.p, type='regular', n=round(ln))
	plot(coastlines.p, lwd=0.2, axes=TRUE)
	points(landpts, col='red', pch=16, cex=0.2)	
	
	plot(coastlines.p, lwd=0.2, axes=TRUE, ylim=c(0, 6e5), xlim=c(-2.4e6, -2e6)) # zoom in on CA coast
	points(landpts, col='red', pch=16, cex=0.2)	

# find nearest landpt to each town
landdist <- gDistance(landpts, towns.nearcoast, byid=TRUE) # takes a couple minutes. rows are towns. cols are landpts

# project back to latlong in order to make a table for NatCap
landpts.ll <- spTransform(landpts, crslatlong)
towns.nearcoast.ll <- spTransform(towns.nearcoast, crslatlong)

# make output table of town and nearest landing point locations
towns.coords <- coordinates(towns.nearcoast.ll)
land.coords <- coordinates(landpts.ll)

n <- numeric(nrow(landdist))
outgrid <- data.frame(ID=1:nrow(landdist), LAT=towns.coords[,2], LONG=towns.coords[,1], TYPE='GRID', LOCATION=towns.nearcoast.ll@data$NAME)
outland <- data.frame(ID=1:nrow(landdist), LAT=n, LONG=n, TYPE='LAND', LOCATION=towns.nearcoast.ll@data$NAME)

for(i in 1:nrow(landdist)){ # fill in nearest landpt for each town
	j <- which.min(landdist[i,]) # index for nearest landpt
	outland$LAT[i] <- land.coords[j,2]
	outland$LONG[i] <- land.coords[j,1]
}

	# make sure it looks OK
	plot(outland$LONG, outland$LAT, pch=16, cex=0.5)
	points(outgrid$LONG, outgrid$LAT, col='red', pch=16, cex=0.5) # the towns
	
# combine town and landings points
out <- rbind(outland, outgrid)

	head(out)
	tail(out)
	
# write out
write.csv(out, file='cmsp_data/LandGridPts_NorthAmerica.csv', row.names=FALSE)