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

# project to North American Albers, but in WGS84 rather than NAD83
climgridp <- st_transform(climgrid, crs='+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=WGS84 +units=m +no_defs') # North America Albers Equal Area Conic from https://epsg.io/102008, except with WGS84
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
crs <- '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs' # North America Albers Equal Area Conic from https://epsg.io/102008
crslatlong = "+init=epsg:4326"

# read in files
coastlines <- st_read('data/natcap/NAmainland_lines.shp') # global coast layer from NatCap, trimmed to North America
#	plot(coastlines['id'], axes=TRUE) # a bit slow
townsUS <- st_read('dataDL/usgs/citiesx020_nt00007/citiesx020.shp') # US National Map cites and towns
    st_crs(townsUS) <- 4269 # set CRS to NAD83
#	plot(townsUS['NAME'], pch=16, cex=0.2, axes=TRUE)
#   plot(st_geometry(townsUS), add=TRUE, col='red') # to add to coastlines plot. Not working for some reason.
townsCan <- st_read('dataDL/statcan/lpc_000b16a_e/lpc_000b16a_e.shp') # Canadian census population centres



# project to North American Albers
coastlinesp <- st_transform(coastlines, crs = crs)
townsUSp <- st_transform(townsUS, crs = crs)
townsCanp <- st_transform(townsCan, crs = crs)
    
# find centroids for Canada
townsCanpc <- st_centroid(townsCanp)
        
    plot(st_geometry(coastlinesp), lwd=0.2, axes=TRUE)
    plot(st_geometry(townsUSp), pch=16, cex=0.1, add=TRUE, col='blue') # plots on top of NA well
    plot(st_geometry(townsCanpc), pch=16, cex=0.1, add=TRUE, col='green') # plots on top of NA well

# trim US towns to those >1000 people (Canada is already trimmed)
townsUS1000 <- townsUSp[which(townsUSp$POP_2000 >= 1000),]
	dim(townsUSp) # 35432
	dim(townsUS1000) # 9992

# combine US and Canada
names(townsCanpc)[names(townsCanpc)=='PCNAME'] <- 'NAME'
towns1000 <- rbind(townsUS1000[,c('NAME', 'geometry')], townsCanpc[,c('NAME', 'geometry')])
		
# trim town to those <50km from the NA coast
nearcoast <- st_is_within_distance(towns1000, coastlinesp, dist = 50*1000, sparse = FALSE) # slow (a few min). units in meters (50km). returns a matrix with columns corresponding to each coastline (includes a few major islands)
nearcoastany <- rowSums(nearcoast) > 0 # sum across rows: we don't care which coast a town is close to
	sum(nearcoast) # 2233
	sum(nearcoastany) # 1873. shows that some towns were close to multiple coastlines
	plot(st_geometry(towns1000[which(nearcoastany),]), col='red', pch=16, cex=1, add=TRUE) # adds to the plot before

towns1000nearcoast <- towns1000[which(nearcoastany),]


# find nearest point on the coastline to each town. So much faster than the old method of sampling!
#coastpoint.near <- st_as_sf(rgeos::gNearestPoints(as(towns1000nearcoast,"Spatial"), as(coastlinesp,"Spatial"))[2,]) # from https://gis.stackexchange.com/questions/288570/find-nearest-point-along-polyline-using-sf-package-in-r. But only returns one point. Would have to implement in a loop
towncoastlines <- st_nearest_points(towns1000nearcoast, coastlinesp) # returns LINESTRINGS from first to second geometry
towncoastlengths <- st_length(towncoastlines) # length of each line, so that we can find the closest coastline to each town
nearestptinds <- aggregate(list(ind = towncoastlengths), by = list(town = rep(1:nrow(towns1000nearcoast), each = nrow(coastlinesp))), FUN = function(x) which.min(x)) # find index of shortest line from town to a coast. Works because st_nearest_points returns a vector where y cycles fastest and x cycles slowest
nearestptinds2 <- nearestptinds$ind + seq(0, length.out = nrow(nearestptinds), by = nrow(coastlinesp)) # convert to an index into towncoastlines
pts <- st_cast(towncoastlines[nearestptinds2], "POINT") # gives all start (towns) & end (coastlines) points, alternating
coastpts <- pts[seq(2, length(pts), 2)] # just the end points (on coastlines)
	length(coastpts)
	
    #plot(st_geometry(towns1000nearcoast[1:10,]), axes = TRUE, col = 'blue') # to plot a few towns
	plot(st_geometry(towns1000nearcoast), axes = TRUE, col = 'blue', xlim = c(-2.3e6, -1.8e6), ylim = c(-4e5, 0)) # to plot a few towns
	#plot(st_geometry(towns1000nearcoast[1:100,]), axes = TRUE, col = 'blue') # to plot many towns
	plot(st_geometry(coastlinesp), add = TRUE)
    plot(st_geometry(coastpts), add = TRUE, col = 'red')

# project back to latlong in order to make a table for NatCap InVEST
coastpts.ll <- st_transform(coastpts, crs = crslatlong)
towns1000nearcoast.ll <- st_transform(towns1000nearcoast, crs = crslatlong)

# make output table of town and nearest landing point locations
# format as specified by NatCap InVEST Wave
towns.coords <- st_coordinates(towns1000nearcoast.ll)
coast.coords <- st_coordinates(coastpts.ll)

outgrid <- data.frame(ID=1:nrow(towns.coords), LAT=towns.coords[,2], LONG=towns.coords[,1], TYPE='GRID', LOCATION=towns1000nearcoast.ll$NAME)
outland <- data.frame(ID=(1+nrow(towns.coords)):(nrow(towns.coords)+nrow(coast.coords)), LAT=coast.coords[,2], LONG=coast.coords[,1], TYPE='LAND', LOCATION=towns1000nearcoast.ll$NAME)

	# make sure it looks OK
	plot(outland$LONG, outland$LAT, pch=16, cex=0.5)
	points(outgrid$LONG, outgrid$LAT, col='red', pch=16, cex=0.5) # the towns
	
    # combine town and landings points
    out <- rbind(outgrid, outland)

    head(out[out$TYPE == 'GRID', ])
    head(out[out$TYPE == 'LAND', ])
    tail(out[out$TYPE == 'GRID', ])
    tail(out[out$TYPE == 'LAND', ])
	
    # write out
    write.csv(out, file='output/landgridpts_northamerica.csv', row.names=FALSE)

# make output table of town and nearest landing point locations
# format as specified by NatCap InVEST Wind
out2 <- out[, c('ID', 'TYPE', 'LAT', 'LONG')]
names(out2)[names(out2)=='LAT'] <- 'LATI'
    head(out2)
    tail(out2)
    
    # write out
    write.csv(out2, file='output/landgridpts_northamerica_wind.csv', row.names=FALSE)
    