# Make a shapefile of our Areas of Interest, for use in the NatCap's InVest software

## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	}
# could add code for Lauren's working directory here

#######################
## Load libraries
#######################
require(maptools)
require(RColorBrewer)
require(rgdal)
require(PBSmapping)
require(rgeos) # for gBuffer


##############################
## Read in and process data
##############################
climgrid <- read.csv('data/climGrid.csv', row.names=1)
#land <- readOGR(dsn='cmsp_data/', layer='global_polygon') # global land layer from NatCap
#	plot(land) # takes too long

# label by area of interest regions
#climgrid$AOI <- 'northamerica'
climgrid$AOI <- NA
climgrid$AOI[climgrid$region %in% c('AFSC_Aleutians', 'AFSC_EBS', 'AFSC_GOA')] <- 'Alaska'
#climgrid$AOI[climgrid$region %in% c('AFSC_Aleutians')] <- 'Alaska_Aleutians'
#climgrid$AOI[climgrid$region %in% c('AFSC_Aleutians') & climgrid$lon > 183] <- 'Alaska_Aleutians_east'
#climgrid$AOI[climgrid$region %in% c('AFSC_Aleutians') & climgrid$lon >= 177 & climgrid$lon <= 183] <- 'Alaska_Aleutians_middle'
#climgrid$AOI[climgrid$region %in% c('AFSC_Aleutians') & climgrid$lon < 177] <- 'Alaska_Aleutians_west'
#climgrid$AOI[climgrid$region %in% c('AFSC_EBS')] <- 'Alaska_EBS'
#climgrid$AOI[climgrid$region %in% c('AFSC_GOA')] <- 'Alaska_GOA'
climgrid$AOI[climgrid$region %in% c('AFSC_WCTri', 'NWFSC_WCAnn', 'SEFSC_GOMex', 'DFO_NewfoundlandFall', 'DFO_NewfoundlandSpring', 'DFO_ScotianShelf', 'DFO_SoGulf', 'NEFSC_NEUSFall', 'NEFSC_NEUSSpring')] <- 'notAlaska'
#climgrid$AOI[climgrid$region %in% c('AFSC_WCTri', 'NWFSC_WCAnn', )] <- 'WestCoast'
#climgrid$AOI[climgrid$region %in% c()] <- 'GoMex'
#climgrid$AOI[climgrid$region %in% c('DFO_NewfoundlandFall', 'DFO_NewfoundlandSpring', 'DFO_ScotianShelf', 'DFO_SoGulf', 'NEFSC_NEUSFall', 'NEFSC_NEUSSpring')] <- 'Northeast'
sum(is.na(climgrid$AOI))
sort(unique(climgrid$AOI))

# draw 1/4deg boxes around all grid cell centers
PIDs <- rep(1:nrow(climgrid), rep(4,nrow(climgrid))) # polygon ids
POSs <- rep(1:4, nrow(climgrid)) # vertex ID within a polygon
Xs <- rep(climgrid$lon, rep(4, nrow(climgrid))) + c(-0.125, 0.125, 0.125, -0.125)
Ys <- rep(climgrid$lat, rep(4, nrow(climgrid))) + c(-0.125, -0.125, 0.125, 0.125) # lower left, upper left, upper right, lower right

gridPS <- as.PolySet(data.frame(PID = PIDs, POS=POSs, X=Xs, Y=Ys), projection='LL')
	plotPolys(gridPS) # appears to cover a somewhat wider area than we project to. that's OK
gridSP <- PolySet2SpatialPolygons(gridPS) # convert to spatial polygon for union

# merge all into multipart polygons based on AOI
gridSP.m <- unionSpatialPolygons(gridSP, climgrid$AOI)	
	plot(gridSP.m, lwd=0.2, axes=TRUE, col=1:length(gridSP.m))
#	plot(gridSP.m, axes=TRUE, ylim=c(40,70), xlim=c(-170,-60)) # doesn't work
	
# project into meters using Lambert Conformal Conic
gridSP.p <- spTransform(gridSP.m, CRS('+proj=lcc +lat_1=32 +lat_2=44 +lat_0=40 +lon_0=-96'))
	plot(gridSP.p, lwd=0.2, axes=TRUE, col=1:4)

# buffer
gridSP.b <- gBuffer(gridSP.p, width=100000, byid=TRUE) # 100km
	plot(gridSP.b, lwd=0.2, axes=TRUE, col=1:4)
	
# convert to spatialpolygonsdataframe for writing
pid <- sapply(slot(gridSP.b, "polygons"), function(x) slot(x, "ID"))
gridSPD <- SpatialPolygonsDataFrame(gridSP.b, data.frame(region=pid, row.names=pid))

# split Alaska along 179 and 181deg longitude (3 pieces)
maskPS <- as.PolySet(data.frame(PID=rep(c(1:3),rep(4,3)), POS=rep(1:4, 3), X=rep(c(150,179,181,220), c(2,4,4,2)), Y=rep(c(45,65,45,65,45,65,45), c(1,2,2,2,2,2,1))), projection='LL')
mask <- PolySet2SpatialPolygons(maskPS)
mask.t <- spTransform(mask, proj4string(gridSPD))
	#plot(mask.t, col=1:length(mask.t))
	plot(gridSP.b, lwd=0.2, axes=TRUE, col=1:4)
	plot(mask.t, add=TRUE)
gridSPAK <- gIntersection(gridSPD[grep('Alaska', gridSPD$region),], mask.t, byid=TRUE)
	plot(gridSPAK, col=1:length(gridSPAK))

	# convert Alaska to spatialpolygonsdataframe for writing
	pid <- sapply(slot(gridSPAK, "polygons"), function(x) slot(x, "ID"))
	regs <- c('Alaska_west', 'Alaska_middle', 'Alaska_east')
	gridSPDAK <- SpatialPolygonsDataFrame(gridSPAK, data.frame(region=regs, row.names=pid))

# write out each region as a separate shapefile
#writePolyShape(merged3, fn='cmsp_data/AOI_northamerica') # won't write .prj
notAK <- grep('Alaska', gridSPD$region, invert=TRUE)
for(i in notAK){ # for all except Alaska
	writeOGR(gridSPD[i,], "./cmsp_data/", paste('AOI', gridSPD$region[i], sep='_'), driver="ESRI Shapefile", overwrite_layer=TRUE) # write the .prj file as part of the shapefile
}
for(i in 1:length(gridSPDAK)){ # for Alaska
	writeOGR(gridSPDAK[i,], "./cmsp_data/", paste('AOI', gridSPDAK$region[i], sep='_'), driver="ESRI Shapefile", overwrite_layer=TRUE) # write the .prj file as part of the shapefile
}
