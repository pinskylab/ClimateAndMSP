# evaluate protected areas against shifts in species distribution
# calculate species gains, losses, turnover, etc.


############
## Flags
############

## choose which runs to use
## runtype refers to how the Species Distribution Models (SDMs) were fit
## projtype refers to how the SDM projections were done
#runtype <- 'test'; projtype=''
#runtype <- ''; projtype=''`
#runtype <- 'testK6noSeas'; projtype='_xreg'
runtype <- 'fitallreg'; projtype='_xreg'

# choose the rcp
rcp <- 85
otherrcp <- 45

# select initial and final timeperiod for these grids
periods <- c('2006-2020', '2081-2100')

####################
## helper functions
####################
require(Hmisc)
require(data.table)
require(dplyr)
require(lme4) # for mixed-effects models
require(car) # for testing ME models
require(lattice)
require(RColorBrewer)


lu <- function(x) return(length(unique(x)))

se <- function(x,na.rm=FALSE){ # standard error
	if(!na.rm){
		return(sd(x, na.rm=FALSE)/sqrt(length(x)))
	}
	if(na.rm){
		return(sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x))))
	}
}

# weighted mean function to use with summarize()
wmean <- function(x){ # values in col 1, weights in col 2
	inds <- !is.na(x[,1]) & !is.na(x[,2])
	return(weighted.mean(x[inds,1], x[inds,2])) 
}

# expects x to have columns 'periodmids' and 'rich' (in that order)
calcrichtrendmids <- function(x){
	mod <- lm(x[,2] ~ x[,1])
	return(mod$coefficients[2])
}

# turnover metrics
# expects region in column 1, lat (2), lon (3), time period in column 4, spp in col 5, pres in col 6
# periods lists the periods that can appear in col4 (in order from earliest to latest)
# only present species included
turnover <- function(x, periods){
	pers <- sort(unique(x[,4]))
	if(length(periods) != 2) stop(paste('expecting 2 time periods:', paste(unique(periods, collapse=','))))
	initcomm <- x[x[,4] == periods[1] & x[,6]==TRUE,5]
	finalcomm <- x[x[,4] == periods[2] & x[,6]==TRUE,5]
	
	a <- length(intersect(initcomm, finalcomm)) # number of shared species
	lost <- length(setdiff(initcomm, finalcomm)) # number of species lost
	gained <- length(setdiff(finalcomm, initcomm))
	return(c(nstart=length(initcomm), nend=length(finalcomm), nlost=lost, ngained=gained, flost=lost/length(initcomm), fgained=gained/length(initcomm), beta_sor=2*a/(2*a+gained+lost)))
}

# same functions, but for use with data.table
nstart <- function(period, spp, pres, periods){
	pers <- sort(unique(period))
	if(length(periods) != 2) stop(paste('expecting 2 time periods:', paste(unique(periods, collapse=','))))
	initcomm <- spp[period == periods[1] & pres==TRUE]
	return(length(initcomm))
}

nend <- function(period, spp, pres, periods){
	pers <- sort(unique(period))
	if(length(periods) != 2) stop(paste('expecting 2 time periods:', paste(unique(periods, collapse=','))))
	finalcomm <- spp[period == periods[2] & pres==TRUE]
	return(length(finalcomm))
}

nlost <- function(period, spp, pres, periods){
	pers <- sort(unique(period))
	if(length(periods) != 2) stop(paste('expecting 2 time periods:', paste(unique(periods, collapse=','))))
	initcomm <- spp[period == periods[1] & pres==TRUE]
	finalcomm <- spp[period == periods[2] & pres==TRUE]
	
	lost <- length(setdiff(initcomm, finalcomm)) # number of species lost
	return(lost)
}

ngained <- function(period, spp, pres, periods){
	pers <- sort(unique(period))
	if(length(periods) != 2) stop(paste('expecting 2 time periods:', paste(unique(periods, collapse=','))))
	initcomm <- spp[period == periods[1] & pres==TRUE]
	finalcomm <- spp[period == periods[2] & pres==TRUE]
	
	gained <- length(setdiff(finalcomm, initcomm))
	return(gained)
}

flost <- function(period, spp, pres, periods){
	pers <- sort(unique(period))
	if(length(periods) != 2) stop(paste('expecting 2 time periods:', paste(unique(periods, collapse=','))))
	initcomm <- spp[period == periods[1] & pres==TRUE]
	finalcomm <- spp[period == periods[2] & pres==TRUE]
	
	lost <- length(setdiff(initcomm, finalcomm)) # number of species lost
	return(lost/length(initcomm))
}

fgained <- function(period, spp, pres, periods){ # divide by number of spp in initial community
	pers <- sort(unique(period))
	if(length(periods) != 2) stop(paste('expecting 2 time periods:', paste(unique(periods, collapse=','))))
	initcomm <- spp[period == periods[1] & pres==TRUE]
	finalcomm <- spp[period == periods[2] & pres==TRUE]
	
	gained <- length(setdiff(finalcomm, initcomm))
	return(gained/length(initcomm))
}

fgainedalt <- function(period, spp, pres, periods){ # divide by number of spp in final community
	pers <- sort(unique(period))
	if(length(periods) != 2) stop(paste('expecting 2 time periods:', paste(unique(periods, collapse=','))))
	initcomm <- spp[period == periods[1] & pres==TRUE]
	finalcomm <- spp[period == periods[2] & pres==TRUE]
	
	gained <- length(setdiff(finalcomm, initcomm))
	return(gained/length(finalcomm))
}


beta_sor <- function(period, spp, pres, periods){
	pers <- sort(unique(period))
	if(length(periods) != 2) stop(paste('expecting 2 time periods:', paste(unique(periods, collapse=','))))
	initcomm <- spp[period == periods[1] & pres==TRUE]
	finalcomm <- spp[period == periods[2] & pres==TRUE]
	
	a <- length(intersect(initcomm, finalcomm)) # number of shared species
	lost <- length(setdiff(initcomm, finalcomm)) # number of species lost
	gained <- length(setdiff(finalcomm, initcomm))
	return(2*a/(2*a+gained+lost))
}

# pick the most common item in a list, or the first if a tie
pickone <- function(x){
	tab <- table(x)
	pick <- which.max(tab)
	return(names(tab)[pick])
}

# pick the first non-NA item in a list, or else return NA
notNA <- function(x){
	x <- x[!is.na(x)]
	if(length(x)>0){ return(as.character(x[1]))}
	if(length(x)==0){ return('NA')}
}

###########################
## Load and set up data
###########################

wdpa <- read.csv('data/wdpa_cov_by_grid0.25.csv', row.names=1) # shows wich MPAs are in which grid cells. each line is a unique grid cell-MPA combination.
	wdpa$lon[wdpa$lon<0] <- wdpa$lon[wdpa$lon<0] + 360 # convert lon to positive

	wdpa <- as.data.table(wdpa)

	# set up MPA networks
	wdpa[sub_loc == 'US-CA' & mang_auth == 'California Department of Fish and Game',network:='mlpa']
	wdpa[(sub_loc %in% c('US-DE', 'US-FL', 'US-GA', 'US-MA', 'US-MD', 'US-ME', 'US-NC', 'US-NJ', 'US-NY', 'US-RI', 'US-SC', 'US-VA')) & (lon > 280), network:='eastcoast'] # by choosing by state, this won't include federal water closures
	wdpa[(sub_loc %in% c('US-AK')), network:='ak'] # by choosing by state, this won't include federal water closures
	
	# choose one and only one region per MPA, assing to a new column
	wdpa[,region.one:=pickone(region),by=wdpapolyID]

	# calculate mpa extent
	wdpa[, ':=' (lat_min=min(lat), lat_max=max(lat)), by=wdpapolyID]


# load presmap for all model runs, a bit slow (large files)
load(paste(presmapbymodfolder, 'presmapbymod_', runtype, projtype, '_rcp', rcp, '.RData', sep='')) # loads presmap data.frame with presence/absence information from each model (slow to load).
	presmapbymod.1 <- as.data.table(presmapbymod) # seems to remove presmapbymod as well?

	setkey(presmapbymod.1,period,pres)
	presmapbymod.1 <- presmapbymod.1[.(periods,TRUE),nomatch=0] # trim to periods and pres=TRUE rows
		presmapbymod.1[,sum(is.na(pres))] # should be zero
		dim(presmapbymod.1)
	presmapbymod.1[,rcp:=rcp] # add label for rcp

load(paste(presmapbymodfolder, 'presmapbymod_', runtype, projtype, '_rcp', otherrcp, '.RData', sep='')) # loads presmap data.frame with presence/absence information from each model (slow to load) (for the other rcp, the one not used for planning)
	presmapbymod <- as.data.table(presmapbymod)

	setkey(presmapbymod,period,pres)
	presmapbymod.3 <- presmapbymod[.(periods,TRUE),nomatch=0] # trim to rows (time periods) relevant to analysi
		presmapbymod.3[,sum(is.na(pres))] # should be zero
		dim(presmapbymod.3)
	presmapbymod.3[,rcp:=otherrcp]

# bind the presmaps together
	presmapbymod <- rbind(presmapbymod.1, presmapbymod.3)
		nrow(presmapbymod.1)
		nrow(presmapbymod.3)
		dim(presmapbymod)

	rm(presmapbymod.1)
	rm(presmapbymod.3)

# trim presmap to species present, locations and time periods relevant to wdpa analysis
	regs <- sort(unique(presmapbymod$region))
	regstokeep <- setdiff(regs, c('NEFSC_NEUSFall', 'AFSC_WCTri', 'DFO_NewfoundlandSpring')) # remove some surveys to avoid duplicates in a region
	setkey(presmapbymod, region)
	thesesppbymod <- presmapbymod[.(regstokeep)] # trim out duplicate regions
	dim(thesesppbymod) # 19,859,747

	rm(presmapbymod)
	
# Add a grid ID
thesesppbymod[,gridID:=paste(lat,lon)]


#########################################################
# Calc stats for pre-defined MPA networks
# use grid cells (doesn't account for overlapping MPAs
# USE temp/wdpaturnbynetbymod.csv.gz FROM 1.2_evalWDPA.r INSTEAD?
#########################################################
climGrid <- as.data.table(read.csv('data/climGrid.csv', row.names=1)) # for temperature climatology and depth

#create list of grids in each region
thesegrids <- thesesppbymod[,.(lat, lon, region)]
setkey(thesegrids, region, lat, lon)
thesegrids <- unique(thesegrids)

# merge in bottom temperatures
setkey(climGrid, region, lat, lon)
	dim(thesegrids)
thesegrids <- climGrid[thesegrids,.(region, lat, lon, bottemp.clim.int)]
	dim(thesegrids)

# aggregate MPA networks by grid cell
wdpanet <- wdpa[,.(lat=max(lat), lon=max(lon), network=notNA(network), prop_grid=sum(prop_grid)), by=.(region, lat, lon)]

# merge in MPA networks
setkey(wdpanet, region, lat, lon)
	dim(thesegrids)
thesegrids <- wdpanet[thesegrids,.(region, lat, lon, bottemp.clim.int, network, prop_grid)]
	dim(thesegrids)

# calc size and temperature range by region
regstats <- thesegrids[,.(temprng=diff(range(bottemp.clim.int)), size=length(unique(cbind(lat,lon)))), by=region]

# calc size and temperature range of each network by region 
netstats <- thesegrids[,.(temprngnet=diff(range(bottemp.clim.int)), sizenet=sum(prop_grid)), by=.(region, network)]

# merge regional and network stats
setkey(regstats, region)
setkey(netstats, region)
stats <- netstats[regstats, .(region, network, temprng, size, temprngnet, sizenet), nomatch=0][network != 'NA',]

# calc fraction of regional ranges covered by networks
stats[,fractemp:=temprngnet/temprng]
stats[,fracsize:=sizenet/size]

stats[,.(region, network, fracsize, fractemp)]

# write out
write.csv(stats, file='cmsp_data/MPA_network_stats.csv')

#########################################################
# Calc stats for pre-defined MPA networks
# use shape files to account for overlapping MPAs
# USE temp/wdpaturnbynetbymod.csv.gz FROM 1.2_evalWDPA.r INSTEAD?
#########################################################
require(maptools)
require(rgeos)
require(rgdal)

# read in old stats (by grid square)
stats <- read.csv('cmsp_data/MPA_network_stats.csv', row.names=1)

# read in the WDPA shapes
wdpa = readShapePoly('../WDPA/NA_marine_MPA/mpinsky-search-1382225374362.shp', proj4string = CRS('+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs'))
	wdpa = wdpa[wdpa$marine == 1,] # trim to only marine (remove 20 of 3023)

# add lat and lon by MPA
getlon <- function(x){
	cds <- x@coords
	return(cds[,1])
}
avelon <- function(x){
	lons <- unlist(sapply(x@Polygons, getlon))
	return(mean(lons))
}
getlat <- function(x){
	cds <- x@coords
	return(cds[,2])
}
avelat <- function(x){
	lats <- unlist(sapply(x@Polygons, getlat))
	return(mean(lats))
}
	
wdpa$lon <- sapply(slot(wdpa, 'polygons'), avelon)
wdpa$lat <- sapply(slot(wdpa, 'polygons'), avelat)


# set up MPA networks
	wdpa$network <- NA
	wdpa$network[wdpa$sub_loc == 'US-CA' & wdpa$mang_auth == 'California Department of Fish and Game'] <- 'mlpa'
	wdpa$network[(wdpa$sub_loc %in% c('US-DE', 'US-MA', 'US-MD', 'US-ME', 'US-NC', 'US-NJ', 'US-NY', 'US-RI', 'US-VA')) & (wdpa$lon > -80) & (wdpa$lat > 35.2)] <- 'eastcoast' # by choosing by state, this won't include federal water closures
	wdpa$network[(wdpa$sub_loc %in% c('US-AK'))] <- 'ak' # by choosing by state, this won't include federal water closures

# trim only to MPAs in networks
wdpatrim <- wdpa[!is.na(wdpa$network),]

# merge all polygons from same network together
wdpamerge <- gUnaryUnion(wdpatrim, id = wdpatrim$network) # do the merge by network name

# project to Albers equal area
wdpaAlb <- spTransform(wdpamerge, CRS('+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs +datum=WGS84')) # from http://spatialreference.org/ref/esri/usa-contiguous-albers-equal-area-conic/proj4js/, except using WGS84

# calculate area by network
areas <- sapply(slot(wdpaAlb, 'polygons'), slot, 'area') # pull out area. would be even better to project to an equal-area projection...
names(areas) <- names(wdpaAlb)
areas

	# checks
		gIsValid(wdpaAlb, reason=TRUE) # fails
	wdpaAlb2 <- gBuffer(wdpaAlb, width=0, byid=TRUE) # fixes for some reason. see http://gis.stackexchange.com/questions/153905/dissolving-unifying-ill-behaved-irregular-polygons-in-r
		gIsValid(wdpaAlb2, reason=TRUE) # now good
	

## make spatial grid for each study area
# Make the grid based on corner points
	gridsz <- 0.25
	setkey(thesesppbymod, region, lat, lon)
	cds <- unique(thesesppbymod[,.(region, lat, lon)]) # pick out all unique lat lons by region
	y = numeric(5*nrow(cds)) # all the Y coords, in order
	# lower right, lower left, upper left, upper right, lower right
	for(i in 1:nrow(cds)){ y[(5*(i-1)+1):(5*(i-1)+5)] = c(cds$lat[i]-gridsz/2, cds$lat[i]-gridsz/2, cds$lat[i]+gridsz/2, cds$lat[i]+gridsz/2, cds$lat[i]-gridsz/2)}

	x = numeric(5*nrow(cds))
	# lower right, lower left, upper left, upper right, lower right: clockwise so that sp sees it as an island, not a hole
	for(i in 1:nrow(cds)){ x[(5*(i-1)+1):(5*(i-1)+5)] = c(cds$lon[i]+gridsz/2, cds$lon[i]-gridsz/2, cds$lon[i]-gridsz/2, cds$lon[i]+gridsz/2, cds$lon[i]+gridsz/2) }

	# Create a SpatialPolygonsDataFrame (a list of "Polygons", each of which is a list of "Polygon")
	pgns = vector('list', nrow(cds))
	for(i in 1:length(pgns)){
		inds2 = (5*(i-1)+1):(5*(i-1)+5)
		pgns[[i]] = Polygons(list(Polygon(cbind(x[inds2],y[inds2]))), i)
	}
	SP <- SpatialPolygons(pgns, proj4string=CRS(proj4string(wdpa)))
	SPdata <- data.frame(gridpolyID = as.numeric(sapply(slot(SP, 'polygons'), slot, 'ID')), lat = cds$lat, lon = cds$lon, region=cds$region) # would be better to put this in with SP as a SpatialPolygonsDataFrame
		length(SP)
		nrow(SPdata)
	SPdf <- SpatialPolygonsDataFrame(SP, SPdata, match.ID = FALSE)


# write out
save(SPdf, file='data/study_regions_grid_SPDF.RData')
		
# read in
load('data/study_regions_grid_SPDF.RData') # load SPdf

# merge grids by region
SPmerge <- gUnaryUnion(SPdf, id = as.character(SPdf$region))

# project grids to Albers equal area
SPAlb <- spTransform(SPmerge, CRS('+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs +datum=WGS84')) # from http://spatialreference.org/ref/esri/usa-contiguous-albers-equal-area-conic/proj4js/, except using WGS84

# calculate area by region
areasSP <- sapply(slot(SPAlb, 'polygons'), slot, 'area') # pull out area. would be even better to project to an equal-area projection...
names(areasSP) <- names(SPAlb)
areasSP

# intersect networks by region
#gI <- gIntersects(SPAlb, wdpaAlb2, byid=TRUE) # determines which regions intersection which networks
#	gI

#akAI <- gIntersection(SPAlb['AFSC_Aleutians'], wdpaAlb2['ak'], byid=FALSE) # calc the intersection just for ak network in AI
#	akAI@polygons[[1]]@area
#	akAI@polygons[[1]]@area/areasSP[1] # 8%
	
SPwdpa_int <- gIntersection(SPAlb, wdpaAlb2, byid=TRUE) # calc the intersection. takes a minute.

# calc areas of network in each region
areasSPwdpa <- data.frame(name=character(length(SPwdpa_int)), region=character(length(SPwdpa_int)), network=character(length(SPwdpa_int)), sizem2=numeric(length(SPwdpa_int)), sizenetm2=numeric(length(SPwdpa_int)))

areasSPwdpa$sizenetm2 <- sapply(slot(SPwdpa_int, 'polygons'), slot, 'area')
areasSPwdpa$name <- names(SPwdpa_int)
areasSPwdpa$region <- unlist(strsplit(areasSPwdpa$name, split=' '))[seq(1,2*nrow(areasSPwdpa),by=2)]
areasSPwdpa$network <- unlist(strsplit(areasSPwdpa$name, split=' '))[seq(2,2*nrow(areasSPwdpa),by=2)]
for(i in 1:nrow(areasSPwdpa)) areasSPwdpa$sizem2[i] <- as.numeric(areasSP[areasSPwdpa$region[i]])

areasSPwdpa

# add on to stats data.frame
stats <- merge(stats, areasSPwdpa[,c('region', 'network', 'sizem2', 'sizenetm2')], all.x=TRUE)

# calculate fraction
stats$fracsizem2 <- stats$sizenetm2/stats$sizem2

stats

# save the area calculations
write.csv(stats, file='cmsp_data/MPA_network_statsGIS.csv')



############################################
# Generate random MPA networks
# Across gradients of area and lat extent
############################################
lats <- seq(0.01,1,length.out=16) # how many steps of latitudinal extent (0 to 100%)
sizes <- seq(0.01,0.5,length.out=16) # which steps of area (could be 0 to 100%). Note that a network constrained to <100% of lat extent can't include all planning sites.
nreps <- 4 # how many repeats at each combination of lat and size

regs <- as.character(thesesppbymod[,sort(unique(region))])
nunits <- thesesppbymod[,.(nunits=length(unique(paste(lat,lon)))), by=region] # slow. calculate number of grids in each region
latext <- thesesppbymod[,.(minlat=min(lat), maxlat=max(lat), latext=max(lat)-min(lat)), by=region] # latitudinal extent by region

# to store results
randMPAs <- expand.grid(region=regs, latstep=lats, sizestep=sizes, repnum=1:nreps, latrng=NA, size=NA, meanturnIndiv=NA, netwrkturn=NA)

# a slow loop, especially on higher values of j and k
for(i in 1:length(regs)){
	print(regs[i])

	# set up the set of planning units
	setkey(thesesppbymod, region)
	punits <- thesesppbymod[regs[i],.(lat=lat,lon=lon,gridID=gridID)]
	punits <- as.data.frame(punits[!duplicated(punits),])

	le <- latext[region==regs[i],latext]

	setkey(thesesppbymod,gridID,region)

	for(j in 1:length(lats)){
		print(paste('j=',j))
		for(k in 1:length(sizes)){
			cat(paste('\tk=',k,': '))
			for(r in 1:nreps){
				cat(paste('r=',r,' '))
	
				# select a first MPA
				punits$sel <- FALSE
				firstind <- sample(1:nrow(punits),1)
				punits$sel[firstind] <- TRUE
				fl <- punits$lat[firstind] # pick a first MPA

				# calculate nearby planning units from first MPA (based on lat extent)

				punits$nearby <- FALSE # initialize "nearby" flag (within a certain fraction of lat extent)
				punits$nearby[abs(punits$lat-fl)/le <= lats[j]] <- TRUE # mark some as nearby
				
				fracwholearea <- 1/nunits[region==regs[i],nunits]
				fracnearbyarea <- 1/sum(punits$nearby)

				# add more MPAs until size criterion met or all nearby units selected
				while(fracwholearea <= sizes[k] & fracnearbyarea <= 1 & sum(punits$nearby & !punits$sel)>0 ) {
					newunit <- sample(which(punits$nearby & !punits$sel),1) # select a new unit
					punits$sel[newunit] <- TRUE # assign as selected
					fracwholearea <- sum(punits$sel)/nunits[region==regs[i],nunits]
					fracnearbyarea <- sum(punits$sel)/sum(punits$nearby)
					punits$nearby <- ((punits$lat - min(punits$lat[punits$sel]))/le <= lats[j] & punits$lat >= min(punits$lat[punits$sel])) | ((max(punits$lat[punits$sel]) - punits$lat)/le < lats[j] & punits$lat <= max(punits$lat[punits$sel]))
				}
				
				# basic calcs for this random network
				ind <- with(randMPAs, which(region==regs[i] & latstep==lats[j] & sizestep==sizes[k] & repnum==r))
				randMPAs$latrng[ind] <- with(punits[punits$sel,], max(lat)-min(lat))/le
				randMPAs$size[ind] <- sum(punits$sel)/nunits[region==regs[i],nunits]
				
				# evaluate ecological turnover for each MPA
				randMPAs$meanturnIndiv[ind] <- thesesppbymod[.(punits$gridID[punits$sel], regs[i]),.(beta_sor=beta_sor(period, sppocean, pres, periods)), by="gridID,rcp,model"][,mean(beta_sor)]

				# collapse to a network
				thesesppbymodnet <- thesesppbymod[.(punits$gridID[punits$sel], regs[i]),.(pres=max(pres)), by="sppocean,period,rcp,model"] #  aggregate function within the selected grids
				randMPAs$netwrkturn[ind] <- thesesppbymodnet[,.(beta_sor=beta_sor(period, sppocean, pres, periods)), by="rcp,model"][,mean(beta_sor)]

			}			
		}
	}
}

# write out
write.csv(randMPAs, file=paste('cmsp_data/randMPAs_', Sys.Date(), '.csv', sep=''))



## PLOTS

#xyplot(netwrkturn ~ latrng | region, data=randMPAs)
#xyplot(netwrkturn ~ size | region, data=randMPAs)
#
#latrngshingle <- equal.count(randMPAs$latrng, number=4, overlap=0.25)
#xyplot(netwrkturn ~ size | latrngshingle + region, data=randMPAs)
#
#sizeshingle <- equal.count(randMPAs$size, number=4, overlap=0)
#xyplot(netwrkturn ~ latrng | sizeshingle + region, data=randMPAs)

# plot netwrkturn as color dots
#colrmp <- colorRamp(brewer.pal(11, name='Spectral'))
#
#par(mfrow=c(3,3))
#for(i in 1:length(regs)){
#	inds <- randMPAs$region==regs[i]
#	plot(randMPAs$size[inds], randMPAs$latrng[inds], pch=16, col=rgb(colrmp(randMPAs$netwrkturn[inds]), maxColorValue=256), main=regs[i])
#}
#

# plot netwrkturn as averages within grid squares
colrmp <- colorRamp(brewer.pal(11, name='Spectral'))

szs <- seq(min(randMPAs$size, na.rm=TRUE), max(randMPAs$size, na.rm=TRUE), length.out=20)
rngs <- seq(min(randMPAs$latrng, na.rm=TRUE), max(randMPAs$latrng, na.rm=TRUE), length.out=20)
szstep <- diff(szs)[1]
latstep <- diff(rngs)[1]

gridave <- expand.grid(region=regs, size=szs, latrng=rngs, ave=NA)
for(i in 1:nrow(gridave)){
	inds <- randMPAs$region==gridave$region[i] & abs(randMPAs$latrng - gridave$latrng[i]) < latstep/2 & abs(randMPAs$size - gridave$size[i]) < szstep
	gridave$ave[i] <- mean(randMPAs$netwrkturn[inds])		
}

par(mfrow=c(3,3))
for(i in 1:length(regs)){
	inds <- gridave$region==regs[i] & !is.na(gridave$ave)
	newz <- gridave$ave[inds] - min(gridave$ave[inds])
	newz <- newz/max(newz)
	plot(gridave$size[inds], gridave$latrng[inds], pch=15, cex=1.5, col=rgb(colrmp(newz), maxColorValue=256), main=regs[i])
}

# a loess smooth to plot a levelplot
#library(geoR)
#
#par(mfrow=c(3,3))
#for(j in 1:length(regs)){
#	i <- randMPAs$region == regs[j]
#	loess = loess(netwrkturn ~ size*latrng, data = randMPAs[i,], span=0.4)
#	szs <- seq(min(randMPAs$size, na.rm=TRUE), max(randMPAs$size, na.rm=TRUE), length.out=50)
#	rngs <- seq(min(randMPAs$latrng, na.rm=TRUE), max(randMPAs$latrng, na.rm=TRUE), length.out=50)
#	fit = expand.grid(list(size = szs, latrng = rngs))
#	z <- predict(loess, newdata=fit)
#
#	# NA in blanks
#	szstep <- diff(lats)[1]
#	latstep <- diff(sizes)[1]
#	for(k in 1:nrow(fit)){
#		inds <- randMPAs$region==regs[j] & abs(randMPAs$latrng - fit$latrng[k]) < latstep & abs(randMPAs$size - fit$size[k]) < szstep
#		if(sum(inds)==0){
#			z[k] <- NA
#		}
#	}
#
#	image(szs, rngs, z, main=regs[j])
#	#points(randMPAs$size[i], randMPAs$latrng[i])
#}
#




##################################################
##################################################
# Generate random MPA networks
# Across gradients of area and temperature extent
##################################################
##################################################
climGrid <- read.csv('data/climGrid.csv', row.names=1) # for temperature climatology and depth

temps <- seq(0.01,1,length.out=10) # how many steps of temperature extent (0 to 100%)
sizes <- seq(0.01,0.5,length.out=10) # which steps of area (could be 0 to 100%). Note that a network constrained to <100% of lat extent can't include all planning sites.
nreps <- 10 # how many repeats at each combination of lat and size

regs <- as.character(thesesppbymod[,sort(unique(region))])

# to store results
randMPAs <- expand.grid(region=regs, tempstep=temps, sizestep=sizes, repnum=1:nreps, temprng=NA, size=NA, meanturnIndiv=NA, netwrkturn=NA)

# a slow loop, especially on higher values of j and k
for(i in 1:length(regs)){
	print(regs[i])

	# set up the set of planning units
	setkey(thesesppbymod, region)
	punits <- thesesppbymod[regs[i],.(lat=lat,lon=lon,gridID=gridID)]
	punits <- as.data.frame(punits[!duplicated(punits),])
	oldrow <- nrow(punits)
	punits <- merge(punits, climGrid[climGrid$region==regs[i], c('lat', 'lon', 'bottemp.clim.int')], all.x=TRUE)
	if(nrow(punits) != oldrow) print(paste('i=', i, ': rows changed during merge with clim'))
	names(punits)[names(punits)=='bottemp.clim.int'] <- 'bt'

	tx <- diff(range(punits$bt, na.rm=TRUE))
	nu <- length(unique(punits$gridID))

	setkey(thesesppbymod,gridID,region)

	for(j in 1:length(temps)){
		cat(paste('\nj=',j, sep=''))
		for(k in 1:length(sizes)){
			cat(paste('\tk=',k,': ', sep=''))
			for(r in 1:nreps){
				cat(paste('r=',r,' '))
	
				# select a first MPA
				punits$sel <- FALSE
				firstind <- sample(1:nrow(punits),1)
				punits$sel[firstind] <- TRUE
				ft <- punits$bt[firstind] # pick a first MPA

				# calculate nearby planning units from first MPA (based on temp extent)
				punits$nearby <- FALSE # initialize "nearby" flag (within a certain fraction of lat extent)
				punits$nearby[abs(punits$bt-ft)/tx <= temps[j]] <- TRUE # mark some as nearby
				
				fracwholearea <- 1/nu
				fracnearbyarea <- 1/sum(punits$nearby)

				# add more MPAs until size criterion met or all nearby units selected
				while(fracwholearea <= sizes[k] & fracnearbyarea <= 1 & sum(punits$nearby & !punits$sel)>0 ) {
					newunit <- sample(which(punits$nearby & !punits$sel),1) # select a new unit
					punits$sel[newunit] <- TRUE # assign as selected
					fracwholearea <- sum(punits$sel)/nu
					fracnearbyarea <- sum(punits$sel)/sum(punits$nearby)
					punits$nearby <- ((punits$bt - min(punits$bt[punits$sel]))/tx <= temps[j] & punits$bt >= min(punits$bt[punits$sel])) | ((max(punits$bt[punits$sel]) - punits$bt)/tx < temps[j] & punits$bt <= max(punits$bt[punits$sel])) # measure from the min temp if new potential temps are > min, and measure from max temp if new potential temps are < max
				}
				
				# basic calcs for this random network
				ind <- with(randMPAs, which(region==regs[i] & tempstep==temps[j] & sizestep==sizes[k] & repnum==r))
				randMPAs$temprng[ind] <- with(punits[punits$sel,], max(bt)-min(bt))/tx
				randMPAs$size[ind] <- sum(punits$sel)/nu
				
				# evaluate ecological turnover for each MPA
				randMPAs$meanturnIndiv[ind] <- thesesppbymod[.(punits$gridID[punits$sel], regs[i]),.(beta_sor=beta_sor(period, sppocean, pres, periods)), by="gridID,rcp,model"][,mean(beta_sor)]

				# collapse to a network
				thesesppbymodnet <- thesesppbymod[.(punits$gridID[punits$sel], regs[i]),.(pres=max(pres)), by="sppocean,period,rcp,model"] #  aggregate function within the selected grids
				randMPAs$netwrkturn[ind] <- thesesppbymodnet[,.(beta_sor=beta_sor(period, sppocean, pres, periods)), by="rcp,model"][,mean(beta_sor)]

			}			
		}
	}
	cat('\n')
}

# write out
write.csv(randMPAs, file=paste('cmsp_data/randMPAs_byBT_', Sys.Date(), '.csv', sep=''))


#################
## PLOTS
#################
randMPAs1 <- read.csv('cmsp_data/randMPAs_byBT_2016-07-19.csv', row.names=1) # 10x10x5
randMPAs2 <- read.csv('cmsp_data/randMPAs_byBT_2016-07-20.csv', row.names=1) # 10x10x10
stats <- read.csv('cmsp_data/MPA_network_statsGIS.csv', row.names=1)

# combine the data
randMPAs <- rbind(randMPAs1, randMPAs2)

regs <- sort(unique(randMPAs$region))


#xyplot(netwrkturn ~ latrng | region, data=randMPAs)
#xyplot(netwrkturn ~ size | region, data=randMPAs)
#
#latrngshingle <- equal.count(randMPAs$latrng, number=4, overlap=0.25)
#xyplot(netwrkturn ~ size | latrngshingle + region, data=randMPAs)
#
#sizeshingle <- equal.count(randMPAs$size, number=4, overlap=0)
#xyplot(netwrkturn ~ latrng | sizeshingle + region, data=randMPAs)

# plot netwrkturn as color dots
colrmp <- colorRamp(brewer.pal(11, name='Spectral'))

par(mfrow=c(3,3))
for(i in 1:length(regs)){
	inds <- randMPAs$region==regs[i]
	plot(randMPAs$size[inds], randMPAs$temprng[inds], pch=16, col=rgb(colrmp(randMPAs$netwrkturn[inds]), maxColorValue=256), main=regs[i])
}


# plot netwrkturn as averages within grid squares
colrmp <- colorRamp(brewer.pal(11, name='Spectral'))
colpal <- colorRampPalette(brewer.pal(11, name='Spectral'))
xlabs <- c('', '', '', '', '', '', 'Fraction of region in network', 'Fraction of region in network', 'Fraction of region in network')
ylabs <- c('Fraction of thermal\nrange in network', '', '', 'Fraction of thermal\nrange in network', '', '', 'Fraction of thermal\nrange in network', '', '') 
regs <- c('AFSC_EBS', 'AFSC_Aleutians', 'AFSC_GOA', 'NWFSC_WCAnn', 'SEFSC_GOMex', 'NEFSC_NEUSSpring', 'DFO_ScotianShelf', 'DFO_SoGulf', 'DFO_NewfoundlandFall') # set plot order
#regsnice = c('EBS', 'AI', 'GoA', 'WC', 'GoM', 'Neast', 'SS', 'SGoSL', 'Newf')
regsnice = c('Eastern Bering Sea', 'Aleutian Islands', 'Gulf of Alaska', 'West Coast U.S.', 'Gulf of Mexico', 'Northeast U.S.', 'Scotian Shelf', 'So. Gulf of St. Lawrence', 'Newfoundland')

szs <- seq(min(randMPAs$size, na.rm=TRUE), max(randMPAs$size, na.rm=TRUE), length.out=10)
rngs <- seq(min(randMPAs$temprng, na.rm=TRUE), max(randMPAs$temprng, na.rm=TRUE), length.out=10)
szstep <- diff(szs)[1]
tempstep <- diff(rngs)[1]

gridave <- expand.grid(region=regs, size=szs, temprng=rngs, ave=NA)
for(i in 1:nrow(gridave)){
	inds <- randMPAs$region==gridave$region[i] & abs(randMPAs$temprng - gridave$temprng[i]) < tempstep/2 & abs(randMPAs$size - gridave$size[i]) < szstep
	gridave$ave[i] <- mean(randMPAs$netwrkturn[inds])		
}

	# using boxes
#quartz(width=6, height=7)
## pdf(width=6, height=7, file='cmsp_figures/randMPAs_byBT.pdf')
#par(mfrow=c(3,3), mgp=c(2,1,0))
#for(i in 1:length(regs)){
#	inds <- gridave$region==regs[i] & !is.na(gridave$ave)
#	newz <- gridave$ave[inds] - min(gridave$ave[inds])
#	newz <- newz/max(newz)
#	plot(gridave$size[inds], gridave$temprng[inds], pch=15, cex=2.5, col=rgb(colrmp(newz), maxColorValue=256), main=regs[i], xlab=xlabs[i], ylab=ylabs[i])
#
#	nas <-  gridave$region==regs[i] & is.na(gridave$ave)
#	points(gridave$size[nas], gridave$temprng[nas], pch=15, cex=2.3, col='grey')
#	text(x=0.4, y=0.05, labels=paste(signif(range(gridave$ave[inds]),2), collapse='-'))
#
#	# add dot for empirical network
#	inds2 <- as.character(stats$region)==regs[i]
#	if(sum(inds2)>0){
#		points(stats$fracsize[inds2], stats$fractemp[inds2], pch=10, cex=2, col='white')
#	}
#}

dev.off()

	# using image
quartz(width=6, height=6)
# pdf(width=6, height=6, file='cmsp_figures/randMPAs_byBT.pdf')
# jpeg(units='in', res=300, width=6, height=6, file='cmsp_figures/randMPAs_byBT.jpg')
par(mfrow=c(3,3), mgp=c(2,0.5,0), mai=c(0.2, 0.3, 0.3, 0.1), omi=c(0.4, 0.4, 0, 0), xpd=NA, tcl=-0.3, las=1)
for(i in 1:length(regs)){
	inds <- gridave$region==regs[i] & !is.na(gridave$ave)
	thisdat <- gridave[inds,]
	thisdat$newz <- thisdat$ave - min(thisdat$ave)
	thisdat$newz <- thisdat$newz/max(thisdat$newz)
	mat <- as.data.frame(dcast.data.table(as.data.table(thisdat), temprng ~ size, value.var='newz'))
	row.names(mat) <- mat$temprng
	mat <- t(as.matrix(mat[,2:ncol(mat)])) # transpose because of the way image handles matrices
	image(z=mat, x=sort(unique(thisdat$size)), y=sort(unique(thisdat$temprng)), col=colpal(100), main=regsnice[i], xlab=xlabs[i], ylab=ylabs[i])

	matna <- mat
	matna[is.na(mat)] <- 1
	matna[!is.na(mat)] <- NA
	image(z=matna, x=sort(unique(thisdat$size)), y=sort(unique(thisdat$temprng)), col='grey', add=TRUE)

	text(x=0.34, y=-0.01, labels=paste('Similarity:', paste(signif(range(gridave$ave[inds]),2), collapse='-')))
	
	# add dot for empirical network
	inds2 <- as.character(stats$region)==regs[i]
	if(sum(inds2)>0 & regs[i] != 'DFO_ScotianShelf'){
		points(stats$fracsizem2[inds2], stats$fractemp[inds2], pch=10, cex=2, col='white')
	}
}

dev.off()


# a loess smooth to plot a levelplot
#library(geoR)
#
#par(mfrow=c(3,3))
#for(j in 1:length(regs)){
#	i <- randMPAs$region == regs[j]
#	loess = loess(netwrkturn ~ size*latrng, data = randMPAs[i,], span=0.4)
#	szs <- seq(min(randMPAs$size, na.rm=TRUE), max(randMPAs$size, na.rm=TRUE), length.out=50)
#	rngs <- seq(min(randMPAs$latrng, na.rm=TRUE), max(randMPAs$latrng, na.rm=TRUE), length.out=50)
#	fit = expand.grid(list(size = szs, latrng = rngs))
#	z <- predict(loess, newdata=fit)
#
#	# NA in blanks
#	szstep <- diff(lats)[1]
#	latstep <- diff(sizes)[1]
#	for(k in 1:nrow(fit)){
#		inds <- randMPAs$region==regs[j] & abs(randMPAs$latrng - fit$latrng[k]) < latstep & abs(randMPAs$size - fit$size[k]) < szstep
#		if(sum(inds)==0){
#			z[k] <- NA
#		}
#	}
#
#	image(szs, rngs, z, main=regs[j])
#	#points(randMPAs$size[i], randMPAs$latrng[i])
#}
#
