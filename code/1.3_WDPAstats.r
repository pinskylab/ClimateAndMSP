# stats about impacts on protected areas from shifting species
# also some basic stats on the area of MPAs and the area and range of temperatures spanned by MPA networks


############
## Flags
############



####################
## helper functions
####################
# require(Hmisc)
require(data.table)
# require(lme4) # for mixed-effects models
# require(car) # for testing ME models


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


###########################################
# Examine MPA size relative to grid size
###########################################
wdpagrid <- fread('gunzip -c output/wdpa_cov_by_grid0.05.csv.gz', drop = 1) # shows which MPAs are in which grid cells. each line is a unique grid cell-MPA combination.


# plot vs. grid size
wdpagrid[!duplicated(WDPA_PID), hist(log10(area_fullwdpa), xlab = 'log10(area in m2)', main = 'MPA vs. grid cell sizes')]
abline(v = log10(unique(wdpagrid$area_grid)), col = '#FF000001')

# number > grid size
wdpagrid[!duplicated(WDPA_PID), .N] # 925
wdpagrid[!duplicated(WDPA_PID), sum(area_fullwdpa >= area_grid)] # 332
wdpagrid[!duplicated(WDPA_PID), sum(area_fullwdpa < area_grid)] # 593
wdpagrid[!duplicated(WDPA_PID), sum(area_fullwdpa < area_grid)/.N] # 64%


#########################################################
# Calc stats for pre-defined MPA networks
# use grid cells (doesn't account for overlapping MPAs
# USE temp/wdpaturnbynetbymod.csv.gz FROM 1.2_evalWDPA.r INSTEAD?
# MOVED FROM 2.0_evalRandomNetworks.r
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
# MOVED FROM 2.0_evalRandomNetworks.r
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



#############################################
# Examine turnover within MPAs
# Examine all climate models individually
#############################################
# read in turnover data
wdpaturnbyMPAbymod <- fread('gunzip -c temp/wdpaturnbyMPAbymod.csv.gz', drop=1) # Turnover in each MPA. Don't read the row numbers

ntk <- wdpaturnbyMPAbymod[,NO_TAKE %in% c('All', 'Part')] # no take reserves (index into wdpaturnbyMPAbymod)
sum(ntk) # 74 no take reserves
wdpaturnbyMPAbymod[ntk, sort(SUB_LOC)]
wdpaturnbyMPAbymod[ntk, sum(grepl('US-CA', SUB_LOC))] # 49 no-take in California

# calculate Sorenson turnover and fractional change
for (r in c(26, 85)) {
    for (m in 1:18) {
        wdpaturnbyMPAbymod[, (paste0('beta_sor.', r, '.', m)) := 
                               2*get(paste0('nshared.', r, '.', m)) /
                               (2*get(paste0('nshared.', r, '.', m)) + get(paste0('ngained.', r, '.', m)) + get(paste0('nlost.', r, '.', m)))]
        wdpaturnbyMPAbymod[, (paste0('flost.', r, '.', m)) := 
                               get(paste0('nlost.', r, '.', m)) /
                               get(paste0('ninit.', r, '.', m))]
        wdpaturnbyMPAbymod[, (paste0('fgained.', r, '.', m)) := 
                               get(paste0('ngained.', r, '.', m)) /
                               get(paste0('nfinal.', r, '.', m))]
    }
}

# convert to long format
wdpalong <- melt(wdpaturnbyMPAbymod, id.vars = c('WDPA_PID', 'NAME', 'MANG_AUTH', 'SUB_LOC', 'lat_min', 'lat_max', 'lon_min', 'lon_max', 'area_fullwdpa', 'network'), 
                 measure.vars = patterns('26|85'), variable.name = 'measure', value.name = 'val')
wdpalong[, c('measure', 'rcp', 'model') := tstrsplit(measure, '.', fixed = TRUE)] # split name of turnover measure apart from RCP and model #
    dim(wdpalong) # 266400 x 14

# read in change in temp by region
# trend1 <- read.csv('data/climTrendbyreg_rcp45.csv', row.names=1)
# 	trend1$rcp <- 45
# trend2 <- read.csv('data/climTrendbyreg_rcp85.csv', row.names=1)
# 	trend2$rcp <- 85
# nms <- c('region', 'rcp', 'delta_surf', 'delta_bott')
# trend <- as.data.table(rbind(trend1[,nms], trend2[,nms]))


# Fraction of species lost (fraction of original community)
    # all MPAs
wdpalong[measure == 'flost', .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','WDPA_PID')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models with MPAs/RCPs, then average across MPAs

wdpalong[measure == 'flost', .(mean = mean(val,na.rm = TRUE), lat = mean(c(lat_min, lat_max)), pac = (lon_min< -100)[1]), 
         by = c('rcp','WDPA_PID')][, plot(lat, mean, col = ifelse(rcp == 26, ifelse(pac, 'light blue', 'pink'),
                                                       ifelse(pac, 'blue', 'red')))] # plot MPA ensemble mean vs. lat. Ocean by color (Pac is blue, Atl is red), RCPs by shade (26 is light, 85 is dark.


	# no-take
wdpalong[ntk & measure == 'flost', .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','WDPA_PID')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models with MPAs/RCPs, then average across MPAs


# Fraction of species gained (fraction of final community)
	# all MPAs
wdpalong[measure == 'fgained', .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','WDPA_PID')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models with MPAs/RCPs, then average across MPAs

wdpalong[measure == 'fgained', .(mean = mean(val,na.rm = TRUE), lat = mean(c(lat_min, lat_max)), pac = (lon_min< -100)[1]), 
         by = c('rcp','WDPA_PID')][, plot(lat, mean, 
                                          col = ifelse(rcp == 26, ifelse(pac, 'light blue', 'pink'),
                                                       ifelse(pac, 'blue', 'red')))] # plot MPA ensemble mean vs. lat. Ocean by color (Pac is blue, Atl is red), RCPs by shade (26 is light, 85 is dark.


	# no-take
wdpalong[ntk & measure == 'fgained', .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','WDPA_PID')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models with MPAs/RCPs, then average across MPAs

                                   
# Similarity
	# all MPAs
wdpalong[measure == 'beta_sor', .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','WDPA_PID')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models with MPAs/RCPs, then average across MPAs

wdpalong[measure == 'beta_sor', .(mean = mean(val,na.rm = TRUE), lat = mean(c(lat_min, lat_max)), pac = (lon_min< -100)[1]), 
         by = c('rcp','WDPA_PID')][, plot(lat, mean, 
                                          col = ifelse(rcp == 26, ifelse(pac, 'light blue', 'pink'),
                                                       ifelse(pac, 'blue', 'red')))] # plot MPA ensemble mean vs. lat. Ocean by color (Pac is blue, Atl is red), RCPs by shade (26 is light, 85 is dark.


	# no-take
wdpalong[ntk & measure == 'beta_sor', .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','WDPA_PID')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models with MPAs/RCPs, then average across MPAs

	# examine similarity vs. MPA size
wdpalong[measure == 'beta_sor', cor.test(area_fullwdpa, val), by = c('rcp', 'model')]

wdpalong[measure == 'beta_sor', .(beta_sor = mean(val,na.rm = TRUE), log10size = mean(log10(area_fullwdpa))), 
         by = c('rcp','WDPA_PID')][, plot(log10size, beta_sor, 
                                          col = ifelse(rcp == 26, 'blue', 'red'))] # plot MPA ensemble mean vs. size. RCP by color.

				
	# examine similarity by MPA latitudinal range
wdpalong[measure == 'beta_sor', cor.test(lat_max - lat_min, val), by = c('rcp', 'model')]

wdpalong[measure == 'beta_sor', .(beta_sor = mean(val,na.rm = TRUE), latsize = mean(log10(lat_max - lat_min + 0.1))), 
         by = c('rcp','WDPA_PID')][, plot(latsize, beta_sor, 
                                          col = ifelse(rcp == 26, 'blue', 'red'))] # plot MPA ensemble mean vs. size. RCP by color.

		
# examine similarity by MPA size and lat rng
wdpalong[measure == 'beta_sor', .(beta_sor = mean(val, na.rm = TRUE), latsize = mean(log10(lat_max - lat_min + 0.1)), log10size = mean(log10(area_fullwdpa))), 
         by = c('rcp','WDPA_PID')][, summary(mod <- lm(beta_sor~rcp*latsize*log10size))] # linear model with 3-way interaction


		require(interplot) # to plot the interaction
		interplot(m=mod, var1='latsize', var2='log10size') +
			xlab('log(area) scaled') + 
			ylab('Coefficient for scaled latitudinal range')

		interplot(m=mod, var2='latrngsc', var1='lrep_areasc') +
			ylab('Coefficient for log(area) scaled') + 
			xlab('Latitudinal range scaled')



################################################################
# Examine change within MPA networks and the component MPAs
# Across each model in the ensemble
################################################################
wdpaturnbynetbymod <- fread('gunzip -c temp/wdpaturnbynetbymod.csv.gz') # network results
# also need wdpalong from above section
		
# Size of networks
	wdpaturnbyMPAbymod[!is.na(network),.(num_mpas=length(unique(WDPA_PID))), by=network]

# calculate Sorenson turnover and fractional change
for (r in c(26, 85)) {
    for (m in 1:18) {
        wdpaturnbynetbymod[, (paste0('beta_sor.', r, '.', m)) := 
                               2*get(paste0('nshared.', r, '.', m)) /
                               (2*get(paste0('nshared.', r, '.', m)) + get(paste0('ngained.', r, '.', m)) + get(paste0('ngained.', r, '.', m)))]
        wdpaturnbynetbymod[, (paste0('flost.', r, '.', m)) := 
                               get(paste0('nlost.', r, '.', m)) /
                               get(paste0('ninit.', r, '.', m))]
        wdpaturnbynetbymod[, (paste0('fgained.', r, '.', m)) := 
                               get(paste0('ngained.', r, '.', m)) /
                               get(paste0('nfinal.', r, '.', m))]
    }
}

# convert to long format
wdpanetlong <- melt(wdpaturnbynetbymod, id.vars = c('network', 'lat', 'lon', 'area'), 
                 measure.vars = patterns('26|85'), variable.name = 'measure', value.name = 'val')
wdpanetlong[, c('measure', 'rcp', 'model') := tstrsplit(measure, '.', fixed = TRUE)] # split name of turnover measure apart from RCP and model #
    dim(wdpanetlong) # 864 x 8

# Fraction of species lost
	# all MPA networks
wdpanetlong[measure == 'flost', .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','network')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models within networks x RCPs, then average across networks

	# within the individual MPAs of these networks
wdpalong[measure == 'flost' & !is.na(network), .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','WDPA_PID')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models within MPAs x RCPs, then average across MPAs

# Fraction of species gained (fraction of final community)
	# all MPA networks
wdpanetlong[measure == 'fgained', .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','network')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models within networks x RCPs, then average across networks

	# within the individual MPAs of these networks
wdpalong[measure == 'fgained' & !is.na(network), .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','WDPA_PID')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models within MPAs x RCPs, then average across MPAs

# Similarity (Sorenson)
	# all MPA networks
wdpanetlong[measure == 'beta_sor', .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','network')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models within networks x RCPs, then average across networks

	# within the individual MPAs of these networks
wdpalong[measure == 'beta_sor' & !is.na(network), .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','WDPA_PID')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models within MPAs x RCPs, then average across MPAs

