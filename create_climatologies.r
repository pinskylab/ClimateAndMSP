###############################################################################
## Create climatologies for each region: gridded averages with interpolation ##
###############################################################################
setwd('/Users/mpinsky/Documents/Rutgers/Range projections/')
load("Output/trawl_allregionsforprojections_2015-02-02.RData") # loads dat

# remove NA lat/lon
	dat = dat[!is.na(dat$lat) & !is.na(dat$lon),]

# output lat/lons from trawl tows for Lauren Roger to match against benthic habitat data
	outlatlons = dat[!duplicated(dat[,c('lat', 'lon', 'depth')]),c('lat', 'lon', 'depth')]
	outlatlons = outlatlons[order(outlatlons$lat, outlatlons$lon, outlatlons$depth),]
	write.csv(outlatlons, file=paste('Output/trawl_latlons_forLauren_', Sys.Date(), '.csv', sep=''), row.names=FALSE)

# Fix lon to only positive to match CMIP5
	dat$lon[dat$lon < 0] = dat$lon[dat$lon < 0] + 360 # fix lons to only positive to match CMIP5 data
	
# Prep temperature data
	temp = dat[complete.cases(dat[,c('lat', 'lon', 'depth')]),c('region', 'stratum', 'lat', 'lon', 'depth', 'year', 'yearsurv', 'month', 'bottemp', 'surftemp')]
		temp = droplevels(temp)
		rm(dat)

	# Add grid indicator for the climatology
	gridsize=0.25 # size of grid of the climate data, in degrees
	temp$latgrid = floor(temp$lat/gridsize)*gridsize + gridsize/2 # round to nearest grid center
	temp$longrid = floor(temp$lon/gridsize)*gridsize + gridsize/2
	
	# Trim to no later than December 16, 2005 (to match historical period in CMIP5)
	temp = temp[temp$year<2006,]
	
	# Label by decade within each region (using start date within region as start of the decades)
	regs = sort(unique(temp$region))
	temp$decade = NA
	for(i in 1:length(regs)){
		inds = temp$region == regs[i]
		yrs = sort(unique(temp$year[inds]))
		decs = seq(min(yrs), max(yrs), by=10) + 5 # center of each decade
		temp$decade[inds] = floor((temp$yearsurv[inds] - min(decs)+5)/10)*10+min(decs)
	}

	# Label by season
	temp$season = NA
	inds = temp$month <= 3
	temp$season[inds] = 1
	inds = temp$month > 3 & temp$month<=6
	temp$season[inds] = 2
	inds = temp$month > 6 & temp$month<=9
	temp$season[inds] = 3
	inds = temp$month > 9
	temp$season[inds] = 4

	# Remove Newfoundland SST data, since there is so little
	temp$surftemp[temp$region %in% c('DFO_NewfoundlandSpring', 'DFO_NewfoundlandFall')] = NA

	
# Average by grid cell and by decade within seasons. lat, lon are now grid centers
	climdec = aggregate(list(bottemp.clim = temp$bottemp, surftemp.clim = temp$surftemp, depth=temp$depth), by=list(lat = temp$latgrid, lon = temp$longrid, decade = temp$decade, region=temp$region, season=temp$season), FUN=mean, na.rm=TRUE)
		dim(climdec) # 20,228
		names(climdec)

# Average across decades (done as a second step to avoid one decade with lots of data swamping the average)
	clim = aggregate(list(bottemp.clim = climdec$bottemp.clim, surftemp.clim = climdec$surftemp.clim, depth=climdec$depth), by=list(lat = climdec$lat, lon = climdec$lon,  region=climdec$region, season=climdec$season), FUN=mean) # don't remove NAs because it would mean dropping a decade
		dim(clim) # 9988
		names(clim)
	
		# plot each region in a separate figure. color axis is relative within each region
		require(lattice)
		cols = colorRampPalette(colors = c('blue', 'white', 'red'))
		pdf(width=30, height=12, file=paste('Figures/climBT_grid_', Sys.Date(), '.pdf', sep=''))
		levelplot(bottemp.clim ~ lon*lat|region*season, scales = list(x='free', y='free'), data=clim, at = seq(0,1, length.out=40), col.regions=cols, panel=function(...,z, subscripts, at){
				panel.fill(col='light grey')
				panel.levelplot(..., z=z, subscripts=subscripts, at=seq(min(z[subscripts], na.rm=TRUE), max(z[subscripts], na.rm=TRUE), length.out=20))})

		dev.off()

		pdf(width=30, height=12, file=paste('Figures/climSST_grid_', Sys.Date(), '.pdf', sep=''))
		levelplot(surftemp.clim ~ lon*lat|region*season, scales = list(x='free', y='free'), data=clim, at = seq(0,1, length.out=40), col.regions=cols, panel=function(...,z, subscripts, at){
				panel.fill(col='light grey')
				panel.levelplot(..., z=z, subscripts=subscripts, at=seq(min(z[subscripts], na.rm=TRUE), max(z[subscripts], na.rm=TRUE), length.out=20))})

		dev.off()

# Choose the seasons within each region to interpolate
table(clim$region, clim$season)

	# keep region/seasons that have enough data for decent interpolation
keep = c('AFSC_Aleutians 2', 'AFSC_Aleutians 3', 'AFSC_EBS 2', 'AFSC_EBS 3', 'AFSC_GOA 3', 'AFSC_WCTri 3', 'DFO_NewfoundlandFall 4', 'DFO_NewfoundlandSpring 2', 'DFO_ScotianShelf 1', 'DFO_ScotianShelf 3', 'DFO_ScotianShelf 4', 'DFO_SoGulf 3', 'NEFSC_NEUSFall 4', 'NEFSC_NEUSSpring 1', 'NEFSC_NEUSSpring 2', 'NWFSC_WCAnn 3', 'SEFSC_GOMex 2', 'SEFSC_GOMex 3')

clim = clim[paste(clim$region, clim$season) %in% keep,] # trim down to these region/season combinations
	dim(clim) # 8046
climold = clim # save a copy

	# Scotian Shelf winter extends far south. Trim this off by only keep gridcells present in at least two seasons
	sskeepgrids = clim[clim$region == 'DFO_ScotianShelf' & !duplicated(clim[,c('region', 'lat', 'lon', 'season')]), c('region', 'lat', 'lon', 'season')] # list of Scotian Shelf gridcells and the season in which they appear
	sskeepgrids2 = reshape(sskeepgrids, direction='wide', timevar='season', v.names='season', idvar = c('region', 'lat', 'lon')) # reshape to wide to easily see which grids cells appear in which seasons
	sskeepgrids3 = sskeepgrids2[rowSums(!is.na(sskeepgrids2[,c('season.1', 'season.3', 'season.4')]))>1,] # only keep those grids present in at least 2 seasons (e.g., not only winter)

	clim = clim[clim$region != 'DFO_ScotianShelf' | (paste(clim$lat, clim$lon) %in% paste(sskeepgrids3$lat, sskeepgrids3$lon)),] # trim out the Scotian Shelf grids that didn't make the cut
	
	nrow(climold)-nrow(clim) # dropped 165 rows
	nrow(sskeepgrids2)-nrow(sskeepgrids3) # also dropped 165 rows from Scotian Shelf: perfect (to make sure I didn't drop anything extra)


# Interpolate temperatures
	# useful functions from http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
	deg2rad <- function(deg) return(deg*pi/180)
	# Haversine formula (hf)
	gcd.hf <- function(long1, lat1, long2, lat2) {
		long1 = deg2rad(long1); long2 = deg2rad(long2); lat1 = deg2rad(lat1); lat2=deg2rad(lat2)
		R <- 6371 # Earth mean radius [km]
		delta.long <- (long2 - long1)
		delta.lat <- (lat2 - lat1)
		a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
		c <- 2 * asin(pmin(1,sqrt(a)))
		return(d = R * c) # Distance in km
	}

	# make sure all grid cells appear in each season
	regs = sort(unique(clim$region))
	for(i in 1:length(regs)){
		grids = clim[!duplicated(clim[,c('lat', 'lon')]) & clim$region == regs[i], c('lat', 'lon', 'depth')] # will arbitrarily pull the first depth value for each grid cell
		seasons = sort(unique(clim$season[clim$region == regs[i]]))
		for(j in seasons){ # for each season
			thesegrids = clim[!duplicated(clim[,c('lat', 'lon', 'season')]) & clim$region == regs[i] & clim$season == j, c('lat', 'lon')]
			newrows = grids[!(paste(grids$lat, grids$lon) %in% paste(thesegrids$lat, thesegrids$lon)),]
			if(nrow(newrows)>0){
				newrows$region = regs[i]
				newrows$season = j
				newrows$bottemp.clim = NA
				newrows$surftemp.clim = NA
				clim = rbind(clim, newrows)
			}
		}
	}
	dim(clim) # 8633

	inds = which(is.na(clim$bottemp.clim)) # points to fill for BT
	clim$bottemp.clim.int = clim$bottemp.clim # interpolated BT
	for(i in 1:length(inds)){
		lt = clim$lat[inds[i]]
		ln = clim$lon[inds[i]]
		nn = clim$lat %in% c(lt, lt-gridsize, lt+gridsize, lt-gridsize*2, lt+gridsize*2, lt-gridsize*3, lt+gridsize*3) & clim$lon %in% c(ln, ln-gridsize, ln+gridsize, ln-gridsize*2, ln+gridsize*2, ln-gridsize*3, ln+gridsize*3) & clim$region == unique(clim$region[inds[i]]) & clim$season == unique(clim$season[inds[i]]) # nearest neighbors from same region and season (including the focal grid cell)
		dist = pmax(1, gcd.hf(ln, lt, clim$lon[nn], clim$lat[nn])) # distances in km, but set 1km as a floor (so no zeros)
		clim$bottemp.clim.int[inds[i]] = weighted.mean(x=clim$bottemp.clim[nn], w=1/dist, na.rm=TRUE) # weight values from 3 grid cells in every direction by 1/distance
	}	
		sum(is.na(clim$bottemp.clim.int)) # = 251: some missing values
		unique(clim$region[is.na(clim$bottemp.clim.int)]) # EBS, NEUSSpring, GOMex
		table(clim$region, is.na(clim$bottemp.clim.int))

	inds = which(is.na(clim$surftemp.clim) & !(clim$region %in% c('DFO_NewfoundlandSpring', 'DFO_NewfundlandFall', 'NWFSC_WCAnn'))) # points to fill for SST
	clim$surftemp.clim.int = clim$surftemp.clim # interpolated ST
	for(i in 1:length(inds)){
		lt = clim$lat[inds[i]]
		ln = clim$lon[inds[i]]
		nn = clim$lat %in% c(lt, lt-gridsize, lt+gridsize, lt-gridsize*2, lt+gridsize*2, lt-gridsize*3, lt+gridsize*3) & clim$lon %in% c(ln, ln-gridsize, ln+gridsize, ln-gridsize*2, ln+gridsize*2, ln-gridsize*3, ln+gridsize*3) & clim$region == unique(clim$region[inds[i]])  & clim$season == unique(clim$season[inds[i]]) # nearest neighbors (including the focal grid cell)
		dist = pmax(1,gcd.hf(ln, lt, clim$lon[nn], clim$lat[nn])) # distances in km (>=1)
		clim$surftemp.clim.int[inds[i]] = weighted.mean(x=clim$surftemp.clim[nn], w=1/dist, na.rm=TRUE)
	}	
		sum(is.na(clim$surftemp.clim.int)) # = 2373: many missing SST values
		unique(clim$region[is.na(clim$surftemp.clim.int)]) # Newfoundland Fall and Spring and West Coast Annual still missing, plus DFO_SoGulf. That's what we'd expect based on missing SST data. Also EBS, GOA, NEUSSpring, GOMex.
		table(clim$region, is.na(clim$surftemp.clim.int))


	inds = which(is.na(clim$depth)) # points to fill for depth (0)
		sum(inds)

	# plot each region in a separate figure (interpolated temp)
		#BT
	require(lattice)
	cols = colorRampPalette(colors = c('blue', 'white', 'red'))
	pdf(width=30, height=12, file=paste('Figures/climBT_grid_interp_', Sys.Date(), '.pdf', sep=''))
	levelplot(bottemp.clim.int ~ lon*lat|region*season, scales = list(x='free', y='free'), data=clim, at = seq(0,1, length.out=40), col.regions = cols, panel=function(...,z, subscripts, at){
			panel.fill(col='light grey')
			panel.levelplot(..., z=z, subscripts=subscripts, at=seq(min(z[subscripts], na.rm=TRUE), max(z[subscripts], na.rm=TRUE), length.out=20))})
	dev.off()

		#SST
	pdf(width=30, height=12, file=paste('Figures/climSST_grid_interp_', Sys.Date(), '.pdf', sep=''))
	levelplot(surftemp.clim.int ~ lon*lat|region*season, scales = list(x='free', y='free'), data=clim, at = seq(0,1, length.out=40), col.regions = cols, panel=function(...,z, subscripts, at){
			panel.fill(col='light grey')
			panel.levelplot(..., z=z, subscripts=subscripts, at=seq(min(z[subscripts], na.rm=TRUE), max(z[subscripts], na.rm=TRUE), length.out=20))})

	dev.off()

		#Depth
	pdf(width=10, height=6, file=paste('Figures/climDepth_grid_interp_', Sys.Date(), '.pdf', sep=''))
	levelplot(depth ~ lon*lat|region, scales = list(x='free', y='free'), data=clim, at = seq(0,1, length.out=40), col.regions = terrain.colors, panel=function(...,z, subscripts, at){
			panel.fill(col='light grey')
			panel.levelplot(..., z=z, subscripts=subscripts, at=quantile(z[subscripts], probs = seq(0,1, length.out=40), na.rm=TRUE))})

	dev.off()

# Add grid indicator for merging with climate deltas (1deg instead of 1/4deg for climatology)
	climgridsize=1 # size of grid of the climate data, in degrees
	clim$latgrid = floor(clim$lat/climgridsize)*climgridsize + climgridsize/2 # round to nearest grid center
	clim$longrid = floor(clim$lon/climgridsize)*climgridsize + climgridsize/2

## Add stratum based on majority within the grid (a slow process since line-by-line, 15+ min)
#	clim$stratum = character(nrow(clim))
#	ties = numeric(0)
#	for(i in 1:nrow(clim)){
#		if(i %% 100 == 0) print(i)
#		strats = temp$stratum[temp$latgrid == clim$lat[i] & temp$longrid == clim$lon[i] & temp$region == clim$region[i]]
#		t = sort(table(strats), decreasing=TRUE) # sort by frequency
#		if(length(t)>1){ # check for ties
#			if(t[1] == t[2]){
#				print(paste('there was a tie for i=', i))
#				ties = c(ties, i)
#			}
#		}
#		clim$stratum[i] = names(t[1])
#		
#		# this would be for nearest neighbor instead
#		#dist = gcd.hf(clim$lon[i], clim$lat[i], temp$lon, temp$lat) # distances in km to all hauls
#		#clim$stratum[i] = as.character(temp$stratum[which.min(dist)])
#	}
#	length(ties) # number of grid cells with ties between which stratum had the most points (34)
#	sum(is.na(clim$stratum)) # 0

# write out climatology
	write.csv(clim, file=paste('Output/climGrid_', Sys.Date(), '.csv', sep=''))
	
# write out all grid cell lat/lons for Lauren Rogers
	outlatlons2 = clim[!duplicated(clim[,c('lat', 'lon')]),c('lat', 'lon', 'depth')]
	outlatlons2 = outlatlons2[order(outlatlons2$lat, outlatlons2$lon),]
	write.csv(outlatlons2, file=paste('Output/projectiongrid_latlons_forLauren_', Sys.Date(), '.csv', sep=''), row.names=FALSE)
