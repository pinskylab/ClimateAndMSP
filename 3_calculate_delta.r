#################################
## Read netCDF files from IPCC ##
#################################
setwd('/Users/mpinsky/Documents/Rutgers/Range projections/')
require(ncdf4) # for reading ncdf files
require(chron) # for converting julian days

nyrs = function(x) as.numeric(as.character(years(x))) # return year from a chron object
mo = function(x) 1+as.numeric(strftime(x, format='%m')) # return month
jmo2mo = function(x, startmonth=1){ # converts from julian month (starting at 1 for startmonth in some year) to month 1:12
	temp = (x %% 12) + startmonth-1
	temp[temp>12] = temp[temp>12]-12
	temp[temp==0] = 12
	return(temp)
}
su = function(x) sort(unique(x)) # sort unique values from a vector
meanbyyear = function(x, yr) aggregate(list(x=x), by=list(yr=yr), FUN=mean)$x # take mean of x by levels of yr and return x
linearsmooth = function(x, yr=1:length(x)){ # smooth a time-series with a linear regresion
	mod = lm(x ~ yr)
	return(predict(mod, new.data = data.frame(yr = seq(min(yr), max(yr)))))
}
linearb = function(x, yr=1:length(x)){ # return slope of x vs. yrs. require at least 3 datapoints
	if(sum(!is.na(x))>2){
		mod = lm(as.numeric(x) ~ yr)
		return(coef(mod)[2])
	} else {
		return(NA)
	}
}
agmean = function(x, by) aggregate(list(x), by=by, FUN=mean)[,2] # aggregate and take the mean of a vector. NOTE: much faster to reshape so that I can use straight up mean
season = function(x){ # convert month of year to season
	out = rep(NA, length(x))
	out[x <= 3 & x>0] = 1
	out[x <= 6 & x>3] = 2
	out[x <= 9 & x>6] = 3
	out[x <= 12 & x>9] = 4
	return(out)
}
roundto = function(x,y){r = which.min(abs(x-y)); return(y[r])} # rounds x to nearest number in y

load("Output/trawl_allregionsforprojections_2014-11-03.RData") # loads dat to get base years for each region
clim = read.csv('Output/climGrid_2014-12-10.csv', row.names=1, stringsAsFactors=FALSE); type='Grid'

# Fix dat to match CMIP5 format
	dat = dat[!is.na(dat$lat) & !is.na(dat$lon),] # remove NA lat/lon
	dat$lon[dat$lon < 0] = dat$lon[dat$lon < 0] + 360 # fix lons to only positive to match climate data

# Get range of years, lats, and lons for each region. Use Jan 1960 as the base date (month 1)
	dat$dates = chron(dates. = paste('01', dat$month, dat$year, sep='/'), format='d/m/y') # don't have day information in dat
	baseyrs = aggregate(list(rng = dat$year), by=list(region = dat$region), FUN=range)
	basedates = aggregate(list(min = dat$dates), by=list(region = dat$region), FUN=min)
		basedates = merge(basedates, aggregate(list(max = dat$dates), by=list(region = dat$region), FUN=max))
		basedates$minmo = (nyrs(basedates$min)-1960)*12+mo(basedates$min) # months since Jan 1960
		basedates$maxmo = (nyrs(basedates$max)-1960)*12+mo(basedates$max)
	basemonths = aggregate(list(months = mo(dat$dates)), by=list(region = dat$region), FUN=su) # list the months
	baselats = aggregate(list(rng = dat$lat), by=list(region = dat$region), FUN=range)
		baselats$rng = cbind(floor(baselats$rng[,1])+0.5, ceiling(baselats$rng[,2])-0.5) # round to nearest degree center
	baselons = aggregate(list(rng = dat$lon), by=list(region = dat$region), FUN=range)
		baselons$rng = cbind(floor(baselons$rng[,1])+0.5, ceiling(baselons$rng[,2])-0.5) # round to nearest degree center
	basedepths = aggregate(list(rng = dat$depth), by=list(region = dat$region), FUN=range, na.rm=T)
		basedepths$rng[,1] = 0 # to make sure I get all the way to the surface
	regs = sort(unique(dat$region))

	rm(dat)

	# write out
	write.csv(basedates, file=paste('Output/basedates_', Sys.Date(), '.csv', sep=''))
	write.csv(baselats, file=paste('Output/baselats_', Sys.Date(), '.csv', sep=''))
	write.csv(baselons, file=paste('Output/baselons_', Sys.Date(), '.csv', sep=''))
	write.csv(basedepths, file=paste('Output/baseddepths_', Sys.Date(), '.csv', sep=''))


# Read in, average, and write out each model's historical period for each region (average)
agency = c('CNRM-CERFACS', 'IPSL', 'IPSL', 'MOHC', 'MPI-M', 'MPI-M', 'MRI', 'NCAR', 'NCC', 'NCC', 'NOAA-GFDL', 'NOAA-GFDL', 'NOAA-GFDL')
model = c('CNRM-CM5', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR', 'HadGem2-CC', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MRI-CGCM3', 'CCSM4', 'NorESM1-M', 'NorESM1-ME', 'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M')

for(j in 1:length(agency)){
	file = paste('Froelicher_data/', agency[j], '/', model[j], '/historical/thetao_all_regrid_malin.nc', sep='')
	print(paste(agency[j], model[j]))
	n = nc_open(file) # Open the netCDF file
	lats = ncvar_get(n, 'LAT')
	lons = ncvar_get(n, 'LON')
	depths = ncvar_get(n, 'LEV')
	times = 1:552 # Jan 1960 to Dec 2005
	if(model[j]=='HadGem2-CC') times = 1:551 # the exception, ends Nov 2005

	# Variable to hold the climate data in the base time period for each region
	basetemps = vector('list', length(regs))
		names(basetemps) = regs

	# Loop through the regions
	for(i in 1:length(regs)){
		regtimes = basedates$minmo[i]:basedates$maxmo[i]
		xtimes = intersect(times, regtimes) # the intersecting times
		start = c(which.min(abs(lons - baselons$rng[i,1])), which.min(abs(lats - baselats$rng[i,1])), which.min(abs(depths - basedepths$rng[i,1])), xtimes[1]) # indices of the start (lon, lat, depth, time)
		end = c(which.min(abs(lons - baselons$rng[i,2])), which.min(abs(lats - baselats$rng[i,2])), which.min(abs(depths - basedepths$rng[i,2])), xtimes[length(xtimes)]) # indices of the end
		temp = ncvar_get(n, 'THETAO_REGRID', start = start, count=end - start+1) # a slice through time and space (lon, lat, depth, time)
		theselons = ncvar_get(n, 'LON', start=start[1], count= end[1]-start[1]+1)
		theselats = ncvar_get(n, 'LAT', start=start[2], count= end[2]-start[2]+1)
		thesedepths = ncvar_get(n, 'LEV', start=start[3], count= end[3]-start[3]+1)
			print(paste(round(min(thesedepths)), 'm', 1960+floor((xtimes[1]-1)/12), '/', jmo2mo(xtimes[1]), 1960+floor((xtimes[length(xtimes)]-1)/12), '/', jmo2mo(xtimes[length(xtimes)]))) # make sure I got the depth and date conversions right. 
	
		# Average across years within each lat/lon/depth/season
		dim(temp)
		temp2 = apply(temp, MARGIN=c(1,2,3), FUN=agmean, by=list(season(jmo2mo(xtimes)))) 
		temp2 = aperm(temp2, c(2, 3, 4, 1)) # move season dimension from first to last
		print(dim(temp2))

		basetemps[[i]] = temp2
		dimnames(basetemps[[i]]) = list(lon=theselons, lat=theselats, depth=thesedepths, season = 1:4)
	}

	# Write out
	outfile = paste('Output/tempshist_', agency[j], '_', model[j], '.RData', sep='')
	save(basetemps, file=outfile)
	
	# close
	nc_close(n)
}


# Identify depth of bottom grid cell in each model: is more than one cell the "bottom" at each lat/lon point? YES (it turns out)
#	agency = c('CNRM-CERFACS', 'IPSL', 'IPSL', 'MOHC', 'MPI-M', 'MPI-M', 'MRI', 'NCAR', 'NCC', 'NCC', 'NOAA-GFDL', 'NOAA-GFDL', 'NOAA-GFDL')
#	model = c('CNRM-CM5', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR', 'HadGem2-CC', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MRI-CGCM3', 'CCSM4', 'NorESM1-M', 'NorESM1-ME', 'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M')
#	regs = sort(unique(clim$region))

	# add depthgrid to clim (nearest depthgrid relevant to each climate model)
	# mark each gridcell that we'll use (sst or bt)
#	for(i in 1:length(agency)){ # for each climate model
#		load(paste('Output/tempshist_', agency[i], '_', model[i], '.RData', sep=''))
#		botindstokeep = vector('list', length(basetemps)) # array with 1 where we will need a value for the bottom temp calculation
#
#		dgnm = paste('depthgrid', i, sep='') # different grid for each model
#		clim[[dgnm]] = NA # to hold depth rounded to the GCM grid
#		for(j in 1:length(regs)){ # for each region
#			inds = clim$region == regs[j]
#			regind = which(names(basetemps) == regs[j]) # index into array for this region (likely same as j)
#
#			#add depth grid specific to each region and each model
#			deps = sort(unique(as.numeric(dimnames(basetemps[[regind]])$depth))) # find the depths used by each model
#			clim[[dgnm]][inds] = sapply(clim$depth[inds], FUN=roundto, y=deps) # round and save
#
#			# create an empty botindstokeep for this region
#			botindstokeep[[regind]] = array(0, dim=dim(basetemps[[regind]]))
#
#			# mark all surface to keep
#			#indstokeep[[regind]][,,1,] = 1
#
#			# mark bottom gridcells as to keep where we will need it ###### workign here####
#			dims = dimnames(basetemps[[regind]])
#			for(k in 1:nrow(clim[inds,])){ # cycle through each row of clim for this region
#				lonind = which(as.numeric(dims$lon) == clim$longrid[inds][k]) # indices into basetemp for this lon/lat/depth (all seasons)
#				latind = which(as.numeric(dims$lat) == clim$latgrid[inds][k]) 
#				depind = which(as.numeric(dims$depth) == clim[[dgnm]][inds][k])
#				lookwider = FALSE # flag to look at adjacent cells
#				if(length(lonind)>0 & length(latind)>0 & length(depind)>0){ # make sure the cell exists in the array
#
#					if(!is.na(basetemps[[regind]][lonind, latind, depind, 1])){ # test if this gridcell has data in the climate model
#						botindstokeep[[regind]][lonind, latind, depind, ] = 1 # if so, mark it as to keep for all seasons (assumes all seasons have data)
#					} else { # look elsewhere in the water column at same lat/lon for an existing value
#						if(any(!is.na(basetemps[[regind]][lonind, latind, , 1]))){
#							jj = which(!is.na(basetemps[[regind]][lonind, latind, , 1])) # find depths that have a temperature value
#							kk = which.min(abs(as.numeric(dims$depth)[jj] - clim[[dgnm]][inds][k]))
#							botindstokeep[[regind]][lonind, latind, jj[kk], ] = 1 #mark value closest to desired depth as "keep"
#						} else {
#							lookwider= TRUE # if no value exists in the water column here, look at a wider set of cells
#						}
#					}
#				}
#				if(lookwider | length(lonind)==0 | length(latind)==0 | length(depind)==0){ # need to look at 3x3 grid of lat/lons for climate values
#					loninds = which(abs(as.numeric(dims$lon) - clim$longrid[inds][k])<=1) # indices into basetemp within 1 degree of this lon
#					latinds = which(abs(as.numeric(dims$lat) - clim$latgrid[inds][k])<=1) # indices into basetemp within 1 degree of this lat
#					if(any(!is.na(basetemps[[regind]][loninds, latinds, depind, 1]))){ # look in adjacent cells at same depth for climate values
#						botindstokeep[[regind]][loninds, latinds, depind, ] = 1 #mark values as "keep"
#					} else {
#						if(any(!is.na(basetemps[[regind]][loninds, latinds, , 1]))){ # look in adjacent cells at any depth
#							for(jj in loninds){ # cycle through all nine cells <= 1 step away
#								for(kk in latinds){
#									iii = which(!is.na(as.numeric(basetemps[[regind]][jj, kk, , 1]))) # find depths that have a temperature value
#									jjj = which.min(abs(as.numeric(dims$depth)[iii] - clim[[dgnm]][inds][k])) # find depth closest to desired depth
#									botindstokeep[[regind]][jj, kk, iii[jjj], ] = 1 #mark value closest to desired depth as "keep" in this cell		
#								}
#							}
#						}	else {
#							print(paste(model[i], regs[j], 'k=', k, 'lat=', clim$lat[inds][k], 'lon=', clim$lon[inds][k], ': missing all <=1 deg'))
#						}									
#					}
#				}
#			}
#		}	
#	}
#
#
#	# does every marked lat/lon have only depth marked?
#	sums = vector('list', length(botindstokeep))
#	for(i in 1:length(sums)){
#		sums[[i]] = apply(botindstokeep[[i]][,,,1], MARGIN=c(1,2), FUN=sum)
#	}
#	sums # NO: multiple depths marked in each lat/lon cell. so I can't simply trim to the single depth that corresponds to true depth


# Read in and write out each model's 2020-2100 period RCP8.5 (by year and season)
agency = c('CNRM-CERFACS', 'IPSL', 'IPSL', 'MOHC', 'MPI-M', 'MPI-M', 'MRI', 'NCAR', 'NCC', 'NCC', 'NOAA-GFDL', 'NOAA-GFDL', 'NOAA-GFDL')
model = c('CNRM-CM5', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR', 'HadGem2-CC', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MRI-CGCM3', 'CCSM4', 'NorESM1-M', 'NorESM1-ME', 'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M')

for(j in 1:length(agency)){ # very slow because of apply() down below. takes 36 hours on Amphiprion.
	print(j)
	file = paste('Froelicher_data/', agency[j], '/', model[j], '/rcp85/thetao_all_regrid_malin.nc', sep='')
	print(paste(agency[j], model[j]))
	n = nc_open(file) # Open the netCDF file
	lats = ncvar_get(n, 'LAT')
	lons = ncvar_get(n, 'LON')
	depths = ncvar_get(n, 'LEV')
	startmonth = 1; startyear = 2006
	if(model[j]=='HadGem2-CC'){	startmonth = 12; startyear = 2005} # HadGEM2-CC model is different for RCP8.5 (see readme.txt), starts Dec 2005
	times = 1:length(ncvar_get(n, 'TIME')) # time, measured in months since the start month/year for this model

	# Variable to hold the climate data in the base time period for each region
	temps2100 = vector('list', length(regs))
		names(temps2100) = regs

	# times to extract for each model: Jan 2020 to Dec 2100 = 960 months. measured in months since the model's start month/year
	keeptimes = seq(from=(2020-startyear)*12+1-startmonth+1, length.out=960) # time, measured in months since model's start month/year
	xtimes = intersect(times, keeptimes) # the intersecting times
	if(keeptimes[1] != xtimes[1]) stop(paste('j=',j, 'did not start in January'))

	# Loop through the regions
	for(i in 1:length(regs)){
		start = c(which.min(abs(lons - baselons$rng[i,1])), which.min(abs(lats - baselats$rng[i,1])), which.min(abs(depths - basedepths$rng[i,1])), xtimes[1]) # indices of the start (lon, lat, depth, time)
		end = c(which.min(abs(lons - baselons$rng[i,2])), which.min(abs(lats - baselats$rng[i,2])), which.min(abs(depths - basedepths$rng[i,2])), xtimes[length(xtimes)]) # indices of the end

		temp = ncvar_get(n, 'THETAO_REGRID', start = start, count=end - start+1) # a slice of time: dims are lon, lat, depth, time
		theselons = ncvar_get(n, 'LON', start=start[1], count= end[1]-start[1]+1)
		theselats = ncvar_get(n, 'LAT', start=start[2], count= end[2]-start[2]+1)
		thesedepths = ncvar_get(n, 'LEV', start=start[3], count= end[3]-start[3]+1)

		# reshape so that month, season and year are in separate dimensions
		dim(temp)
		temp2 = array(temp, dim=c(dim(temp)[1:3], 3, 4, floor(length(xtimes)/12))) # lon, lat, depth, month, season, year
		dim(temp2)

		# Average across months within each lat/lon/depth/season/year
		dim(temp2)
		temp3 = apply(temp2, MARGIN=c(1,2,3,5,6), FUN=mean)
		print(dim(temp3))

		temps2100[[i]] = temp3
		dimnames(temps2100[[i]]) = list(lon=theselons, lat=theselats, depth=thesedepths, yr = 2020:(2019+floor(length(xtimes)/12)), season=1:4)
	}

	# Write out
	outfile = paste('Output/tempsRCP85_2020-2100_', agency[j], '_', model[j], '.RData', sep='')
	save(temps2100, file=outfile)

	# close
	nc_close(n)
}


# Read in, average, and write out each model's control run drift (as degrees per year)
agency = c('CNRM-CERFACS', 'IPSL', 'IPSL', 'MOHC', 'MPI-M', 'MPI-M', 'MRI', 'NCAR', 'NCC', 'NCC', 'NOAA-GFDL', 'NOAA-GFDL', 'NOAA-GFDL')
model = c('CNRM-CM5', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR', 'HadGem2-CC', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MRI-CGCM3', 'CCSM4', 'NorESM1-M', 'NorESM1-ME', 'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M')

for(j in 1:length(agency)){
	file = paste('Froelicher_data/', agency[j], '/', model[j], '/piControl/thetao_all_regrid_malin.nc', sep='')
	print(paste(agency[j], model[j]))
	n = nc_open(file) # Open the netCDF file
	lats = ncvar_get(n, 'LAT')
	lons = ncvar_get(n, 'LON')
	depths = ncvar_get(n, 'LEV')
	startmonth = 1; startyear = 2006; seasonlist = c(1,1,1,2,2,2,3,3,3,4,4,4) # seasons Jan-Mar, Apr-Jun, Jul-Sep, Oct-Dec
	if(model[j]=='HadGem2-CC'){	 # HadGEM2-CC model is different for RCP8.5 (see readme.txt), starts Dec 2005
		startmonth = 12; startyear = 2005
		seasonlist = c(4,1,1,1,2,2,2,3,3,3,4,4)
	}
	times = 1:length(ncvar_get(n, 'TIME')) # time, measured in months since the start month/year for this model

	# Variable to hold the smoothed climate data for the control runs (linear regression)
	tempscontrol = vector('list', length(regs))
		names(tempscontrol) = regs

	# Loop through the regions
	for(i in 1:length(regs)){
		start = c(which.min(abs(lons - baselons$rng[i,1])), which.min(abs(lats - baselats$rng[i,1])), which.min(abs(depths - basedepths$rng[i,1])), 1) # indices of the start (lon, lat, depth, time)
		end = c(which.min(abs(lons - baselons$rng[i,2])), which.min(abs(lats - baselats$rng[i,2])), which.min(abs(depths - basedepths$rng[i,2])), start[4]-2) # indices of the end

		temp = ncvar_get(n, 'THETAO_REGRID', start = start, count=end - start+1) # a slice of time
		theselons = ncvar_get(n, 'LON', start=start[1], count= end[1]-start[1]+1)
		theselats = ncvar_get(n, 'LAT', start=start[2], count= end[2]-start[2]+1)
		thesedepths = ncvar_get(n, 'LEV', start=start[3], count= end[3]-start[3]+1)

		# reshape so that month and year are in separate dimensions
		dim(temp)
		temp2 = array(temp, dim=c(dim(temp)[1:3], 3, 4, floor(length(xtimes)/12))) # lon, lat, depth, month, season, year
		dim(temp2)

		# Take slope of temperature across years within each lat/lon/depth/season
		dim(temp2)
		temp3 = apply(temp2, MARGIN=c(1,2,3,5), FUN=linearb, yr=rep(1:dim(temp2)[6], each=3))
		print(dim(temp3))
	
		tempscontrol[[i]] = temp3
		dimnames(tempscontrol[[i]]) = list(lon=theselons, lat=theselats, depth=thesedepths, season=1:4)
	}

	# Write out
	outfile = paste('Output/tempscontrol_', agency[j], '_', model[j], '.RData', sep='')
	save(tempscontrol, file=outfile)

	# close
	nc_close(n)
}

# next: read in rcp45 runs for 2020-2100




#####################################
## Calculate climate change deltas ##
#####################################
setwd('/Users/mpinsky/Documents/Rutgers/Range projections/')
agency = c('CNRM-CERFACS', 'IPSL', 'IPSL', 'MOHC', 'MPI-M', 'MPI-M', 'MRI', 'NCAR', 'NCC', 'NCC', 'NOAA-GFDL', 'NOAA-GFDL', 'NOAA-GFDL')
model = c('CNRM-CM5', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR', 'HadGem2-CC', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MRI-CGCM3', 'CCSM4', 'NorESM1-M', 'NorESM1-ME', 'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M')

basedates = read.csv('Output/basedates_2014-11-22.csv')
require(chron) # for converting julian days
nyrs = function(x) as.numeric(as.character(years(x))) # return year from a chron object

# Read in data (historical, rcp85 2020-2100, and control)
histl = vector('list', length(agency)) # "l" at end so that it matches tempscontroll and temps2100l
	for(i in 1:length(agency)){
		print(i)
		infile = paste('Output/tempshist_', agency[i], '_', model[i], '.RData', sep='')
		load(infile)
		histl[[i]] = basetemps
	}
	names(histl) = paste(agency, model, sep='_')
		
tempscontroll = vector('list', length(agency)) # "l" at end so that it doesn't get overwritten by load()
	for(i in 1:length(agency)){
		print(i)
		infile = paste('Output/tempscontrol_', agency[i], '_', model[i], '.RData', sep='')
		load(infile)
		tempscontroll[[i]] = tempscontrol
	}
	names(tempscontroll) = paste(agency, model, sep='_')

temps2100l = vector('list', length(agency)) # "l" at end so that it doesn't get overwritten by load()
	for(i in 1:length(agency)){
		print(i)
		infile = paste('Output/tempsRCP85_2020-2100_', agency[i], '_', model[i], '.RData', sep='')
		load(infile)
		temps2100 = lapply(temps2100, FUN=aperm, perm=c(1,2,3,5,4)) # transpose so that year is last dimension
		temps2100l[[i]] = temps2100
	}
	names(temps2100l) = paste(agency, model, sep='_')

regs = names(temps2100l[[1]])	
	
# Calculate deltas: offsets from historical periods
	# RCP8.5 deltas for each season in each year in each region in each model
	delta.raw2100 = vector('list', length(agency)) # each item is a model
	for(i in 1:length(agency)){ # for each model
		delta.raw2100[[i]] = vector('list', length(histl[[i]])) # each item is a region (in a model)
		for(j in 1:length(histl[[i]])){ # for each region
			a = histl[[i]][[j]] # dims are lon, lat, depth, season
			dim(a) = c(dim(a), 1) # add a fake year dimension
			# make sure lon/lat/depth dimensions of histl and temps2100 match
			a = a[,,,,rep(1,dim(temps2100l[[i]][[j]])[5])] # expand histl to match dimensions of temps2100l
			delta.raw2100[[i]][[j]] = temps2100l[[i]][[j]] - a
		}
		names(delta.raw2100[[i]]) = names(temps2100l[[i]])
	}
	names(delta.raw2100) = names(temps2100l)
	rm(a) # remove our copy of histl

	# write out
	save(delta.raw2100, file=paste('Output/delta.raw2100_', Sys.Date(), '.Rdata', sep=''))


		# Plot raw deltas for each model in each region (averaged across all depths and lon/lats)
		require(RColorBrewer)
		cols = c('#000000', brewer.pal(12, 'Set3')) 

		quartz(width=11, height=8.5)
		# pdf(width=11, height=8.5, file=paste('Figures/delta.raw_', Sys.Date(), '.pdf', sep=''))
		par(mfrow=c(4,3), mai=c(0.5, 0.5, 0.3, 0.1), mgp = c(2.3,0.8,0))
		ylims = c(-1,7.5)
		xlims=c(2020,2100)
		for(i in 1:length(delta.raw2100[[1]])){ # for each region
			print(regs[i])
			plot(1,1, col='white', xlim=xlims, ylim=ylims, xlab='Year', ylab='Delta (°C)', main=names(delta.raw2100[[1]])[i])
			for(j in 1:length(delta.raw2100)){ # for each model
				y = apply(delta.raw2100[[j]][[i]], MARGIN=5, FUN=mean, na.rm=TRUE) # mean by year within region across all lon/lat/depth
				x = (2020:2100)[1:length(y)]
				lines(x, y, col=cols[j]) # for 2020-2100
			}
		}
		legend('topleft', col=cols, legend=model, cex=0.8, lty=1, ncol=2, bty='n')
	
		dev.off()

		# Plot control deltas for each model in each region
		require(RColorBrewer)
		cols = c('#000000', brewer.pal(12, 'Set3')) 
		ymat = sapply(tempscontroll, FUN=function(x){ sapply(x, FUN=mean, na.rm=TRUE)}) # take mean by region and model

		quartz(width=6, height=4)
		# pdf(width=6, height=4, file=paste('Figures/delta.cont_', Sys.Date(), '.pdf', sep=''))
		par(mai=c(1.3, 0.7, 0.3, 0.1), mgp = c(2.3,0.8,0))
		plot(0,0, xlim=c(0,13), ylim=c(-0.017, 0.017), col='white', xlab='', ylab='Control Drift (°C/year)', xaxt='n')
		axis(1, at=1:length(regs), labels=regs, cex.axis=0.5, las=2)
		for(i in 1:length(tempscontroll)){ # for each model
			points(1:length(regs), ymat[,i], col=cols[i])
		}
		abline(h=0,col='grey')
		legend('topleft', col=cols, legend=model, cex=0.4, pch=1, ncol=2, bty='n')
	
		dev.off()

		# Compare magnitude of raw and control deltas (both in degC per year)
		compare = data.frame(agency = rep(agency, each=length(regs)), model = rep(model, each=length(regs)), region = rep(regs, times=length(model)), rawdelt = numeric(length(regs)*length(model)), controldelt = numeric(length(regs)*length(model)))
		for(i in 1:length(delta.raw2100[[1]])){ # for each region
			for(j in 1:length(delta.raw2100)){ # for each model
				y = apply(delta.raw2100[[j]][[i]], MARGIN=5, FUN=mean, na.rm=TRUE) # mean within region across years # is year in dim 1??
				x = (2020:2100)[1:length(y)]
				ind = compare$region==names(delta.raw2100[[1]])[i] & paste(as.character(compare$agency), as.character(compare$model), sep='_') == names(delta.raw2100)[j]
				compare$rawdelt[ind] = coef(lm(y ~ x))[2]
				compare$controldelt[ind] = mean(tempscontroll[[j]][[i]], na.rm=TRUE)
			}
		}
		compare$prop = compare$rawdelt/(compare$rawdelt + compare$controldelt)
		compare$ratio = compare$controldelt/compare$rawdelt

	# write out
	write.csv(compare, file=paste('Tables/compare_delta_control_', Sys.Date(), '.csv', sep=''))

# Calculate total deltas (raw delta - drift*year)
	require(plyr) # for aaply() function
	delta2100 = vector('list', length(agency))
	for(i in 1:length(agency)){ # for each model
		print(agency[i])
		delta2100[[i]] = vector('list', length(histl[[i]]))
		for(j in 1:length(histl[[i]])){ # for each region
			print(regs[j])
			print(as.character(basedates$region[j]))
			baseyear = mean(c(nyrs(basedates$min[j]), nyrs(basedates$max[j])))
			a = tempscontroll[[i]][[j]]
			dim(a) = c(dim(a), 1) # add a year dimension
			a = a[,,,,rep(1,dim(temps2100l[[i]][[j]])[5])] # expand a to match dimensions of temps2100l	(80 years)
			a = aaply(a, .margins=c(1,2,3,4), .fun= function(x){return(x*((2020:2099)-baseyear))}) # multiply degC/yr by years to get degC of drift
			delta2100[[i]][[j]] = delta.raw2100[[i]][[j]] - a
		}
		names(delta2100[[i]]) = names(delta.raw2100[[i]])
	}
	names(delta2100) = names(delta.raw2100)

	# write out
	save(delta2100, file=paste('Output/delta2100_', Sys.Date(), '.Rdata', sep=''))

# Plot deltas for each model in each region
	require(RColorBrewer)
	cols = c('#000000', brewer.pal(12, 'Set3')) 

	quartz(width=11, height=8.5)
	# pdf(width=11, height=8.5, file=paste('Figures/deltas_', Sys.Date(), '.pdf', sep=''))
	par(mfrow=c(4,3), mai=c(0.5, 0.5, 0.3, 0.1), mgp = c(2.3,0.8,0))
	ylims = c(-1,7.5)
	xlims=c(2020,2100)
	for(i in 1:length(delta2100[[1]])){ # for each region
		plot(1,1, col='white', xlim=xlims, ylim=ylims, xlab='Year', ylab='Delta (°C)', main=names(delta2100[[1]])[i])
		for(j in 1:length(delta2100)){ # for each model
			y = apply(delta2100[[j]][[i]], MARGIN=5, FUN=mean, na.rm=TRUE)
			x = as.numeric(names(y))
			lines(x, y, col=cols[j]) # for 2020-2100
		}
	}
	legend('topleft', col=cols, legend=model, cex=0.8, lty=1, ncol=2, bty='n')

	dev.off()


###################################	
# convert deltas to long format   #
# (for each region of each model) #
# takes ~ 1 hour                  #
###################################
setwd('/Users/mpinsky/Documents/Rutgers/Range projections/')
load('Output/delta2100_2014-11-22.RData')
require(reshape2) # for melt()

delta2100long = vector('list', length(delta2100)) # one for each model
names(delta2100long) =names(delta2100)
an = numeric(0) # a placeholder
ac = character(0)
for(i in 1:length(delta2100)){ # for each climate model
	print(i)
	aa = data.frame(region = ac, year = an, season = an, lat = an, lon=an, depth=an, delta = an) # holds all regions for this model 2100
	for(j in 1:length(delta2100[[i]])){ # for each region
		cat(paste(j, ' ', sep=''))
		b = melt(delta2100[[i]][[j]]) # very fast
		b$region = names(delta2100[[i]])[j]
		names(b)[names(b) == 'yr'] = 'year'
		names(b)[names(b) == 'value'] = 'delta'
		aa = rbind(aa,b)	
	}


	delta2100long[[i]] = aa
	cat('\n')
}	


# but the same lat/lon/depth can appear in multiple regions. trim those duplicates out!
for(i in 1:length(delta2100long)){
	print(nrow(delta2100long[[i]]))
	rem = duplicated(delta2100long[[i]][,c('lat', 'lon', 'depth', 'year', 'delta')])
	delta2100long[[i]] = delta2100long[[i]][!rem, c('lat', 'lon', 'depth', 'year', 'delta')]
	print(nrow(delta2100long[[i]]))

}

# write out
save(delta2100long, file=paste('Output/delta2100long_', Sys.Date(), '.RData', sep=''))
