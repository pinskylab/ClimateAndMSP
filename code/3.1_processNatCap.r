# Summarize the NatCap InVest output onto our analysis grid

##############################
## Process NatCap Wind data ##
##############################
require(raster)
require(sf)
require(data.table)

gridsize = 0.05 # size of grid of the climate data, in degrees


# read in data
wind_west <- raster('../NatCap_temp/westcoastwind/output/npv_US_millions.tif')
	# image(wind_west)
wind_east <- raster('../NatCap_temp/eastcoastwind/output/npv_US_millions.tif')
	# image(wind_east)

grid <- readRDS('temp/SPsf2.rds') # the analysis grid

# project to LL
wind_east.t <- projectRaster(wind_east, crs=crs(grid))
wind_west.t <- projectRaster(wind_west, crs=crs(grid))

# extract from raster to the grid cells: VERY SLOW approach
# wind_east_df <- extract(x=wind_east.t, y=as(grid[grid$lon > -100,], 'Spatial'), method='bilinear', fun=mean, na.rm=TRUE) # get raster values by climate grid cell
# wind_west_df <- extract(x=wind_west.t, y=as(grid[grid$lon < -100,], 'Spatial'), method='bilinear', fun=mean, na.rm=TRUE)
# wind_east_df <- cbind(npv = wind_east_df, grid[grid$lon > -100, c('lon', 'lat')])
# wind_west_df <- cbind(npv = wind_west_df, grid[grid$lon < -100, c('lon', 'lat')])

# extract from raster to the grid cells: fast approach
wind_east_dt <- data.table(cbind(coordinates(wind_east.t), npv=extract(x=wind_east.t, y=extent(wind_east.t)))) # get raster values by raster grid cell
wind_west_dt <- data.table(cbind(coordinates(wind_west.t), npv=extract(x=wind_west.t, y=extent(wind_west.t))))
wind_dt <- rbind(wind_east_dt, wind_west_dt) # concatenate

wind_dt[, latgrid := floor(y/gridsize)*gridsize + gridsize/2] # round to nearest climate grid center
wind_dt[, longrid := floor(x/gridsize)*gridsize + gridsize/2]

wind_sum <- wind_dt[, .(npv = mean(npv, na.rm = TRUE)), by = c('latgrid', 'longrid')] # average by climate grid cell

# plot to make sure it worked
wind_sum[, plot(longrid, latgrid, col=c('red', 'blue')[1+(npv > 0)], pch=16, cex=0.05)] # plots of <> a threshold


# mark which grids are in climate grid
grid$latgrid <- floor(grid$lat/gridsize)*gridsize + gridsize/2 # round to nearest climate grid center (to fix some rounding errors)
grid$longrid <- floor(grid$lon/gridsize)*gridsize + gridsize/2

wind_sum[, keep := FALSE] # set up a column to mark the ones to keep
wind_sum[paste(latgrid, longrid) %in% paste(grid$latgrid, grid$longrid), keep := TRUE] # keep if in the climate grid
    wind_sum[, sum(keep)]
    wind_sum[, sum(!keep)]

    wind_sum[keep == TRUE, plot(longrid, latgrid, col=c('red', 'blue')[1+(npv > 0)], pch=16, cex=0.05)] # plots of <> a threshold
    
# remove grids not in clim
wind.out <- wind_sum[keep == TRUE, .(lat = latgrid, lon = longrid, npv = npv)]

# convert NAs to lowest value (too deep)
minnpv <- wind.out[!is.na(npv), min(npv)]
wind.out[is.na(npv), npv := minnpv]

    wind.out[, plot(lon, lat, col=c('red', 'blue')[1+(npv > 0)], pch=16, cex=0.05)] # plots of <> a threshold

# are all climate grid cells in the wind object?
missing <- !(paste(grid$latgrid, grid$longrid) %in% wind.out[, paste(lat, lon)])
sum(missing) # 0 
    
# write out
saveRDS(wind.out, 'output/wind_npv.rds')


##############################
## Process NatCap Wave data ##
##############################
require(raster)
require(sf)
require(data.table)

gridsize = 0.05 # size of grid of the climate data, in degrees

# read in data
wave_west <- raster('../NatCap_temp/westcoastwave/output/npv_usd.tif')
    # image(wave_west)
wave_east <- raster('../NatCap_temp/eastcoastwave/output/npv_usd.tif')
    # image(wave_east)

grid <- readRDS('temp/SPsf2.rds') # the analysis grid

# project to LL
wave_east.t <- projectRaster(wave_east, crs=crs(grid))
wave_west.t <- projectRaster(wave_west, crs=crs(grid)) # took 15 min. why? because has to grid the whole globe from -180 to 180.

# split west into east and west of -180
# otherwise R doesn't have enough memory to do the next step all at once
wave_west.t1 <- crop(wave_west.t, extent(165, 180, 40, 65))
wave_west.t2 <- crop(wave_west.t, extent(-180, -100, 0, 90))

# extract from raster to the grid cells: fast approach
wave_east_dt <- data.table(cbind(coordinates(wave_east.t), npv=extract(x=wave_east.t, y=extent(wave_east.t)))) # get raster values by raster grid cell
wave_west_dt1 <- data.table(cbind(coordinates(wave_west.t1), npv=extract(x=wave_west.t1, y=extent(wave_west.t1))))
wave_west_dt2 <- data.table(cbind(coordinates(wave_west.t2), npv=extract(x=wave_west.t2, y=extent(wave_west.t2))))
wave_dt <- rbind(wave_east_dt, wave_west_dt1, wave_west_dt2) # concatenate
    nrow(wave_east_dt)
    nrow(wave_west_dt1)
    nrow(wave_west_dt2) # very big
    nrow(wave_dt) # very big

wave_dt[, latgrid := floor(y/gridsize)*gridsize + gridsize/2] # round to nearest climate grid center
wave_dt[, longrid := floor(x/gridsize)*gridsize + gridsize/2]

wave_sum <- wave_dt[, .(npv = mean(npv, na.rm = TRUE)), by = c('latgrid', 'longrid')] # average by climate grid cell
    nrow(wave_sum) # more reasonable

# plot to make sure it worked
wave_sum[, plot(longrid, latgrid, col=c('red', 'blue')[1+(npv > 0)], pch=16, cex=0.05)] # plots of <> a threshold


# mark which grids are in climate grid
grid$latgrid <- floor(grid$lat/gridsize)*gridsize + gridsize/2 # round to nearest climate grid center (to fix some rounding errors)
grid$longrid <- floor(grid$lon/gridsize)*gridsize + gridsize/2

wave_sum[, keep := FALSE] # set up a column to mark the ones to keep
wave_sum[paste(latgrid, longrid) %in% paste(grid$latgrid, grid$longrid), keep := TRUE] # keep if in the climate grid
wave_sum[, sum(keep)]
wave_sum[, sum(!keep)]

wave_sum[keep == TRUE, plot(longrid, latgrid, col=c('red', 'blue')[1+(npv > 0)], pch=16, cex=0.05)] # plots of <> a threshold

# remove grids not in clim
wave.out <- wave_sum[keep == TRUE, .(lat = latgrid, lon = longrid, npv = npv)]

# convert NAs to lowest value (too deep)
minnpv <- wave.out[!is.na(npv), min(npv)]
wave.out[is.na(npv), npv := minnpv]

wave.out[, plot(lon, lat, col=c('red', 'blue')[1+(npv > 0)], pch=16, cex=0.05)] # plots of <> a threshold

# are all climate grid cells in the wave object?
missing <- !(paste(grid$latgrid, grid$longrid) %in% wave.out[, paste(lat, lon)])
sum(missing) # 0 

# write out
saveRDS(wave.out, 'output/wave_npv.rds')

