# Calculate area and fraction of each grid cell that is in a PA

library(sf)

# Read in gridded WDPA data
wdpagrid <- readRDS('temp/wdpa_by_grid0.05_intersect.rds'); gridsz <- 0.05
SPsf2 <- readRDS('temp/SPsf2.rds') # the analysis grid
	
# Calc fraction of each grid covered by a PA
	# Merge together wdpa pieces in the same grid cell
	out2 <- aggregate(x=wdpagrid[,'gridpolyID'], by=list(ID=wdpagrid$gridpolyID), FUN=unique, na.rm=TRUE) # do the merge by grid cell ID (=rowname in SP)

	# convert to equal-area projection
	out3 <- st_transform(out2, crs='+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs') # North America Albers Equal Area Conic from https://epsg.io/102008
    SPsf3 <- st_transform(SPsf2, crs=st_crs(out3))
	
	# Calculate area of the grid intersected by PAs
	SPsf3$area_grid <- st_area(SPsf3)
    out3$area_wdpa=st_area(out3) # area of all intersected wdpa pieces
    wdpaprop <- merge(st_drop_geometry(out3[,c('gridpolyID', 'area_wdpa')]), st_drop_geometry(SPsf3[,c('gridpolyID', 'lat', 'lon', 'area_grid')]), by='gridpolyID')
    wdpaprop$prop_wdpa = as.numeric(wdpaprop$area_wdpa/wdpaprop$area_grid) # turn to a proportion of the grid cell
        with(wdpaprop, sum(area_wdpa > area_grid)) # a few are too big
        max(wdpaprop$prop_wdpa)
    wdpaprop$prop_wdpa[wdpaprop$prop_wdpa>1] <- 1 # fix the cases that are slightly >1
	
    
	# add grid cells with 0 coverage
	missinggridpolyIDs = setdiff(SPsf3$gridpolyID, wdpaprop$gridpolyID) # all the rows in out grid (SP) that didn't find a matching protected area
	new = st_drop_geometry(SPsf3[SPsf3$gridpolyID %in% missinggridpolyIDs, c('gridpolyID', 'lat', 'lon', 'area_grid')])
    new$area_wdpa <- new$prop_wdpa <- 0 # add PA area and proportion
    units(new$area_wdpa) <- units(wdpaprop$area_wdpa) # set the units
	
	wdpacov = rbind(wdpaprop, new)		
		dim(wdpacov)
		nrow(SPsf3) # should match
		hist(wdpacov$prop_wdpa, xlab='Proportion of grid cell covered by WDPA', main='') # basically all or nothing, since grid size so small
		hist(wdpacov$prop_wdpa[!(wdpacov$prop_wdpa %in% c(0,1))], xlab='Proportion of grid cell covered by WDPA', main='Not 0 or 1') # the rest
		
	# re-order
	wdpacov = wdpacov[order(wdpacov$gridpolyID),]

	# Write out
	saveRDS(wdpacov, file=paste('output/grid', gridsz, '_cov_by_wdpa.rds', sep=''))

# Calculate the fraction of each PA in each grid cell
    # get area of each WDPA piece in each grid cell
	wdpagrid2 <- st_transform(wdpagrid, crs='+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs') # North America Albers Equal Area Conic from https://epsg.io/102008
	wdpagrid2$area_wdpa=st_area(wdpagrid2) # area of all intersected wdpa pieces
	
	# get area of each WDPA PA
	wdpa = st_read('dataDL/WDPA/WDPA_Aug2019_marine-shapefile/WDPA_Aug2019_marine-shapefile-polygons.shp')
	wdpa2 <- st_transform(wdpa, crs='+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs') # North America Albers Equal Area Conic from https://epsg.io/102008
	wdpa2$area_fullwdpa <- st_area(wdpa2)
	wdpa.by.grid <- merge(st_drop_geometry(wdpagrid2), st_drop_geometry(wdpa2[,c('WDPA_PID', 'area_fullwdpa')]))
		dim(wdpa.by.grid)
		dim(wdpagrid2) # should match
		head(wdpa.by.grid)
	wdpa.by.grid$prop_grid <- as.numeric(wdpa.by.grid$area_wdpa/wdpa.by.grid$area_fullwdpa)
		summary(wdpa.by.grid$prop_grid) # max is very slightly above 1, which must be a mistake
	wdpa.by.grid$prop_grid[wdpa.by.grid$prop_grid>1] <- 1 # fix the cases that are slightly >1
		
	# re-order
	wdpa.by.grid <- wdpa.by.grid[order(wdpa.by.grid$WDPA_PID, wdpa.by.grid$gridpolyID),]
		head(wdpa.by.grid)
	
# write out
	saveRDS(wdpa.by.grid, file=paste('output/wdpa_cov_by_grid', gridsz, '.rds', sep=''))
	