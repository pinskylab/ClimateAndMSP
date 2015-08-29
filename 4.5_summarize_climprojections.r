# basic summary statistics and plots of the BT and SST projections

# Read in temperature fields and models, then make range projections
# This could probably be sped up by switching from data.frames to data.tables

require(mgcv)
require(Hmisc)
require(parallel) # for multi-core calculations


## Set working directory
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	}
# could add code for Lauren's working directory here


########################
## Prep projection data
########################
load('data/climGrid.proj2.RData') # loads clim

# reshape clim to long format
require(reshape2)
names(clim)[names(clim)=='bottemp.clim.int'] = 'bottemp.proj_0' # "model" 0 is the climatology
names(clim)[names(clim)=='surftemp.clim.int'] = 'surftemp.proj_0'

clim2 = melt(clim[, !(names(clim) %in% c(paste('depthgrid', 1:13, sep=''), 'bottemp.clim', 'surftemp.clim', 'latgrid', 'longrid'))], id = c('region', 'season', 'lat', 'lon', 'depth', 'year')) # fast (a few sec)! drop bottemp.clim and surftemp.clim.
	rm(clim)
	newcols = strsplit(as.character(clim2$variable), split='_') # a few minutes
	# newcols = colsplit(clim2$variable, names = c('variable', 'model'), pattern='_') # 10 min+. takes too long
	newcols2 = unlist(newcols)
	var1 = newcols2[seq(1,length(newcols2), by=2)]
	var2 = newcols2[seq(2,length(newcols2), by=2)]
	clim2$variable = var1 # slow
	clim2$model = as.numeric(var2) # slow

	clim2 = clim2[!(clim2$model == 0 & clim2$year > 2006),] # remove climatology entries for all years after 2006. they are just duplicates and not needed in this long format

	save(clim2, file='data/climGrid.proj2long_10deginterp.RData') # slow

	rm(newcols, newcols2, var1, var2)

###########################################
## Prep delta data
## Run this on Amphiprion, takes 60G RAM
###########################################
require(reshape2)
load('data/delta2100long.RData') # slow (1.3G)

	# reshape into one long dataframe
delta2100long[[1]]$model = 1
delta2100long2 = delta2100long[[1]]

for(i in 2:length(delta2100long)){ # takes 10 min or so on Amphiprion
	print(paste(i, Sys.time()))
	delta2100long[[i]]$model = i
	delta2100long2 = rbind(delta2100long2, delta2100long[[i]])
}

save(delta2100long2, file='data/delta2100long2.RData') # slow


###########################################
## summarize raw delta data by time period
###########################################
# uses 80GB of RAM on Amphiprion
# load delta2100long2.RData if needed. don't reload if already in memory, since very slow

delta2100long3 = melt(delta2100long2, id = c('lon', 'lat', 'depth', 'season', 'year', 'region', 'model')) # melted version for manipulation
	dim(delta2100long3) # 309,792,340 rows

# identify the surface temperatures
delta2100long3$surf = FALSE
mods = sort(unique(delta2100long3$model))
for(i in 1:length(mods)){ # for each model, since each model has a different depth grid
	print(i)
	k = delta2100long3$model == i # subset to this model
	mindep = min(delta2100long3$depth[k])
	delta2100long3$surf[k & delta2100long3$depth == mindep] = TRUE
}

# define three time periods
delta2100long3$period <- NA 
delta2100long3$period[delta2100long3$year <= 2020] <- 1
delta2100long3$period[delta2100long3$year > 2020 & delta2100long3$year <= 2060] <- 2
delta2100long3$period[delta2100long3$year > 2060] <- 3

# take mean SST by time period
deltaSSTperiods = dcast(data=delta2100long3[delta2100long3$surf == TRUE,], region + lat + lon + season + period + model ~ variable, fun.aggregate = mean, na.rm=TRUE) # average SSTdeltas across all years in each period (before/after 2020 and 2060)
	dim(deltaSSTperiods) # 354,588
	head(deltaSSTperiods)
	sort(unique(deltaSSTperiods$model))
	sort(unique(deltaSSTperiods$period))
	deltaSSTperiods$period[deltaSSTperiods$period == 1] = '2006-2020'
	deltaSSTperiods$period[deltaSSTperiods$period == 2] = '2021-2060'
	deltaSSTperiods$period[deltaSSTperiods$period == 3] = '2061-2100'

save(deltaSSTperiods, file='data/deltaSSTperiods.RData') # fast

## DO THE SAME FOR BT?


###################################
## Plots                       ####
###################################
load('data/climGrid.proj2long_10deginterp.RData') # loads clim2 (long format)
load('data/deltaSSTperiods.RData') # average SST delta by time period for each grid cell. using all data from the climate models.

# Plot maps of raw (un-interpolated) deltas in each region/season, across three time periods. Each region on a separate page

		#SST
	require(lattice)
	source('packet.panel.bycolumn.R') # so that I can plot across multiple pages
	cols = colorRampPalette(colors = c('blue', 'white', 'red'))
	pdf(width=30, height=6, file=paste('figures/deltaSST_nointerp.pdf', sep=''))
	regs = sort(unique(deltaSSTperiods$region))
	seasonnames = c('winter', 'spring', 'summer', 'fall')
	for(i in 1:length(regs)) {
		print(i)
		seasons = sort(unique(deltaSSTperiods$season[deltaSSTperiods$region == regs[i]]))
		for(j in seasons){
			p <- levelplot(delta ~ lon*lat|model*period, scales = list(x='free', y='free'), data=deltaSSTperiods[deltaSSTperiods$region == regs[i] & deltaSSTperiods$season == j,], col.regions=cols, par.strip.text=list(cex=0.5), strip=function(...,strip.levels) strip.default(..., strip.levels=TRUE), main=paste(regs[i], seasonnames[j]), panel = function(...){
				panel.fill(col='light grey')
				panel.levelplot(...)
			})
			plot(p, packet.panel=packet.panel.bycolumn) # need to do it this way in order to plot lattice across multiple pages
		}
	}

	dev.off()


# Plot maps of all model projections in a region/season, across three time periods. Each region on a separate page.
	# uses clim2
	# first summarize projections by time period
	clim2$period <- NA 
	clim2$period[clim2$year <= 2020] <- 1
	clim2$period[clim2$year > 2020 & clim2$year <= 2060] <- 2
	clim2$period[clim2$year > 2060] <- 3

	clim2$variable[clim2$variable == 'bottemp.proj'] = 'bottemp' # now the projections are distinguished from the climatology by model (>0 vs. 0)
	clim2$variable[clim2$variable == 'surftemp.proj'] = 'surftemp'	
	clim2periods = dcast(data=clim2, region + lat + lon + season + period + model ~ variable, fun.aggregate = mean, na.rm=TRUE) # average across all years in each period (before/after 2020 and 2060)
		head(clim2periods)
		clim2periods$period[clim2periods$period == 1] = '2006-2020'
		clim2periods$period[clim2periods$period == 2] = '2021-2060'
		clim2periods$period[clim2periods$period == 3] = '2061-2100'


		#BT
	require(lattice)
	source('packet.panel.bycolumn.R') # so that I can plot across multiple pages
	cols = colorRampPalette(colors = c('blue', 'white', 'red'))
	pdf(width=30, height=6, file=paste('figures/climBTproj_10deginterp.pdf', sep=''))
	regs = sort(unique(clim2periods$region))
	seasonnames = c('winter', 'spring', 'summer', 'fall')
	for(i in 1:length(regs)) {
		print(i)
		seasons = sort(unique(clim2periods$season[clim2periods$region == regs[i]]))
		for(j in seasons){
			p <- levelplot(bottemp ~ lon*lat|model*period, scales = list(x='free', y='free'), data=clim2periods[clim2periods$region == regs[i] & clim2periods$season == j,], col.regions=cols, par.strip.text=list(cex=0.5), strip=function(...,strip.levels) strip.default(..., strip.levels=TRUE), main=paste(regs[i], seasonnames[j]), panel = function(...){
				panel.fill(col='light grey')
				panel.levelplot(...)
			})
			plot(p, packet.panel=packet.panel.bycolumn) # need to do it this way in order to plot lattice across multiple pages
		}
	}

	dev.off()

		#SST
	require(lattice)
	source('packet.panel.bycolumn.R') # so that I can plot across multiple pages
	cols = colorRampPalette(colors = c('blue', 'white', 'red'))
	pdf(width=30, height=6, file=paste('figures/climSSTproj_10deginterp.pdf', sep=''))
	regs = sort(unique(clim2periods$region))
	seasonnames = c('winter', 'spring', 'summer', 'fall')
	for(i in 1:length(regs)) {
		print(i)
		seasons = sort(unique(clim2periods$season[clim2periods$region == regs[i]]))
		for(j in seasons){
			p <- levelplot(surftemp ~ lon*lat|model*period, scales = list(x='free', y='free'), data=clim2periods[clim2periods$region == regs[i] & clim2periods$season == j,], col.regions=cols, par.strip.text=list(cex=0.5), strip=function(...,strip.levels) strip.default(..., strip.levels=TRUE), main=paste(regs[i], seasonnames[j]), panel = function(...){
				panel.fill(col='light grey')
				panel.levelplot(...)
			})
			plot(p, packet.panel=packet.panel.bycolumn) # need to do it this way in order to plot lattice across multiple pages
		}
	}

	dev.off()



# Plot biplot of historical and future temperatures by region and model (annual data)
## OLD CODE. WON'T WORK UNTIL UPDATED
#	data = read.csv('../../Princeton/Trawl Data/Output/goodhauls_allregions_2012-11-19.csv', row.names=1)
#	data$lon[data$lon < 0] = data$lon[data$lon < 0] + 360 # fix lons to only positive to match CMIP5 data
#	temp = data[complete.cases(data[,c('lat', 'lon', 'depth')]) & data$region %in% c('AFSC_EBS', 'AFSC_Aleutians', 'AFSC_GOA', 'WestCoast_Tri', 'SEFSC_GOMex', 'NEFSC_Spring', 'DFO_ScotianShelf', 'DFO_SoGulf', 'DFO_Newfoundland_Fall') ,c('region', 'stratum', 'lat', 'lon', 'yearsurv', 'depth', 'bottemp', 'surftemp')]
#		temp = droplevels(temp)	
#	temp = temp[temp$year<2006,] # Trim to no later than December 16, 2005 (to match historical period in CMIP5)
#	gridsize=0.25 # size of grid of the climate data, in degrees
#	temp$latgrid = floor(temp$lat/gridsize)*gridsize + gridsize/2 # round to nearest grid center
#	temp$longrid = floor(temp$lon/gridsize)*gridsize + gridsize/2
#
#	regs = sort(unique(clim$region))
#	quartz(height=9, width=9)
#	pdf(height = 9, width = 9, file=paste('Figures/nonanalog.pdf', sep=''))
#
#	for(i in 1:13){ # for each climate model
#		print(i)
#		#bmp(height = 9, width = 9, res=300, units='in', file=paste('Figures/nonanalog_', i, '_', Sys.Date(), '.bmp', sep=''))
#		par(mfrow=c(3,3), oma=c(0,0,1.5,0), mai= c(0.5, 0.5, 0.3, 0.1), mgp=c(2,1,0))
#		for(j in 1:length(regs)){
#			inds = clim$region == regs[j]
#			inds1 = temp$region == regs[j]
#			xlims = range(c(temp$bottemp[inds1], clim[[paste('bottemp.proj_', i, sep='')]][inds]), na.rm=TRUE)
#			
#			if(regs[j] == 'DFO_Newfoundland_Fall'){
#				inds2 = clim$year %in% (2061:2100)
#				plot(clim[[paste('bottemp.proj_', i, sep='')]][inds & inds2], jitter(rep(1, sum(inds&inds2))), col='red', xlab='Bottom temperature (°C)', ylab='not an axis', main=regs[j], xlim=xlims)
#				inds2 = clim$year %in% (2020:2060)
#				points(clim[[paste('bottemp.proj_', i, sep='')]][inds & inds2], jitter(rep(1, sum(inds&inds2))), col='blue')
#				points(temp$bottemp[inds1],  jitter(rep(1, sum(inds1))))				
#			} else {
#				ylims = range(c(temp$surftemp[inds1], clim[[paste('surftemp.proj_', i, sep='')]][inds]), na.rm=TRUE)
#			
#				inds2 = clim$year %in% (2061:2100)
#				plot(clim[[paste('bottemp.proj_', i, sep='')]][inds&inds2], clim[[paste('surftemp.proj_', i, sep='')]][inds&inds2], col='red', xlab='Bottom temperature (°C)', ylab='Surface temperature (°C)', main=regs[j], xlim=xlims, ylim=ylims)
#				inds2 = clim$year %in% (2020:2060)
#				points(clim[[paste('bottemp.proj_', i, sep='')]][inds&inds2], clim[[paste('surftemp.proj_', i, sep='')]][inds&inds2], col='blue')
#				points(temp$bottemp[inds1], temp$surftemp[inds1])
#			}
#		}
#		mtext(names(delta2060long)[i], outer=TRUE, side=3)
#		#dev.off()
#	}
#	
#	dev.off()
#
#
#
## Plot maps of where nonanalog future temperatures will occur, by region and model (using climatologies) 
## OLD CODE. WON'T WORK UNTIL UPDATED
#	regs = sort(unique(clim$region))
#	pch =16
#	cex = 0.7
#
#	pdf(height = 9, width = 9, file=paste('figures/nonanalog_maps.pdf', sep=''))
#
#	for(i in 1:13){ # for each climate model
#		par(mfrow=c(3,3), oma=c(0,0,1.5,0), mai= c(0.5, 0.5, 0.3, 0.1), mgp=c(2,1,0))
#		for(j in 1:length(regs)){
#			inds = clim$region == regs[j]
#			inds1 = temp$region == regs[j]
#			maxb = max(temp$bottemp[inds1], na.rm=TRUE)
#
#			if(regs[j] == 'DFO_Newfoundland_Fall'){
#				nonb = clim[[paste('bottemp.2100_', i, sep='')]][inds] > maxb # BT exceeds historical max
#				plot(clim$lon[inds], clim$lat[inds], col='grey', xlab='Lon', ylab='Lat', main=regs[j], pch=pch, cex=cex)
#				points(clim$lon[inds][nonb], clim$lat[inds][nonb], col='blue', pch=pch, cex=cex)
#			} else {
#				maxs = max(temp$surftemp[inds1], na.rm=TRUE)
#				nonb = ((clim[[paste('bottemp.proj_', i, sep='')]][inds] > maxb) & !(clim[[paste('surftemp.proj_', i, sep='')]][inds] > maxs)) # BT but not SST exceeds historical max
#				nons = (!(clim[[paste('bottemp.proj_', i, sep='')]][inds] > maxb) & (clim[[paste('surftemp.proj_', i, sep='')]][inds] > maxs)) # SST but not BT exceeds historical max
#				non2 = 2 == ((clim[[paste('bottemp.proj_', i, sep='')]][inds] > maxb) + (clim[[paste('surftemp.proj_', i, sep='')]][inds] > maxs)) # both of BT and SST exceeds historical max
#			
#				plot(clim$lon[inds], clim$lat[inds], col='grey', xlab='Lon', ylab='Lat', main=regs[j], pch=pch, cex=cex)
#				points(clim$lon[inds][nonb], clim$lat[inds][nonb], col='blue', pch=pch, cex=cex)
#				points(clim$lon[inds][nons], clim$lat[inds][nons], col='green', pch=pch, cex=cex)
#				points(clim$lon[inds][non2], clim$lat[inds][non2], col='red', pch=pch, cex=cex)
#			}
#		}
#		mtext(names(delta2060long)[i], outer=TRUE, side=3)
#	}
#	
#	dev.off()
#
## Plot maps of where nonanalog future temperatures will occur by region (for BT and for SST, summarized across models)
## OLD CODE. WON'T WORK UNTIL UPDATED
#	regs = sort(unique(clim$region))
#	pch =16
#	cex = 0.7
#	cols = c(gray(0.5), rev(heat.colors(13)))
#
#	pdf(height = 9, width = 9, file=paste('figures/nonanalog_mapsummary.pdf', sep=''))
#
#	# for BT
#	par(mfrow=c(3,3), oma=c(0,0,1.5,0), mai= c(0.5, 0.5, 0.3, 0.1), mgp=c(2,1,0))
#	for(j in 1:length(regs)){
#		inds = temp$region == regs[j]
#		maxb = max(temp$bottemp[inds], na.rm=TRUE)
#
#		inds1 = clim$region == regs[j]	
#		ct = rep(0, length(unique(paste(clim$lat[inds1], clim$lon[inds1])))) # counts how many models exceed maxb at each lat/lon
#		for(i in 1:13){ # for each model
#			modb = aggregate(list(bottemp = clim[[paste('bottemp.proj_', i, sep='')]][inds1]), by=list(lat = clim$lat[inds1], lon = clim$lon[inds1]), FUN=max, na.rm=TRUE) # find max temp at each lat/lon in the model projections
#			modb = modb[order(modb$lat, modb$lon),]
#			ct = ct + (modb$bottemp > maxb)
#		}
#			
#		plot(modb$lon, modb$lat, col=cols[ct+1], xlab='Lon', ylab='Lat', main=regs[j], pch=pch, cex=cex)
#	}
#	mtext("Bottom Temperature 2060-2100", outer=TRUE, side=3)
#	legend('bottomleft', col=cols, pch=16, legend=0:13, cex=0.9)
#
#	# for SST
#	par(mfrow=c(3,3), oma=c(0,0,1.5,0), mai= c(0.5, 0.5, 0.3, 0.1), mgp=c(2,1,0))
#	for(j in 1:length(regs)){
#		inds = temp$region == regs[j]
#		inds1 = clim$region == regs[j]
#		if(regs[j] == 'DFO_Newfoundland_Fall'){
#			plot(clim$lon[inds1], clim$lat[inds1], col='white', xlab='Lon', ylab='Lat', main=regs[j], pch=pch, cex=cex)
#		} else {
#			maxb = max(temp$surftemp[inds], na.rm=TRUE)
#	
#			ct = rep(0, length(unique(paste(clim$lat[inds1], clim$lon[inds1])))) # counts how many models exceed maxb at each lat/lon
#			for(i in 1:13){ # for each model
#				modb = aggregate(list(surftemp = clim[[paste('surftemp.proj_', i, sep='')]][inds1]), by=list(lat = clim$lat[inds1], lon = clim$lon[inds1]), FUN=max, na.rm=TRUE) # find max temp at each lat/lon in the model projections
#				modb = modb[order(modb$lat, modb$lon),]
#				ct = ct + (modb$surftemp > maxb)
#			}
#			
#			plot(modb$lon, modb$lat, col=cols[ct+1], xlab='Lon', ylab='Lat', main=regs[j], pch=pch, cex=cex)
#		}
#	}
#	mtext("Surface Temperature 2060-2100", outer=TRUE, side=3)
#	legend('bottomleft', col=cols, pch=16, legend=0:13, cex=0.9)
#
#	
#	dev.off()
#
