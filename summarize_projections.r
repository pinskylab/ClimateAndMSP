# basic summary stats and plots of the BT and SST projections

#######################
## Prep data
#######################
setwd('/Users/mpinsky/Documents/Rutgers/Range projections/')
# read back in
	load('Output/climGrid.proj2_2015-02-10.RData')
	# clim = read.csv('Output/climGrid.proj2_2013-11-08.csv', row.names=1); type='Grid'

# reshape clim to long format
require(reshape2)
#clim2 = reshape(clim, direction='long', varying = list(bottemp.delta = grep('bottemp.delta', names(clim), value=TRUE), bottemp.proj = grep('bottemp.proj', names(clim), value=TRUE), surftemp.delta = grep('surftemp.delta', names(clim), value=TRUE), surftemp.proj = grep('surftemp.proj', names(clim), value=TRUE)), timevar = 'model', times = 1:13, idvar = c('region', 'latgrid', 'longrid', 'season', 'lat', 'lon', 'bottemp.clim', 'surftemp.clim', 'depth')) # failed after using 5GB of memory
names(clim)[names(clim)=='bottemp.clim.int'] = 'bottemp.proj_0' # "model" 0 is the climatology
names(clim)[names(clim)=='surftemp.clim.int'] = 'surftemp.proj_0'

clim2 = melt(clim[, !(names(clim) %in% c(paste('depthgrid', 1:13, sep=''), 'bottemp.clim', 'surftemp.clim', 'latgrid', 'longrid'))], id = c('region', 'season', 'lat', 'lon', 'depth', 'year')) # fast (15 sec?). drop bottemp.clim and surftemp.clim.
	rm(clim)
	newcols = strsplit(as.character(clim2$variable), split='_') # a few minutes
	# newcols = colsplit(clim2$variable, names = c('variable', 'model'), pattern='_') # 10 min+. takes too long
	newcols2 = unlist(newcols)
	var1 = newcols2[seq(1,length(newcols2), by=2)]
	var2 = newcols2[seq(2,length(newcols2), by=2)]
	clim2$variable = var1 # slow
	clim2$model = as.numeric(var2) # slow

	clim2 = clim2[!(clim2$model == 0 & clim2$year > 2020),] # remove climatology entries for all years after 2020. they are just duplicates

	save(clim2, file=paste('Output/climGrid.proj2long_10deginterp_', Sys.Date(), '.RData', sep='')) # slow

	rm(newcols, newcols2, var1, var2)

########
## Prep delta data
## Run this on Amphiprion, takes 60G RAM
########
require(reshape2)
load('Output/delta2100long_2014-12-11.RData') # slow (1.3G)

	# reshape into one long dataframe
delta2100long[[1]]$model = 1
delta2100long2 = delta2100long[[1]]

for(i in 2:length(delta2100long)){
	print(i)
	delta2100long[[i]]$model = i
	delta2100long2 = rbind(delta2100long2, delta2100long[[i]])
}

save(delta2100long2, file=paste('Output/delta2100long2_', Sys.Date(), '.RData', sep='')) # slow

##############################
## summarize raw delta data
##############################

delta2100long3 = melt(delta2100long2, id = c('lon', 'lat', 'depth', 'season', 'year', 'region', 'model')) # melted version for manipulation

	delta2100long3$surf = FALSE
	mods = sort(unique(delta2100long3$model))
	for(i in 1:length(mods)){ # for each model, since each model has a different depth grid
		print(i)
		k = delta2100long3$model == i # subset to this model
		mindep = min(delta2100long3$depth[k])
		delta2100long3$surf[k & delta2100long3$depth == mindep] = TRUE
	}

	delta2100long3$period = delta2100long3$year > 2060
	deltaSSTperiods = dcast(data=delta2100long3[delta2100long3$surf == TRUE,], region + lat + lon + season + period + model ~ variable, fun.aggregate = mean) # average SSTdeltas across all years in each period (before/after 2060)
		dim(deltaSSTperiods)
		head(deltaSSTperiods)
		sort(unique(deltaSSTperiods$model))
		deltaSSTperiods$period[deltaSSTperiods$period == FALSE] = '2020-2060'
		deltaSSTperiods$period[deltaSSTperiods$period == 'TRUE'] = '2060-2100'

save(deltaSSTperiods, file=paste('Output/deltaSSTperiods_', Sys.Date(), '.RData', sep='')) # fast


#######################
## summarize from clim (projected data)
#######################

# summarize average deltas and future temperatures by model and region
	delbt20202060 = array(NA, dim=c(13, length(unique(clim$region)), 4), dimnames=list(model = paste('Mod', 1:13, sep=''), region = sort(unique(clim$region)), season=1:4)) # deltas of BT. rows are regions, columns are models. levels are season. For 2020-2060.
	delst20202060 = array(NA, dim=c(13, length(unique(clim$region)), 4), dimnames=list(model = paste('Mod', 1:13, sep=''), region = sort(unique(clim$region)), season=1:4)) # deltas of BT. rows are regions, columns are models. levels are season
	delbt20602100 = array(NA, dim=c(13, length(unique(clim$region)), 4), dimnames=list(model = paste('Mod', 1:13, sep=''), region = sort(unique(clim$region)), season=1:4)) # deltas of BT. rows are regions, columns are models. levels are season. For 2060-2100.
	delst20602100 = array(NA, dim=c(13, length(unique(clim$region)), 4), dimnames=list(model = paste('Mod', 1:13, sep=''), region = sort(unique(clim$region)), season=1:4)) # deltas of BT. rows are regions, columns are models. levels are season
	scaff = data.frame(region = sort(unique(clim$region))) # scaffold on which to add regional means
	for(i in 1:13){ # for each model
		print(i)
		for(j in 1:4){ # for each season
			dnm = paste('bottemp.delta_', i, sep='')
			inds = clim$year %in% (2020:2060) & clim$season == j
			temp = merge(scaff, aggregate(clim[[dnm]][inds], by=list(region=clim$region[inds]), FUN=mean, na.rm=TRUE), all.x=TRUE)
			delbt20202060[i,,j] = temp$x
			dnm = paste('surftemp.delta_', i, sep='')
			temp = merge(scaff, aggregate(clim[[dnm]][inds], by=list(region=clim$region[inds]), FUN=mean, na.rm=TRUE), all.x=TRUE)
			delst20202060[i,,j]  = temp$x

			inds = clim$year %in% (2060:2100) & clim$season == j
			temp = merge(scaff, aggregate(clim[[dnm]][inds], by=list(region=clim$region[inds]), FUN=mean, na.rm=TRUE), all.x=TRUE)
			delbt20602100[i,,j] = temp$x
			dnm = paste('surftemp.delta_', i, sep='')
			temp = merge(scaff, aggregate(clim[[dnm]][inds], by=list(region=clim$region[inds]), FUN=mean, na.rm=TRUE), all.x=TRUE)
			delst20602100[i,,j]  = temp$x

		}
	}
	delbt20202060 = melt(delbt20202060)
		names(delbt20202060)[names(delbt20202060)=='value'] = 'delta'
	delst20202060 = melt(delst20202060)
		names(delst20202060)[names(delst20202060)=='value'] = 'delta'
	delbt20602100 = melt(delbt20602100)
		names(delbt20602100)[names(delbt20602100)=='value'] = 'delta'
	delst20602100 = melt(delst20602100)
		names(delst20602100)[names(delst20602100)=='value'] = 'delta'

	write.csv(delbt20202060, file=paste('Tables/delta_BT2020_2060_', Sys.Date(), '.csv', sep=''))
	write.csv(delst20202060, file=paste('Tables/delta_SST2020_2060_', Sys.Date(), '.csv', sep=''))
	write.csv(delbt20602100, file=paste('Tables/delta_BT2060_2100_', Sys.Date(), '.csv', sep=''))
	write.csv(delst20602100, file=paste('Tables/delta_SST2060_2100_', Sys.Date(), '.csv', sep=''))

	# future temps (NOT YET UPDATED)
	bt0 = matrix(NA, ncol=length(delta2060long), nrow = length(unique(clim$region)))
		rownames(delbt60) = sort(unique(clim$region)); colnames(delbt60) = names(delta2060long)
	bt60 = matrix(NA, ncol=length(delta2060long), nrow = length(unique(clim$region)))
		rownames(delbt60) = sort(unique(clim$region)); colnames(delbt60) = names(delta2060long)
	bt100 = matrix(NA, ncol=length(delta2060long), nrow = length(unique(clim$region)))
		rownames(delbt100) = sort(unique(clim$region)); colnames(delbt100) = names(delta2060long)
	st0 = matrix(NA, ncol=length(delta2060long), nrow = length(unique(clim$region)))
		rownames(delst60) = sort(unique(clim$region)); colnames(delst60) = names(delta2060long)
	st60 = matrix(NA, ncol=length(delta2060long), nrow = length(unique(clim$region)))
		rownames(delst60) = sort(unique(clim$region)); colnames(delst60) = names(delta2060long)
	st100 = matrix(NA, ncol=length(delta2060long), nrow = length(unique(clim$region)))
		rownames(delst100) = sort(unique(clim$region)); colnames(delst100) = names(delta2060long)
	for(i in 1:length(delta2060long)){ # for each model
		dnm = 'bottemp.clim.int'
		bt0[,i] = aggregate(clim[[dnm]], by=list(region=clim$region), FUN=mean)[,2]
		dnm = paste('bottemp.proj_', i, sep='')
		inds = clim$year %in% (2020:2060)
		bt60[,i] = aggregate(clim[[dnm]][inds], by=list(region=clim$region[inds]), FUN=mean)[,2]
		inds = clim$year %in% (2060:2100)
		bt100[,i] = aggregate(clim[[dnm]][inds], by=list(region=clim$region[inds]), FUN=mean)[,2]
		dnm = 'surftemp.clim.int'
		st0[,i] = aggregate(clim[[dnm]][inds], by=list(region=clim$region[inds]), FUN=mean)[,2]
		dnm = paste('surftemp.proj_', i, sep='')
		inds = clim$year %in% (2020:2060)
		st60[,i] = aggregate(clim[[dnm]][inds], by=list(region=clim$region[inds]), FUN=mean)[,2]
		inds = clim$year %in% (2060:2100)
		st100[,i] = aggregate(clim[[dnm]][inds], by=list(region=clim$region[inds]), FUN=mean)[,2]
	}
	bt0
	bt60
	bt100
	st0
	st60
	st100

	write.csv(bt0, file=paste('Tables/BTclim_', Sys.Date(), '.csv', sep=''))
	write.csv(bt60, file=paste('Tables/BT2060_', Sys.Date(), '.csv', sep=''))
	write.csv(bt100, file=paste('Tables/BT2100_', Sys.Date(), '.csv', sep=''))
	write.csv(st0, file=paste('Tables/SSTclim_', Sys.Date(), '.csv', sep=''))
	write.csv(st60, file=paste('Tables/SST2060_', Sys.Date(), '.csv', sep=''))
	write.csv(st100, file=paste('Tables/SST2100_', Sys.Date(), '.csv', sep=''))


###################################
## Plots                       ####
###################################
load('Output/climGrid.proj2long_10deginterp_2015-02-10.RData') # loads clim2 (long format)
load('Output/deltaSSTperiods_2015-02-07.RData') # average SST delta by time period for each grid cell. using all data from the climate models.

# Plot maps of deltas in each region/season, across three time periods. Each region on a separate page

		#BT
	require(lattice)
	source('packet.panel.bycolumn.R') # so that I can plot across multiple pages
	cols = colorRampPalette(colors = c('blue', 'white', 'red'))
	pdf(width=30, height=6, file=paste('Figures/deltaSST_nointerp_', Sys.Date(), '.pdf', sep=''))
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
	clim2$period = clim2$year > 2060
	clim2$variable[clim2$variable == 'bottemp.proj'] = 'bottemp' # now the projections are distinguished from the climatology by model (>0 vs. 0)
	clim2$variable[clim2$variable == 'surftemp.proj'] = 'surftemp'	
	clim2periods = dcast(data=clim2, region + lat + lon + season + period + model ~ variable, fun.aggregate = mean) # average across all years in each period (before/after 2060)
		head(clim2periods)
		clim2periods$period[clim2periods$period == FALSE] = '2020-2060'
		clim2periods$period[clim2periods$period == 'TRUE'] = '2060-2100'


		#BT
	require(lattice)
	source('packet.panel.bycolumn.R') # so that I can plot across multiple pages
	cols = colorRampPalette(colors = c('blue', 'white', 'red'))
	pdf(width=30, height=6, file=paste('Figures/climBTproj_10deginterp_', Sys.Date(), '.pdf', sep=''))
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
	pdf(width=30, height=6, file=paste('Figures/climSSTproj_10deginterp_', Sys.Date(), '.pdf', sep=''))
	regs = sort(unique(clim2periods$region))
		regs = regs[!grepl('Newfoundland|WCAnn', regs)]
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



# Plot maps of each region across three time periods on a separate page
	require(gridExtra)

	# BT
	quartz(width=10, height=30)
	# pdf(width = 10, height=13*3, file=paste('Figures/climBTproj_grid_interp_', Sys.Date(), '.pdf', sep=''))
	regs = sort(unique(clim$region))
	for(i in 1:length(regs)){
		print(i)
		inds = clim$region == regs[i]
		rng = range(c(clim[inds,'bottemp.clim.int'], clim[inds,paste('bottemp.2060_', 1:13, sep='')], clim[inds,paste('bottemp.2100_', 1:13, sep='')]), na.rm=TRUE)
		plots = vector('list', 13*3)
		for(mod in 1:13){ # for each model
			bt1 = paste('bottemp.2060_', mod, sep='') # column names in clim
			bt2 = paste('bottemp.2100_', mod, sep='')
			cols = colorRampPalette(c('purple', 'blue', 'red1'), interpolate='linear')
			plots[[(mod-1)*3+1]] = levelplot(bottemp.clim.int ~ lon*lat, data=clim[inds,], at=seq(rng[1], rng[2], length.out=20), colorkey=FALSE, col.regions=cols(100), main='Climatology', ylab=names(delta2060long)[mod])
			plots[[(mod-1)*3+2]] = levelplot(formula(paste(bt1, '~lon*lat')), data=clim[inds,], at=seq(rng[1], rng[2], length.out=20), colorkey=FALSE, col.regions=cols(100), main='2060', ylab='')
			plots[[(mod-1)*3+3]] = levelplot(formula(paste(bt2, '~lon*lat')), data=clim[inds,], at=seq(rng[1], rng[2], length.out=20), col.regions=cols(100), main='2100', ylab='')
			
		}
		args = paste(paste(paste('plots[[', 1:39, ']]', sep=''), collapse = ', '), ', ncol=3', sep='') # the arguments
		fnc = 'grid.arrange'
		eval(parse(text=paste(fnc, "(", args, ")"))) # plot all on one page
	}

	dev.off()


	# SST
	quartz(width=10, height=30)
	# pdf(width = 10, height=13*3, file=paste('Figures/climSSTproj_grid_interp_', Sys.Date(), '.pdf', sep=''))
	regs = sort(unique(clim$region))
		regs = setdiff(regs, 'DFO_Newfoundland_Fall')
	for(i in 1:length(regs)){
		print(i)
		inds = clim$region == regs[i]
		rng = range(c(clim[inds,'surftemp.clim.int'], clim[inds,paste('surftemp.2060_', 1:13, sep='')], clim[inds,paste('surftemp.2100_', 1:13, sep='')]), na.rm=TRUE)
		plots = vector('list', 13*3)
		for(mod in 1:13){ # for each model
			bt1 = paste('surftemp.2060_', mod, sep='') # column names in clim
			bt2 = paste('surftemp.2100_', mod, sep='')
			cols = colorRampPalette(c('purple', 'blue', 'red1'), interpolate='linear')
			plots[[(mod-1)*3+1]] = levelplot(surftemp.clim.int ~ lon*lat, data=clim[inds,], at=seq(rng[1], rng[2], length.out=20), colorkey=FALSE, col.regions=cols(100), main='Climatology', ylab=names(delta2060long)[mod])
			plots[[(mod-1)*3+2]] = levelplot(formula(paste(bt1, '~lon*lat')), data=clim[inds,], at=seq(rng[1], rng[2], length.out=20), colorkey=FALSE, col.regions=cols(100), main='2060', ylab='')
			plots[[(mod-1)*3+3]] = levelplot(formula(paste(bt2, '~lon*lat')), data=clim[inds,], at=seq(rng[1], rng[2], length.out=20), col.regions=cols(100), main='2100', ylab='')
			
		}
		args = paste(paste(paste('plots[[', 1:39, ']]', sep=''), collapse = ', '), ', ncol=3', sep='') # the arguments
		fnc = 'grid.arrange'
		eval(parse(text=paste(fnc, "(", args, ")"))) # plot all on one page
	}

	dev.off()

# Plot biplot of historical and future temperatures by region and model (annual data)
	data = read.csv('../../Princeton/Trawl Data/Output/goodhauls_allregions_2012-11-19.csv', row.names=1)
	data$lon[data$lon < 0] = data$lon[data$lon < 0] + 360 # fix lons to only positive to match CMIP5 data
	temp = data[complete.cases(data[,c('lat', 'lon', 'depth')]) & data$region %in% c('AFSC_EBS', 'AFSC_Aleutians', 'AFSC_GOA', 'WestCoast_Tri', 'SEFSC_GOMex', 'NEFSC_Spring', 'DFO_ScotianShelf', 'DFO_SoGulf', 'DFO_Newfoundland_Fall') ,c('region', 'stratum', 'lat', 'lon', 'yearsurv', 'depth', 'bottemp', 'surftemp')]
		temp = droplevels(temp)	
	temp = temp[temp$year<2006,] # Trim to no later than December 16, 2005 (to match historical period in CMIP5)
	gridsize=0.25 # size of grid of the climate data, in degrees
	temp$latgrid = floor(temp$lat/gridsize)*gridsize + gridsize/2 # round to nearest grid center
	temp$longrid = floor(temp$lon/gridsize)*gridsize + gridsize/2

	regs = sort(unique(clim$region))
	quartz(height=9, width=9)
	pdf(height = 9, width = 9, file=paste('Figures/nonanalog_', Sys.Date(), '.pdf', sep=''))

	for(i in 1:13){ # for each climate model
		print(i)
		#bmp(height = 9, width = 9, res=300, units='in', file=paste('Figures/nonanalog_', i, '_', Sys.Date(), '.bmp', sep=''))
		par(mfrow=c(3,3), oma=c(0,0,1.5,0), mai= c(0.5, 0.5, 0.3, 0.1), mgp=c(2,1,0))
		for(j in 1:length(regs)){
			inds = clim$region == regs[j]
			inds1 = temp$region == regs[j]
			xlims = range(c(temp$bottemp[inds1], clim[[paste('bottemp.proj_', i, sep='')]][inds]), na.rm=TRUE)
			
			if(regs[j] == 'DFO_Newfoundland_Fall'){
				inds2 = clim$year %in% (2061:2100)
				plot(clim[[paste('bottemp.proj_', i, sep='')]][inds & inds2], jitter(rep(1, sum(inds&inds2))), col='red', xlab='Bottom temperature (°C)', ylab='not an axis', main=regs[j], xlim=xlims)
				inds2 = clim$year %in% (2020:2060)
				points(clim[[paste('bottemp.proj_', i, sep='')]][inds & inds2], jitter(rep(1, sum(inds&inds2))), col='blue')
				points(temp$bottemp[inds1],  jitter(rep(1, sum(inds1))))				
			} else {
				ylims = range(c(temp$surftemp[inds1], clim[[paste('surftemp.proj_', i, sep='')]][inds]), na.rm=TRUE)
			
				inds2 = clim$year %in% (2061:2100)
				plot(clim[[paste('bottemp.proj_', i, sep='')]][inds&inds2], clim[[paste('surftemp.proj_', i, sep='')]][inds&inds2], col='red', xlab='Bottom temperature (°C)', ylab='Surface temperature (°C)', main=regs[j], xlim=xlims, ylim=ylims)
				inds2 = clim$year %in% (2020:2060)
				points(clim[[paste('bottemp.proj_', i, sep='')]][inds&inds2], clim[[paste('surftemp.proj_', i, sep='')]][inds&inds2], col='blue')
				points(temp$bottemp[inds1], temp$surftemp[inds1])
			}
		}
		mtext(names(delta2060long)[i], outer=TRUE, side=3)
		#dev.off()
	}
	
	dev.off()



# Plot maps of where nonanalog future temperatures will occur, by region and model (using climatologies) NOT YET CONVERTED TO ANNUAL DATA
	regs = sort(unique(clim$region))
	pch =16
	cex = 0.7

	pdf(height = 9, width = 9, file=paste('Figures/nonanalog_maps_', Sys.Date(), '.pdf', sep=''))

	for(i in 1:13){ # for each climate model
		par(mfrow=c(3,3), oma=c(0,0,1.5,0), mai= c(0.5, 0.5, 0.3, 0.1), mgp=c(2,1,0))
		for(j in 1:length(regs)){
			inds = clim$region == regs[j]
			inds1 = temp$region == regs[j]
			maxb = max(temp$bottemp[inds1], na.rm=TRUE)

			if(regs[j] == 'DFO_Newfoundland_Fall'){
				nonb = clim[[paste('bottemp.2100_', i, sep='')]][inds] > maxb # BT exceeds historical max
				plot(clim$lon[inds], clim$lat[inds], col='grey', xlab='Lon', ylab='Lat', main=regs[j], pch=pch, cex=cex)
				points(clim$lon[inds][nonb], clim$lat[inds][nonb], col='blue', pch=pch, cex=cex)
			} else {
				maxs = max(temp$surftemp[inds1], na.rm=TRUE)
				nonb = ((clim[[paste('bottemp.proj_', i, sep='')]][inds] > maxb) & !(clim[[paste('surftemp.proj_', i, sep='')]][inds] > maxs)) # BT but not SST exceeds historical max
				nons = (!(clim[[paste('bottemp.proj_', i, sep='')]][inds] > maxb) & (clim[[paste('surftemp.proj_', i, sep='')]][inds] > maxs)) # SST but not BT exceeds historical max
				non2 = 2 == ((clim[[paste('bottemp.proj_', i, sep='')]][inds] > maxb) + (clim[[paste('surftemp.proj_', i, sep='')]][inds] > maxs)) # both of BT and SST exceeds historical max
			
				plot(clim$lon[inds], clim$lat[inds], col='grey', xlab='Lon', ylab='Lat', main=regs[j], pch=pch, cex=cex)
				points(clim$lon[inds][nonb], clim$lat[inds][nonb], col='blue', pch=pch, cex=cex)
				points(clim$lon[inds][nons], clim$lat[inds][nons], col='green', pch=pch, cex=cex)
				points(clim$lon[inds][non2], clim$lat[inds][non2], col='red', pch=pch, cex=cex)
			}
		}
		mtext(names(delta2060long)[i], outer=TRUE, side=3)
	}
	
	dev.off()

# Plot maps of where nonanalog future temperatures will occur by region (for BT and for SST, summarized across models)
# updated to use annual temperatures
	regs = sort(unique(clim$region))
	pch =16
	cex = 0.7
	cols = c(gray(0.5), rev(heat.colors(13)))

	pdf(height = 9, width = 9, file=paste('Figures/nonanalog_mapsummary_', Sys.Date(), '.pdf', sep=''))

	# for BT
	par(mfrow=c(3,3), oma=c(0,0,1.5,0), mai= c(0.5, 0.5, 0.3, 0.1), mgp=c(2,1,0))
	for(j in 1:length(regs)){
		inds = temp$region == regs[j]
		maxb = max(temp$bottemp[inds], na.rm=TRUE)

		inds1 = clim$region == regs[j]	
		ct = rep(0, length(unique(paste(clim$lat[inds1], clim$lon[inds1])))) # counts how many models exceed maxb at each lat/lon
		for(i in 1:13){ # for each model
			modb = aggregate(list(bottemp = clim[[paste('bottemp.proj_', i, sep='')]][inds1]), by=list(lat = clim$lat[inds1], lon = clim$lon[inds1]), FUN=max, na.rm=TRUE) # find max temp at each lat/lon in the model projections
			modb = modb[order(modb$lat, modb$lon),]
			ct = ct + (modb$bottemp > maxb)
		}
			
		plot(modb$lon, modb$lat, col=cols[ct+1], xlab='Lon', ylab='Lat', main=regs[j], pch=pch, cex=cex)
	}
	mtext("Bottom Temperature 2060-2100", outer=TRUE, side=3)
	legend('bottomleft', col=cols, pch=16, legend=0:13, cex=0.9)

	# for SST
	par(mfrow=c(3,3), oma=c(0,0,1.5,0), mai= c(0.5, 0.5, 0.3, 0.1), mgp=c(2,1,0))
	for(j in 1:length(regs)){
		inds = temp$region == regs[j]
		inds1 = clim$region == regs[j]
		if(regs[j] == 'DFO_Newfoundland_Fall'){
			plot(clim$lon[inds1], clim$lat[inds1], col='white', xlab='Lon', ylab='Lat', main=regs[j], pch=pch, cex=cex)
		} else {
			maxb = max(temp$surftemp[inds], na.rm=TRUE)
	
			ct = rep(0, length(unique(paste(clim$lat[inds1], clim$lon[inds1])))) # counts how many models exceed maxb at each lat/lon
			for(i in 1:13){ # for each model
				modb = aggregate(list(surftemp = clim[[paste('surftemp.proj_', i, sep='')]][inds1]), by=list(lat = clim$lat[inds1], lon = clim$lon[inds1]), FUN=max, na.rm=TRUE) # find max temp at each lat/lon in the model projections
				modb = modb[order(modb$lat, modb$lon),]
				ct = ct + (modb$surftemp > maxb)
			}
			
			plot(modb$lon, modb$lat, col=cols[ct+1], xlab='Lon', ylab='Lat', main=regs[j], pch=pch, cex=cex)
		}
	}
	mtext("Surface Temperature 2060-2100", outer=TRUE, side=3)
	legend('bottomleft', col=cols, pch=16, legend=0:13, cex=0.9)

	
	dev.off()

