## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	projfolder = '../CEmodels_proj' # holds model projections (outside Git)
	modfolder = '../CEModels' # holds the models (outside Git)
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	projfolder = 'CEmodels_proj'
	modfolder = 'CEmodels'
	}
# could add code for Lauren's working directory here



##########################################################
## Calculate mean lat/depth and sum of biomass by year  ##
##########################################################
require(mgcv)
require(Hmisc)

## choose which run to use
runtype <- 'test'

# load all projections from this run
files <- list.files(path = projfolder, pattern=paste('summproj_', runtype, '_', sep=''))

# weighted mean function to use with summarize()
wmean <- function(x){ # values in col 1, weights in col 2
	inds <- !is.na(x[,1]) & !is.na(x[,2])
	return(weighted.mean(x[inds,1], x[inds,2])) 
}

# save observed and predicted positions for lat, lon, and depth
meanpos <- list(0)


# set up dataframes
# Don't know how many regions for each taxon, so can't pre-allocate the rows (but could guess high... might speed this up)
# could also calc mean depth, but would need to pull from rugosity file
biomasssum = data.frame(sppocean = character(0), region = character(0), year = numeric(0)) # sum of average wtcpue across the region (2020-2099) for each survey in each model (columns)
	for(i in 1:13) biomasssum[[paste('summwtcpue', i, sep='_')]] <- numeric(0)
meanlat = data.frame(sppocean = character(0), region = character(0), year = numeric(0)) # biomass-weighted mean latitude across the region (2020-2099) for each survey in each model (columns)
	for(i in 1:13) meanlat[[paste('lat', i, sep='_')]] <- numeric(0)
meanlon = data.frame(sppocean = character(0), region = character(0), year = numeric(0))
	for(i in 1:13) meanlon[[paste('lon', i, sep='_')]] <- numeric(0)


options(warn=1) # print warnings as they occur

# loop through all files
for(i in 1:length(files)){ # takes a while (a couple hours ?)
	# load data for this species
	load(paste(projfolder, '/', files[i], sep='')) # load summproj for this taxon
	myregions <- sort(unique(summproj$region))
	mysppocean <- gsub('.Rdata', '', gsub(paste('summproj_', runtype, '_', sep=''), '', files[i]))

	print(paste(i, 'of', length(files), mysppocean, myregions, Sys.time()))

	# set up dataframes for this taxon
	mybiomasssum = data.frame(sppocean = mysppocean, region = rep(myregions, rep(80,length(myregions))), year = rep(2020:2099, length(myregions)))
	mymeanlat = data.frame(sppocean = mysppocean, region = rep(myregions, rep(80,length(myregions))), year = rep(2020:2099, length(myregions)))
	mymeanlon = data.frame(sppocean = mysppocean, region = rep(myregions, rep(80,length(myregions))), year = rep(2020:2099, length(myregions)))

	# Summarize projections
	for(j in 1:13){
		snm <- paste('wtcpue.proj', j, sep='_')
		temp <- aggregate(summproj[,snm], by=list(year = summproj$year, region = summproj$region), FUN=sum, na.rm=TRUE) # remove NAs: assumes that NA gridcells are constant through time, which they should be, since determined by the GCM grid
			names(temp)[3] <- paste('summwtcpue', j, sep='_')
		mybiomasssum <- merge(mybiomasssum, temp)

		temp <- summarize(summproj[,c('lat', snm)], by=list(region = summproj$region, year = summproj$year), FUN=wmean)
			names(temp)[3] <- paste('lat', j, sep='_')
		mymeanlat <- merge(mymeanlat, temp)

		temp <- summarize(summproj[,c('lon', snm)], by=list(region = summproj$region, year = summproj$year), FUN=wmean)
			names(temp)[3] <- paste('lon', j, sep='_')
		mymeanlon <- merge(mymeanlon, temp)

	}

	biomasssum <- rbind(biomasssum, mybiomasssum) # an inefficient way to do this: better to pre-allocate
	meanlat <- rbind(meanlat, mymeanlat)
	meanlon <- rbind(meanlon, mymeanlon)

}


### Save meanpos and biomasssum
save(biomasssum, meanlat, meanlon, file = paste('data/meanlat,lon,biomass_', runtype, '.RData', sep=''))




###########################################################
## Plot time-series of summarized spp biomass and ranges ##
###########################################################
require(Hmisc)

wmean <- function(x){
	i <- !is.na(x[,1]) & !is.na(x[,2])
	if(sum(i) == 0) return(NA)
	else return(weighted.mean(x[i,1], x[i,2])) # weighted mean function for summarize()
}

## load the data
runtype <- 'test'
load(paste('data/meanlat,lon,biomass_', runtype, '.RData', sep='')) # biomasssum, meanlat, meanlon: projected biomass by year for each taxon in each region

## plots of change in biomass
	sppregions <- sort(unique(paste(biomasssum$region, biomasssum$sppocean)))
	length(sppregions)

	# quartz(width=10, height=8)
	pdf(file=paste('figures/biomasssum_proj_', runtype, '.pdf', sep=''), width=10, height=8)
	par(mfrow = c(6,6), mai=c(0.3, 0.3, 0.2, 0.05), cex.main=0.7, cex.axis=0.8, omi=c(0,0.2,0.1,0), mgp=c(2.8, 0.7, 0), font.main=3)

	rc = 1 # row counter
	cc = 0 # column counter

	options(warn=1) # print warnings as they occur
	for(i in 1:length(sppregions)){
		if(i %% 100 == 0) print(i)
		inds <- paste(biomasssum$region, biomasssum$sppocean) == sppregions[i]
		thisreg <- unique(biomasssum$region[inds])
		thisspp <- unique(biomasssum$sppocean[inds])

		# increment row and column counters as needed
		cc = cc+1
		if(cc == 7){ cc = 1; rc = rc + 1}
		if(rc == 7){ cc = 1; rc = 1}
	
		if(i>1){ if(thisreg != oldreg){  # switch to a new page when I get to a new region
				par(mfrow = c(6,6), mai=c(0.3, 0.3, 0.2, 0.05), cex.main=0.7, cex.axis=0.8, omi=c(0,0.2,0.1,0), mgp=c(2.8, 0.7, 0), font.main=3)
				rc = 1; cc = 1
		}}
		
		# ylims
		ylims <- c(0, max(biomasssum[inds, grep('summwtcpue', names(biomasssum))], na.rm=TRUE)) # will warn if all values are NA
		if(is.infinite(ylims[2])) ylims = c(0,1)

		# plot data for each GCM
		plot(biomasssum$year[inds], biomasssum$summwtcpue_1[inds], col='grey', las=1, type='l', ylim=ylims, main=thisspp)
		for(j in 2:13) lines(biomasssum$year[inds], biomasssum[[paste('summwtcpue', j, sep='_')]][inds], col='grey')
		
		# plot ensemble mean
		agg <- cbind(biomasssum$year[inds], rowMeans(biomasssum[inds, grep('summwtcpue', names(biomasssum))]))
		lines(agg[,1], agg[,2], col='black')
	
		if(cc==1) mtext(text='Biomass index', side=2, line=2.3, cex=0.6) # add y label on left of each row
		if(rc==1) mtext(text=thisreg, side=3, line=1.3, cex=0.6) # add region header on top of page

		oldreg = thisreg # save the previous region to see if we need a new page on the next round
	}

	dev.off()


## plots of change in latitude
	sppregions <- sort(unique(paste(meanlat$region, meanlat$sppocean)))
	length(sppregions)

	# quartz(width=10, height=8)
	pdf(file=paste('figures/meanlat_proj_', runtype, '.pdf', sep=''), width=10, height=8)
	par(mfrow = c(6,6), mai=c(0.3, 0.3, 0.2, 0.05), cex.main=0.7, cex.axis=0.8, omi=c(0,0.2,0.1,0), mgp=c(2.8, 0.7, 0), font.main=3)

	rc = 1 # row counter
	cc = 0 # column counter

	options(warn=1) # print warnings as they occur
	for(i in 1:length(sppregions)){
		if(i %% 100 == 0) print(i)
		inds <- paste(meanlat$region, meanlat$sppocean) == sppregions[i]
		thisreg <- unique(meanlat$region[inds])
		thisspp <- unique(meanlat$sppocean[inds])

		# increment row and column counters as needed
		cc = cc+1
		if(cc == 7){ cc = 1; rc = rc + 1}
		if(rc == 7){ cc = 1; rc = 1}
	
		if(i>1){ if(thisreg != oldreg){  # switch to a new page when I get to a new region
				par(mfrow = c(6,6), mai=c(0.3, 0.3, 0.2, 0.05), cex.main=0.7, cex.axis=0.8, omi=c(0,0.2,0.1,0), mgp=c(2.8, 0.7, 0), font.main=3)
				rc = 1; cc = 1
		}}
		
		# ylims
		ylims <- range(meanlat[inds, grep('lat', names(meanlat))], na.rm=TRUE) # set ylims for this species in this region
#		ylims <- range(meanlat[meanlat$region == thisreg, grep('lat', names(meanlat))])	# set ylims for the whole region
		if(is.infinite(ylims[1])) ylims = c(0,1)
	
		# plot data for each GCM
		plot(meanlat$year[inds], meanlat$lat_1[inds], col='grey', las=1, type='l', ylim=ylims, main=thisspp)
		for(j in 2:13) lines(meanlat$year[inds], meanlat[[paste('lat', j, sep='_')]][inds], col='grey')
		
		# plot ensemble mean
		agg <- cbind(meanlat$year[inds], rowMeans(meanlat[inds, grep('lat', names(meanlat))]))
		lines(agg[,1], agg[,2], col='black')
	
		if(cc==1) mtext(text='Mean latitude (°N)', side=2, line=2.3, cex=0.6) # add y label on left of each row
		if(rc==1) mtext(text=thisreg, side=3, line=1.3, cex=0.6) # add region header on top of page

		oldreg = thisreg # save the previous region to see if we need a new page on the next round
	}

	dev.off()



#######################################
## Summarize distributions by decade
#######################################

# need to write


##########################################################
## Plot maps of spp projections from each climate model by decade
##########################################################
## NOTE: THIS IS OLD CODE. NEEDS TO BE RE-WRITTEN
## NOTE: not set up for annual data yet
setwd('/Users/mpinsky/Documents/Rutgers/Range projections/')
require(lattice)
require(gridExtra)
load('Output/delta2060long_2013-09-06.RData') # for climate model names
dataclim = read.csv('Output/dataclim_2013-10-23.csv', row.names=1) # to get observed distributoins
	spp = sort(unique(dataclim$taxon))
	regs = sort(unique(dataclim$region))

folderdate = '2013-11-08'
folder = paste('Output/proj_annual_', folderdate, sep='')
freeaxes = 'freeaxes_' # set this to '' to use common z-scale for all plots on a page. set to 'freeaxes_' to allow z-scale to vary for each plot
	
files=list.files(folder, pattern='proj_')
m1 = c('2020-2060', '', '', '', '', '', '', '', '', '', '', '', '') # graph titles for left column
m2 = c('2060-2100', '', '', '', '', '', '', '', '', '', '', '', '') # right column

options(warn=1) # print warnings as they occur
for(s in 1:length(files)){ # NOTE: takes 30 min or so
	print(paste(s, 'of', length(files)))
	thisspp = spp[sapply(spp, FUN=grepl, x=files[s], fixed=TRUE)] # extract spp name from file name
	thisreg = regs[sapply(regs, FUN=grepl, x=files[s])]
	vars = c('surftemp', 'bottemp')
	if(thisreg == 'DFO_Newfoundland_Fall') vars = 'bottemp'

	if(length(thisspp)>1){ # if there were multiple matches, take the one that was longest
		i = which.max(nchar(as.character(thisspp)))
		print(paste('picked', thisspp[i], 'from', paste(thisspp, collapse=', '), 'for', files[s]))
		thisspp = thisspp[i]
	}

	thisproj = read.csv(paste(folder, '/', files[s], sep=''))
	
	#average thisproj across climate models and by time period
	observ = dataclim[dataclim$region==thisreg & dataclim$taxon == thisspp,]
	histor = aggregate(list(wtcpue = thisproj$wtcpue.hist), by=list(lat=thisproj$lat, lon=thisproj$lon, depth=thisproj$depth), FUN=mean, na.rm=TRUE) # historical time period
	projave = rowMeans(thisproj[,grepl('proj_', names(thisproj))], na.rm=TRUE) # average across climate models
	inds = thisproj$year %in% (2020:2060)
	fut2060 = aggregate(list(wtcpue = projave[inds]), by=list(lat=thisproj$lat[inds], lon=thisproj$lon[inds], depth=thisproj$depth[inds]), FUN=mean, na.rm=TRUE) # 2020-2060
	inds = thisproj$year %in% (2060:2100)
	fut2100 = aggregate(list(wtcpue = projave[inds]), by=list(lat=thisproj$lat[inds], lon=thisproj$lon[inds], depth=thisproj$depth[inds]), FUN=mean, na.rm=TRUE) # 2020-2060

	# Make a set of plots
	# quartz(width=10, height=30)
	pdf(width =8, height=8, file=paste(folder, '/Figures/biomass_proj_', freeaxes, thisreg, '_', thisspp, '_', Sys.Date(), '.pdf', sep=''))
	rng = range(c(observ$biomass.clim, histor$wtcpue, fut2060$wtcpue, fut2100$wtcpue), na.rm=TRUE)
	cols = colorRampPalette(c('blue', 'purple', 'red1'), interpolate='linear')
	plots = vector('list', 4) # two time points for future, plus observed and climatology

	if(freeaxes==''){
		plots[[1]] = levelplot(biomass.clim ~ lon*lat, data=observ, at=seq(rng[1], rng[2], length.out=20), colorkey=FALSE, col.regions=cols(100), main='Observed', ylab='lat', xlab='', scales=list(cex=0.5)) # observed averaged biomass
		plots[[2]] = levelplot(wtcpue ~ lon*lat, data=histor, at=seq(rng[1], rng[2], length.out=20), col.regions=cols(100), main='Hindcast', ylab='', xlab='', scales=list(cex=0.5), colorkey=list(axis.text=list(cex=0.5))) # predicted biomass based on climatology
		plots[[3]] = levelplot(wtcpue~lon*lat, data=fut2060, at=seq(rng[1], rng[2], length.out=20), colorkey=FALSE, col.regions=cols(100), main='2020-2060', ylab='lat', scales=list(cex=0.5)) # predicted to 2020-2060
		plots[[4]] = levelplot(wtcpue~lon*lat, data=fut2100, at=seq(rng[1], rng[2], length.out=20), col.regions=cols(100), main='2060-2100', ylab='', scales=list(cex=0.5), colorkey=FALSE) # predicted 2060-2100
	}
	if(freeaxes=='freeaxes_'){
		plots[[1]] = levelplot(biomass.clim ~ lon*lat, data=observ, colorkey=list(axis.text=list(cex=0.5)), col.regions=cols(100), main='Observed', ylab='lat', xlab='', scales=list(cex=0.5)) # observed averaged biomass
		plots[[2]] = levelplot(wtcpue ~ lon*lat, data=histor, colorkey=list(axis.text=list(cex=0.5)), col.regions=cols(100), main='Hindcast', ylab='', xlab='', scales=list(cex=0.5)) # predicted biomass based on climatology
		plots[[3]] = levelplot(wtcpue~lon*lat, data=fut2060, colorkey=list(axis.text=list(cex=0.5)), col.regions=cols(100), main='2020-2060', ylab='lat', scales=list(cex=0.5)) # predicted to 2020-2060
		plots[[4]] = levelplot(wtcpue~lon*lat, data=fut2100, colorkey=list(axis.text=list(cex=0.5)), col.regions=cols(100), main='2060-2100', ylab='', scales=list(cex=0.5)) # predicted 2060-2100
	}
	args = paste(paste(paste('plots[[', 1:4, ']]', sep=''), collapse = ', '), ', ncol=2', sep='') # the arguments
	eval(parse(text=paste('grid.arrange(', args, ')'))) # plot all on one page

	dev.off()
	
}

