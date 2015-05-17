## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	projfolder <- '../CEmodels_proj' # holds model projections (outside Git)
	modfolder <- '../CEModels' # holds the models (outside Git)
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	projfolder <- 'CEmodels_proj'
	modfolder <- 'CEmodels'
	}
# could add code for Lauren's working directory here



##########################################################
## Calculate mean lat/depth and sum of biomass by year  ##
## averages over seasons                                ##
##########################################################
require(mgcv)
require(Hmisc)

## choose which run to use
runtype <- 'testseason'

# list all projections from this run
files <- list.files(path = projfolder, pattern=paste('summproj_', runtype, '_', sep=''))

# weighted mean function to use with summarize()
wmean <- function(x){ # values in col 1, weights in col 2
	inds <- !is.na(x[,1]) & !is.na(x[,2])
	return(weighted.mean(x[inds,1], x[inds,2])) 
}

# save predicted positions for lat, lon, and depth
meanpos <- list(0)


# set up dataframes
# Don't know how many regions for each taxon, so can't pre-allocate the rows (but could guess high... might speed this up)
# could also calc mean depth, but would need to pull from rugosity file
biomasssum <- data.frame(sppocean = character(0), region = character(0), year = numeric(0)) # sum of average wtcpue across the region (2020-2099) for each survey in each model (columns)
	for(i in 1:13) biomasssum[[paste('summwtcpue', i, sep='_')]] <- numeric(0)
meanlat <- data.frame(sppocean = character(0), region = character(0), year = numeric(0)) # biomass-weighted mean latitude across the region (2020-2099) for each survey in each model (columns)
	for(i in 1:13) meanlat[[paste('lat', i, sep='_')]] <- numeric(0)
meanlon <- data.frame(sppocean = character(0), region = character(0), year = numeric(0))
	for(i in 1:13) meanlon[[paste('lon', i, sep='_')]] <- numeric(0)


options(warn=1) # print warnings as they occur

# loop through all files
for(i in 1:length(files)){ # takes a while (a couple hours ?)
	# load data for this species
	load(paste(projfolder, '/', files[i], sep='')) # load summproj for this taxon
	myregions <- sort(unique(summproj$region))
	mysppocean <- gsub('.Rdata', '', gsub(paste('summproj_', runtype, '_', sep=''), '', files[i]))

	print(paste(i, 'of', length(files), mysppocean, paste(myregions, collapse=', '), Sys.time()))

	# set up dataframes for this taxon
	mybiomasssum <- data.frame(sppocean = mysppocean, region = rep(myregions, rep(80,length(myregions))), year = rep(2020:2099, length(myregions)))
	mymeanlat <- data.frame(sppocean = mysppocean, region = rep(myregions, rep(80,length(myregions))), year = rep(2020:2099, length(myregions)))
	mymeanlon <- data.frame(sppocean = mysppocean, region = rep(myregions, rep(80,length(myregions))), year = rep(2020:2099, length(myregions)))

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
## averages over seasons                                 ##
###########################################################
require(Hmisc)

wmean <- function(x){
	i <- !is.na(x[,1]) & !is.na(x[,2])
	if(sum(i) == 0) return(NA)
	else return(weighted.mean(x[i,1], x[i,2])) # weighted mean function for summarize()
}

## load the data
runtype <- 'testseason'
load(paste('data/meanlat,lon,biomass_', runtype, '.RData', sep='')) # biomasssum, meanlat, meanlon: projected biomass by year for each taxon in each region

## plots of change in biomass
	sppregions <- sort(unique(paste(biomasssum$region, biomasssum$sppocean)))
	length(sppregions)

	# quartz(width=10, height=8)
	pdf(file=paste('figures/biomasssum_proj_', runtype, '.pdf', sep=''), width=10, height=8)
	par(mfrow = c(6,6), mai=c(0.3, 0.3, 0.2, 0.05), cex.main=0.7, cex.axis=0.8, omi=c(0,0.2,0.1,0), mgp=c(2.8, 0.7, 0), font.main=3)

	rc <- 1 # row counter
	cc <- 0 # column counter

	options(warn=1) # print warnings as they occur
	for(i in 1:length(sppregions)){
		if(i %% 100 == 0) print(i)
		inds <- paste(biomasssum$region, biomasssum$sppocean) == sppregions[i]
		thisreg <- unique(biomasssum$region[inds])
		thisspp <- unique(biomasssum$sppocean[inds])

		# increment row and column counters as needed
		cc <- cc+1
		if(cc == 7){ cc <- 1; rc <- rc + 1}
		if(rc == 7){ cc <- 1; rc <- 1}
	
		if(i>1){ if(thisreg != oldreg){  # switch to a new page when I get to a new region
				par(mfrow = c(6,6), mai=c(0.3, 0.3, 0.2, 0.05), cex.main=0.7, cex.axis=0.8, omi=c(0,0.2,0.1,0), mgp=c(2.8, 0.7, 0), font.main=3)
				rc <- 1; cc <- 1
		}}
		
		# ylims
		ylims <- c(0, max(biomasssum[inds, grep('summwtcpue', names(biomasssum))], na.rm=TRUE)) # will warn if all values are NA
		if(is.infinite(ylims[2])) ylims <- c(0,1)

		# plot data for each GCM
		plot(biomasssum$year[inds], biomasssum$summwtcpue_1[inds], col='grey', las=1, type='l', ylim=ylims, main=thisspp)
		for(j in 2:13) lines(biomasssum$year[inds], biomasssum[[paste('summwtcpue', j, sep='_')]][inds], col='grey')
		
		# plot ensemble mean
		agg <- cbind(biomasssum$year[inds], rowMeans(biomasssum[inds, grep('summwtcpue', names(biomasssum))]))
		lines(agg[,1], agg[,2], col='black')
	
		if(cc==1) mtext(text='Biomass index', side=2, line=2.3, cex=0.6) # add y label on left of each row
		if(rc==1) mtext(text=thisreg, side=3, line=1.3, cex=0.6) # add region header on top of page

		oldreg <- thisreg # save the previous region to see if we need a new page on the next round
	}

	dev.off()


## plots of change in latitude
	sppregions <- sort(unique(paste(meanlat$region, meanlat$sppocean)))
	length(sppregions)

	# quartz(width=10, height=8)
	pdf(file=paste('figures/meanlat_proj_', runtype, '.pdf', sep=''), width=10, height=8)
	par(mfrow = c(6,6), mai=c(0.3, 0.3, 0.2, 0.05), cex.main=0.7, cex.axis=0.8, omi=c(0,0.2,0.1,0), mgp=c(2.8, 0.7, 0), font.main=3)

	rc <- 1 # row counter
	cc <- 0 # column counter

	options(warn=1) # print warnings as they occur
	for(i in 1:length(sppregions)){
		if(i %% 100 == 0) print(i)
		inds <- paste(meanlat$region, meanlat$sppocean) == sppregions[i]
		thisreg <- unique(meanlat$region[inds])
		thisspp <- unique(meanlat$sppocean[inds])

		# increment row and column counters as needed
		cc <- cc+1
		if(cc == 7){ cc <- 1; rc <- rc + 1}
		if(rc == 7){ cc <- 1; rc <- 1}
	
		if(i>1){ if(thisreg != oldreg){  # switch to a new page when I get to a new region
				par(mfrow = c(6,6), mai=c(0.3, 0.3, 0.2, 0.05), cex.main=0.7, cex.axis=0.8, omi=c(0,0.2,0.1,0), mgp=c(2.8, 0.7, 0), font.main=3)
				rc <- 1; cc <- 1
		}}
		
		# ylims
		ylims <- range(meanlat[inds, grep('lat', names(meanlat))], na.rm=TRUE) # set ylims for this species in this region
#		ylims <- range(meanlat[meanlat$region == thisreg, grep('lat', names(meanlat))])	# set ylims for the whole region
		if(is.infinite(ylims[1])) ylims <- c(0,1)
	
		# plot data for each GCM
		plot(meanlat$year[inds], meanlat$lat_1[inds], col='grey', las=1, type='l', ylim=ylims, main=thisspp)
		for(j in 2:13) lines(meanlat$year[inds], meanlat[[paste('lat', j, sep='_')]][inds], col='grey')
		
		# plot ensemble mean
		agg <- cbind(meanlat$year[inds], rowMeans(meanlat[inds, grep('lat', names(meanlat))]))
		lines(agg[,1], agg[,2], col='black')
	
		if(cc==1) mtext(text='Mean latitude (°N)', side=2, line=2.3, cex=0.6) # add y label on left of each row
		if(rc==1) mtext(text=thisreg, side=3, line=1.3, cex=0.6) # add region header on top of page

		oldreg <- thisreg # save the previous region to see if we need a new page on the next round
	}

	dev.off()



###################################################
## Summarize distributions by two-decade periods
###################################################
require(mgcv)
require(Hmisc)

# weighted mean function to use with summarize()
wmean <- function(x){ # values in col 1, weights in col 2
	inds <- !is.na(x[,1]) & !is.na(x[,2])
	return(weighted.mean(x[inds,1], x[inds,2])) 
}


## choose which run and time periods to use
runtype <- 'test'
timeperiods <- data.frame(year = 2006:2100, period = c(rep('2006-2020', 15), rep('2021-2040', 20), rep('2041-2060', 20), rep('2061-2080', 20), rep('2081-2100', 20)))
periods <- sort(unique(timeperiods$period))
nt <- length(unique(periods))

# list all projections from this run
files <- list.files(path = projfolder, pattern=paste('summproj_', runtype, '_', sep=''))


# set up dataframes
# Don't know how many regions for each taxon, so can't pre-allocate the rows (but could guess high... might speed this up)
# could also calc mean depth, but would need to pull from rugosity file
biomassavemap <- data.frame(sppocean = character(0), region = character(0), period = character(0), season = numeric(0), lat = numeric(0), lon = numeric(0), wtcpue.proj = numeric(0))

options(warn=1) # print warnings as they occur

# loop through all files
for(i in 1:length(files)){ # takes a while (a couple hours ?)
	# load data for this species
	load(paste(projfolder, '/', files[i], sep='')) # load summproj for this taxon
	myregions <- sort(unique(summproj$region))
	mysppocean <- gsub('.Rdata', '', gsub(paste('summproj_', runtype, '_', sep=''), '', files[i]))

	print(paste(i, 'of', length(files), mysppocean, paste(myregions, collapse=', '), Sys.time()))

	summproj <- merge(summproj, timeperiods)

	# set up dataframe for this taxon
	inds <- !duplicated(summproj[,c('region', 'lat', 'lon', 'season')]) # unique locations/seasons to record
	mymap <- data.frame(sppocean = rep(mysppocean, sum(inds)*nt), region = rep(summproj$region[inds], nt), period = rep(sort(unique(timeperiods$period)), rep(sum(inds), nt)), season = rep(summproj$season[inds], nt), lat = rep(summproj$lat[inds], nt), lon = rep(summproj$lon[inds], nt), wtcpue.proj = NA)

	mymap <- mymap[order(mymap$region, mymap$period, mymap$season, mymap$lat, mymap$lon),]
	summproj <- summproj[order(summproj$region, summproj$year, summproj$season, summproj$lat, summproj$lon),]

	# Summarize projections across all models within time periods
	cols <- grep('wtcpue.proj', names(summproj))
	for(j in 1:nt){
		inds <- summproj$period == periods[j]
		inds2 <- mymap$period == periods[j]
		temp <- apply(summproj[inds,cols], MARGIN=1, FUN=mean) # average across models
		temp2 <- aggregate(list(wtcpue.proj = temp), by=list(region = summproj$region[inds], period = summproj$period[inds], season = summproj$season[inds], lat = summproj$lat[inds], lon = summproj$lon[inds]), FUN=mean, na.rm=TRUE)
		temp2 <- temp2[order(temp2$region, temp2$period, temp2$season, temp2$lat, temp2$lon),]

		if(all(temp2$region == mymap$region & temp2$period == mymap$period & temp2$season == mymap$season & temp2$lat == mymap$lat & temp2$lon == mymap$long)){
			mymap$wtcpue.proj[inds2] <- temp2$wtcpue.proj
		} else {
			warning(paste('rows do not match on j=', j))		
		}
	}

	biomassavemap <- rbind(biomassavemap, mymap) # an inefficient way to do this: better to pre-allocate

}

summary(biomassavemap)
dim(biomassavemap)


### Save meanpos and biomasssum
save(biomassavemap, file = paste('data/biomassavemap_', runtype, '.RData', sep=''))


################################################################
## Plot ensemble mean maps of spp projections by 20-year block
################################################################
require(lattice)
require(gridExtra)

runtype <- 'test' # which run type to use
	
load(paste('data/biomassavemap_', runtype, '.RData', sep=''))
	
sppregseas <- biomassavemap[!duplicated(biomassavemap[,c('sppocean', 'region', 'season')]), c('sppocean', 'region', 'season')] # each spp/region combination to make maps for
	sppregseas <- sppregseas[order(sppregseas$region, sppregseas$season, sppregseas$sppocean),]
	nrow(sppregseas)

# Make a set of plots on separate pages, one for each spp/region/season combination
options(warn=1) # print warnings as they occur
cols <- colorRampPalette(c('grey80', 'blue', 'purple', 'red1'), interpolate='linear')
periods <- sort(unique(biomassavemap$period))

pdf(width=10, height=3, file=paste('figures/biomass_proj_maps_', runtype, '.pdf', sep=''))

for(i in 1:nrow(sppregseas)){
	print(paste(i, 'of', nrow(sppregseas), Sys.time()))
	inds <- biomassavemap$sppocean == sppregseas$sppocean[i] & biomassavemap$region == sppregseas$region[i]
	maintitle <- paste(sppregseas$region[i], sppregseas$season[i], sppregseas$sppocean[i])

	rng <- c(0, 1.01*max(biomassavemap$wtcpue.proj[inds], na.rm=TRUE)) # slightly expanded to capture all values

	thisplot <- levelplot(wtcpue.proj ~ lon*lat|period, data=biomassavemap[inds,], at=seq(rng[1], rng[2], length.out=20), colorkey=list(axis.text=list(cex=0.5)), col.regions=cols(100), main=list(label=maintitle, cex=1), ylab=list(label='lat', cex=0.5), xlab=list(label='lon', cex=0.5), scales=list(cex=0.5), layout = c(length(periods), 1)) # observed averaged biomass

	grid.arrange(thisplot) # plot on one page
	
}

dev.off()