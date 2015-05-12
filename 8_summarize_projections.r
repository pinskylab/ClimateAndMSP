## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	projfolder = '../CEmodels_proj' # holds model projections
	modfolder = '../CEModels' # holds the models
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	projfolder = 'CEmodels_proj'
	modfolder = 'CEmodels'
	}
# could add code for Lauren's working directory here



##################################################
## Calculate mean lat/depth and sum of biomass  ##
##################################################
require(mgcv)
require(Hmisc)

## choose which run to use
runtype <- 'test'

# load all projections from this run
files <- list.files(path = projfolder, pattern=paste('summproj_', runtype, '_', sep=''))

# weighted mean function to use with summarize()
wmean <- function(x) return(weighted.mean(x[,1], x[,2])) # values in col 1, weights in col 2

# save observed and predicted positions for lat, lon, and depth
meanpos <- list(0)


# set up dataframes
# Don't know how many regions for each taxon, so can't pre-allocate the rows (but could guess high... might speed this up)
# could also calc mean depth, but would need to pull from rugosity file
biomassave = data.frame(sppocean = character(0), region = character(0), year = numeric(0)) # sum of average wtcpue across the region (2020-2099) for each survey in each model (columns)
	for(i in 1:13) biomassave[[paste('summwtcpue', i, sep='_')]] <- numeric(0)
meanlat = data.frame(sppocean = character(0), region = character(0), year = numeric(0)) # biomass-weighted mean latitude across the region (2020-2099) for each survey in each model (columns)
	for(i in 1:13) meanlat[[paste('lat', i, sep='_')]] <- numeric(0)
meanlon = data.frame(sppocean = character(0), region = character(0), year = numeric(0))
	for(i in 1:13) meanlon[[paste('lon', i, sep='_')]] <- numeric(0)


options(warn=1) # print warnings as they occur

# loop through all files
for(i in 1:length(files)){ # takes a while (a couple hours ?)
	# load data for this species
	print(paste(i, 'of', length(files), Sys.time()))
	load(paste(projfolder, '/', files[i], sep='')) # load summproj for this taxon
	myregions <- sort(unique(summproj$region))
	mysppocean <- gsub('.Rdata', '', gsub(paste('summproj_', runtype, '_', sep=''), '', files[i]))

	# set up dataframes for this taxon
	mybiomassave = data.frame(sppocean = mysppocean, region = rep(myregions, rep(80,length(myregions))), year = rep(2020:2099, length(myregions)))
	mymeanlat = data.frame(sppocean = mysppocean, region = rep(myregions, rep(80,length(myregions))), year = rep(2020:2099, length(myregions)))
	mymeanlon = data.frame(sppocean = mysppocean, region = rep(myregions, rep(80,length(myregions))), year = rep(2020:2099, length(myregions)))

	# Summarize projections
	for(j in 1:13){
		snm <- paste('wtcpue.proj', j, sep='_')
		temp <- aggregate(summproj[,snm], by=list(year = summproj$year, region = summproj$region), FUN=sum)
			names(temp)[3] <- paste('summwtcpue', j, sep='_')
		mybiomassave <- merge(mybiomassave, temp)

		temp <- summarize(summproj[,c('lat', snm)], by=list(region = summproj$region, year = summproj$year), FUN=wmean)
			names(temp)[3] <- paste('lat', j, sep='_')
		mymeanlat <- merge(mymeanlat, temp)

		temp <- summarize(summproj[,c('lon', snm)], by=list(region = summproj$region, year = summproj$year), FUN=wmean)
			names(temp)[3] <- paste('lon', j, sep='_')
		mymeanlon <- merge(mymeanlon, temp)

	}

	biomassave <- rbind(biomassave, mybiomassave) # an inefficient way to do this: better to pre-allocate
	meanlat <- rbind(meanlat, mymeanlat)
	meanlon <- rbind(meanlon, mymeanlon)

}


### Save meanpos and biomassave
save(biomassave, meanlat, meanlon, file = paste('data/meanlat,lon,biomass_', runtype, '.RData', sep=''))




####################################################
## Plot summary data about spp biomass and ranges ##
####################################################
require(Hmisc)

wmean <- function(x){
	i <- !is.na(x[,1]) & !is.na(x[,2])
	if(sum(i) == 0) return(NA)
	else return(weighted.mean(x[i,1], x[i,2])) # weighted mean function for summarize()
}

## load the data
runtype <- 'test'
load(paste('data/meanlat,lon,biomass_', runtype, '.RData', sep='')) # biomassave, meanlat, meanlon: projected biomass by year for each taxon in each region

## plots of change in biomass
	sppregions <- sort(unique(paste(biomassave$region, biomassave$sppocean)))

	# quartz(width=10, height=8)
	pdf(file=paste('figures/biomassave_proj_', runtype, Sys.Date(), '.pdf', sep=''), width=10, height=8)
	par(mfrow = c(6,6), mai=c(0.3, 0.3, 0.2, 0.05), cex.main=0.7, cex.axis=0.8, omi=c(0,0.2,0.1,0), mgp=c(2.8, 0.7, 0), font.main=3)

	rc = 1 # row counter
	cc = 0 # column counter

	options(warn=1) # print warnings as they occur
	for(i in 1:length(sppregions)){
		inds <- paste(biomassave$region, biomassave$sppocean) == sppregions[i]
		thisreg <- unique(biomassave$region[inds])
		thisspp <- unique(biomassave$sppocean[inds])

		# increment row and column counters as needed
		cc = cc+1
		if(cc == 7){ cc = 1; rc = rc + 1}
		if(rc == 7){ cc = 1; rc = 1}
	
		if(i>1){ if(thisreg != oldreg){  # switch to a new page when I get to a new region
				par(mfrow = c(6,6), mai=c(0.3, 0.3, 0.2, 0.05), cex.main=0.7, cex.axis=0.8, omi=c(0,0.2,0.1,0), mgp=c(2.8, 0.7, 0), font.main=3)
				rc = 1; cc = 1
		}}
		
		# plot data for each GCM
		plot(biomassave$year[inds], biomassave$summwtcpue_1[inds], col='grey', las=1, type='l')
		for(j in 2:13) lines(biomassave$year[inds], biomassave[[paste('summwtcpue', j, sep='_')]][inds], col='grey')
		
		# plot ensemble mean
		agg <- cbind(biomassave$year, rowMeans(biomassave[inds, grep('sumwtcpue', names(biomassave))])
		lines(agg$year, agg[,2], col='black')
	
		if(cc==1) mtext(text='Biomass', side=2, line=2.3, cex=0.6) # add y label on left of each row
		if(rc==1) mtext(text=thisreg, side=3, line=1.3, cex=0.6) # add region header on top of page

		oldreg = thisreg # save the previous region to see if we need a new page on the next round
	}

	dev.off()


## plots of change in latitude
	#load('Output/meanpos_allyrsonlytemp_2013-11-10.RData')

	# quartz(width=10, height=8)
	pdf(file=paste('Figures/avelat_proj_', Sys.Date(), '.pdf', sep=''), width=10, height=8)
	par(mfrow = c(6,6), mai=c(0.3, 0.3, 0.2, 0.05), cex.main=0.7, cex.axis=0.8, omi=c(0,0.2,0.1,0), mgp=c(2.8, 0.7, 0), font.main=3)

	rc = 1 # row counter
	cc = 0 # column counter
	for(s in 1:length(files)){
		print(paste(s, 'of', length(files)))

		thisspp = spp[sapply(spp, FUN=grepl, x=files[s], fixed=TRUE)] # extract spp name from file name
		if(length(thisspp)>1){ # if there were multiple matches, take the one that was longest
			i = which.max(nchar(as.character(thisspp)))
			print(paste('picked', thisspp[i], 'from', paste(thisspp, collapse=', '), 'for', files[s]))
			thisspp = thisspp[i]
		}
		thisreg = regs[sapply(regs, FUN=grepl, x=files[s])]

		#j = which(grepl(thisspp, names(meanpos), fixed=TRUE) & grepl(thisreg, names(meanpos))) # calculate the correct index into meanpos for this spp and region
		thisproj = read.csv(paste(folder, '/', files[s], sep=''), row.names=1)
		meanlat = summarize(thisproj[,c('lat', 'wtcpue.hist')], by=thisproj$year, FUN=wmean)
		for(i in 1:13)	meanlat = cbind(meanlat, summarize(thisproj[,c('lat', paste('wtcpue.proj_', i, sep=''))], by=thisproj$year, FUN=wmean)[,2])
		names(meanlat) = c('year', 'hist', paste('m', 1:13, sep=''))

		i = data$spp == as.character(thisspp) & data$region == as.character(thisreg) & complete.cases(data[,c('surftemp', 'bottemp')])
		if(thisreg=='DFO_Newfoundland_Fall') i = data$spp == as.character(thisspp) & data$region == as.character(thisreg) & complete.cases(data[,c('bottemp')])
		obs = summarize(data[i,c('lat', 'wtcpue')], by=data$year[i], FUN=wmean)

		# increment row and column counters as needed
		cc = cc+1
		if(cc == 7){ cc = 1; rc = rc + 1}
		if(rc == 7){ cc = 1; rc = 1}
	
		if(s>1){ if(thisreg != oldreg){  # switch to a new page when I get to a new region
				par(mfrow = c(6,6), mai=c(0.3, 0.3, 0.2, 0.05), cex.main=0.7, cex.axis=0.8, omi=c(0,0.2,0.1,0), mgp=c(2.8, 0.7, 0), font.main=3)
				rc = 1; cc = 1
		}}

		#boxplot(as.data.frame(meanpos[[j]]$meanlat), names=c('Obs', 'Hist', '2060', '2100'), las=1, main=thisspp, col='grey')
		ylims = range(c(obs[,2], meanlat[,2:15]), na.rm=TRUE)
		matplot(meanlat$year, meanlat[,2:15], type='l', col='grey', lty=1, ylim=ylims, las=1, main=thisspp)
		lines(obs[,1], obs[,2], col='black')		
		lines(meanlat$year, rowMeans(meanlat[,3:15], na.rm=TRUE), col='blue')
		if(cc==1) mtext(text='Latitude (°N)', side=2, line=2.3, cex=0.6) # add y label on left of each row
		if(rc==1) mtext(text=thisreg, side=3, line=1.3, cex=0.6) # add region header on top of page

		oldreg = thisreg # save the previous region to see if we need a new page on the next round
	}

	dev.off()



	



##################################################################################
## Summarize position shifts within each assemblage for 2020-2060 and 2060-2100 ##
##################################################################################
setwd('/Users/mpinsky/Documents/Rutgers/Range projections/')
load('Output/meanpos_allyrsonlytemp_2013-11-10_2.RData')
folder = 'Output/proj_annual_2013-11-08/'
dataclim = read.csv('Output/dataclim_2013-10-23.csv', row.names=1) # just to get spp names
	spp = sort(unique(dataclim$taxon))
	regs = sort(unique(dataclim$region))
	
files = list.files(folder, pattern='biomassave_')
	
# find median position projections for each taxon
df = data.frame(region = character(0), spp = character(0), clim=numeric(0), hist=numeric(0), X2020=numeric(0), X2060=numeric(0)) # 2020-2060 and 2060-2100
medianpos = list(medianlat = df, medianlon = df, mediandepth = df, medianbio = df)

for(s in 1:length(meanpos)){
	thisspp = spp[sapply(spp, FUN=grepl, x=names(meanpos)[s], fixed=TRUE)] # extract spp name from meanpos name
	if(length(thisspp)>1){ # if there were multiple matches, take the one that was longest
		i = which.max(nchar(as.character(thisspp)))
		print(paste('picked', thisspp[i], 'from', paste(thisspp, collapse=', '), 'for', names(meanpos)[s]))
		thisspp = thisspp[i]
	}
	thisreg = regs[sapply(regs, FUN=grepl, x=names(meanpos)[s])]

	newlat = data.frame(t(apply(meanpos[[s]]$meanlat, MARGIN=2, FUN=median, na.rm=TRUE)))
		newlat$spp = thisspp
		newlat$region = thisreg
	newlon = data.frame(t(apply(meanpos[[s]]$meanlon, MARGIN=2, FUN=median, na.rm=TRUE)))
		newlon$spp = thisspp
		newlon$region = thisreg
	newdepth = data.frame(t(apply(meanpos[[s]]$meandepth, MARGIN=2, FUN=median, na.rm=TRUE)))
		newdepth$spp = thisspp
		newdepth$region = thisreg
	
	medianpos$medianlat = rbind(medianpos$medianlat, newlat)
	medianpos$medianlon = rbind(medianpos$medianlon, newlon)
	medianpos$mediandepth = rbind(medianpos$mediandepth, newdepth)
	
	# biomass
	biomassave = read.csv(paste(folder, files[grepl(paste(thisreg, '_', thisspp, sep=''), files, fixed=TRUE)], sep=''), row.names=1)
	newbio = data.frame(t(apply(biomassave, MARGIN=2, FUN=median, na.rm=TRUE)))
		names(newbio)[1:2] = c('obs', 'clim')
		newbio$spp = thisspp
		newbio$region = thisreg
	medianpos$medianbio = rbind(medianpos$medianbio, newbio)
}


# summarize into position shifts from historical to 2020-2060 or to 2060-2100
df = data.frame(region = medianpos$medianlat$region, spp = medianpos$medianlat$spp) # shifts from historical projection to future time periods
posshift = list(latshift = df, lonshift = df, depshift = df)
posshift$latshift$X2020 = medianpos$medianlat$X2020 - medianpos$medianlat$clim
posshift$latshift$X2060 = medianpos$medianlat$X2060 - medianpos$medianlat$clim
posshift$lonshift$X2020 = medianpos$medianlon$X2020 - medianpos$medianlon$clim
posshift$lonshift$X2060 = medianpos$medianlon$X2060 - medianpos$medianlon$clim
posshift$depshift$X2020 = medianpos$mediandepth$X2020 - medianpos$mediandepth$clim
posshift$depshift$X2060 = medianpos$mediandepth$X2060 - medianpos$mediandepth$clim

# t-tests of mean shift over all species in the assemblage
regs = sort(unique(posshift$latshift$region))
ttestlat2020 = matrix(NA, nrow=6, ncol=9)
	rownames(ttestlat2020) = c('mean', 'u95', 'l95', 't', 'df', 'p'); colnames(ttestlat2020) = regs
for(i in 1:9){
	mod = t.test(posshift$latshift$X2020[posshift$latshift$region==regs[i]])
	ttestlat2020[,i] = c(mod$estimate, mod$conf.int, mod$statistic, mod$parameter, mod$p.value)		
} 

ttestlon2020 = matrix(NA, nrow=6, ncol=9)
	rownames(ttestlon2020) = c('mean', 'u95', 'l95', 't', 'df', 'p'); colnames(ttestlon2020) = regs
for(i in 1:9){
	mod = t.test(posshift$lonshift$X2020[posshift$lonshift$region==regs[i]])
	ttestlon2020[,i] = c(mod$estimate, mod$conf.int, mod$statistic, mod$parameter, mod$p.value)		
} 

ttestdepth2020 = matrix(NA, nrow=6, ncol=9)
	rownames(ttestdepth2020) = c('mean', 'u95', 'l95', 't', 'df', 'p'); colnames(ttestdepth2020) = regs
for(i in 1:9){
	mod = t.test(posshift$depshift$X2020[posshift$depshift$region==regs[i]])
	ttestdepth2020[,i] = c(mod$estimate, mod$conf.int, mod$statistic, mod$parameter, mod$p.value)		
} 

ttestlat2060 = matrix(NA, nrow=6, ncol=9)
	rownames(ttestlat2060) = c('mean', 'u95', 'l95', 't', 'df', 'p'); colnames(ttestlat2060) = regs
for(i in 1:9){
	mod = t.test(posshift$latshift$X2060[posshift$latshift$region==regs[i]])
	ttestlat2060[,i] = c(mod$estimate, mod$conf.int, mod$statistic, mod$parameter, mod$p.value)		
} 

ttestlon2060 = matrix(NA, nrow=6, ncol=9)
	rownames(ttestlon2060) = c('mean', 'u95', 'l95', 't', 'df', 'p'); colnames(ttestlon2060) = regs
for(i in 1:9){
	mod = t.test(posshift$lonshift$X2060[posshift$lonshift$region==regs[i]])
	ttestlon2060[,i] = c(mod$estimate, mod$conf.int, mod$statistic, mod$parameter, mod$p.value)		
} 

ttestdepth2060 = matrix(NA, nrow=6, ncol=9)
	rownames(ttestdepth2060) = c('mean', 'u95', 'l95', 't', 'df', 'p'); colnames(ttestdepth2060) = regs
for(i in 1:9){
	mod = t.test(posshift$depshift$X2060[posshift$depshift$region==regs[i]])
	ttestdepth2060[,i] = c(mod$estimate, mod$conf.int, mod$statistic, mod$parameter, mod$p.value)		
} 


# plot position shifts (lat/lon/depth)
plotstar = function(ttest){
	x = (1:9)[ttest[6,]<0.05]
	lims = par('usr')
	y = rep(lims[4] - 0.1*(lims[4]-lims[3]), length(x))
	pchs = (c('-', '+')[(ttest[1,]>0)+1])[ttest[6,]<0.05]
	points(x, y, pch=pchs, col='red', cex=2)
}

quartz(width=9, height=6)
# pdf(width=9, height=6, file=paste('Figures/assemblage_shift_boxplots_', Sys.Date(), '.pdf', sep=''))
par(mfrow=c(2,3), las=2, mgp = c(2,0.8,0), mai=c(0.85,0.4, 0.3, 0.1))
regnms = c('Aleutians', 'E. Bering', 'Gulf AK', 'Newfoundland', 'Scotian', 'Gulf St. Law.', 'Northeast', 'Gulf Mex.', 'West Coast')
boxplot(X2020 ~ region, data = posshift$latshift, names=regnms, ylab='Change in lat (°)', main='Latitude to 2020-2060', range=0)
	abline(h=0, col='grey'); plotstar(ttestlat2020)
boxplot(X2020 ~ region, data = posshift$lonshift, names=regnms, ylab='Change in lon (°)', main='Longitude to 2020-2060', range=0)
	abline(h=0, col='grey'); plotstar(ttestlon2020)
boxplot(X2020 ~ region, data = posshift$depshift, names=regnms, ylab='Change in depth (m)', main='Depth to 2020-2060', range=0)
	abline(h=0, col='grey'); plotstar(ttestdepth2020)
boxplot(X2060 ~ region, data = posshift$latshift, names=regnms, ylab='Change in lat (°)', main='Latitude to 2060-2100', range=0)
	abline(h=0, col='grey'); plotstar(ttestlat2060)
boxplot(X2060 ~ region, data = posshift$lonshift, names=regnms, ylab='Change in lon (°)', main='Longitude to 2060-2100', range=0)
	abline(h=0, col='grey'); plotstar(ttestlon2060)
boxplot(X2060 ~ region, data = posshift$depshift, names=regnms, ylab='Change in depth (m)', main='Depth to 2060-2100', range=0)
	abline(h=0, col='grey'); plotstar(ttestdepth2060)
	
dev.off()




##########################################################
## Plot maps of spp projections from each climate model 
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

