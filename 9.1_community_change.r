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

####################
## helper functions
####################
require(Hmisc) # for summarize()
require(RColorBrewer)

lu <- function(x) return(length(unique(x)))

# expects x to have columns 'period' and 'rich' (in that order)
# each as used below (period has a dash, and we use the midpoint here)
calcrichtrend <- function(x){
	pds <- as.numeric(unlist(strsplit(as.character(x[,1]), split='-')))
	dim(pds) <- c(2,nrow(x))
	mids <- colMeans(pds)
	mod <- lm(x[,2] ~ mids)
	return(mod$coefficients[2])
}

# normalize to 0-1
norm01 <- function(x){
	mn <- min(x)
	mx <- max(x)
	return((x-mn)/(mx-mn))
}

# turnover metrics
# for use with Hmisc::summarize() but doesn't work for some reason
# expects region in column 1, lat (2), lon (3), time period in column 4, spp in col 5, pres in col 6
# only present species included
turnover <- function(x){
	pers <- sort(unique(x[,4]))
	if(length(pers) != 2) stop(paste('expecting 2 time periods:', unique(x[,1]), unique(x[,2]), unique(x[,3]), paste(unique(x[,4], collapse=','))))
	initcomm <- x[x[,4] == pers[1] & x[,6]==TRUE,5]
	finalcomm <- x[x[,4] == pers[2] & x[,6]==TRUE,5]
	
	a <- length(intersect(initcomm, finalcomm)) # number of shared species
	lost <- length(setdiff(initcomm, finalcomm)) # number of species lost
	gained <- length(setdiff(finalcomm, initcomm))
	return(c(nstart=length(initcomm), nend=length(finalcomm), nlost=lost, ngained=gained, flost=lost/length(initcomm), fgained=gained/length(initcomm), beta_sor=2*a/(2*a+gained+lost)))
}

############################################
## Summarize community change by grid cell
############################################

load('data/biomassavemap_testK6noSeas.RData') # loads biomassavemap data.frame
	dim(biomassavemap)

# summarize across seasons within each region
#biomassavemap.maxseas <- aggregate(list(wtcpue.proj = biomassavemap$wtcpue.proj), by=list(sppocean=biomassavemap$sppocean, region=biomassavemap$region, period=biomassavemap$period, lat=biomassavemap$lat, lon=biomassavemap$lon), FUN=max) # max across any season. SLOW (a few minutes)
#	dim(biomassavemap.maxseas)

biomassavemap.meanseas <- aggregate(list(wtcpue.proj = biomassavemap$wtcpue.proj), by=list(sppocean=biomassavemap$sppocean, region=biomassavemap$region, period=biomassavemap$period, lat=biomassavemap$lat, lon=biomassavemap$lon), FUN=mean) # mean across all seasons. SLOW (a few minutes)
	dim(biomassavemap.meanseas)
	summary(biomassavemap.meanseas) # NAs, presumably where projections only available in one season
	
	# trim out NAs
	biomassavemap.meanseas <- biomassavemap.meanseas[!is.na(biomassavemap.meanseas$wtcpue.proj),]

# find abundance threshold to count as present for each taxon
# use cumulative 5% of wtcpue from earliest timeperiod, by region
# Quite slow: 10 min or so
i <- !duplicated(biomassavemap[,c('region', 'sppocean')])
taxthresh <- biomassavemap[i,c('region', 'sppocean')]
	taxthresh$thresh<-NA
for(i in 1:nrow(taxthresh)){
	print(paste(i, 'of', nrow(taxthresh)))
	wts <- sort(biomassavemap.meanseas$wtcpue.proj[biomassavemap.meanseas$sppocean == taxthresh$sppocean[i] & biomassavemap.meanseas$region == taxthresh$region[i] & biomassavemap.meanseas$period == '2006-2020'])
	ind<-max(which(cumsum(wts)/sum(wts)<0.05))
	taxthresh$thresh[i] <- wts[ind]
}
	
	save(taxthresh, file='data/taxthresh.RData')

# apply threshold to projections to determine presence/absence
presmap <- merge(biomassavemap.meanseas, taxthresh) # surprisingly fast (a minute or so)
presmap$pres <- presmap$wtcpue.proj > presmap$thresh
	
	# examine the results
	table(presmap$sppocean, presmap$pres) # how many marked present vs. absent

	i = presmap$sppocean == 'gadus morhua_Atl' & presmap$period == '2006-2020'
	i = presmap$sppocean == 'gadus morhua_Atl' & presmap$period == '2081-2100'
	i = presmap$sppocean == 'theragra chalcogramma_Pac' & presmap$period == '2006-2020'
	i = presmap$sppocean == 'theragra chalcogramma_Pac' & presmap$period == '2081-2100'
	plot(presmap$lon[i], presmap$lat[i], col=c('black', 'red')[presmap$pres[i]+1])

	# write out
	save(presmap, file='data/presmap.RData')
	
# richness by grid cell by time period (across all seasons)
i <- presmap$pres # only count where a spp is present
rich <- aggregate(list(rich = presmap$sppocean[i]), by=list(region=presmap$region[i], period=presmap$period[i], lat=presmap$lat[i], lon=presmap$lon[i]), FUN=lu) # calculate richness (# taxa) by grid cell in each time period

	# examine
	nrich <- aggregate(list(nrich = rich$rich), by=list(region=rich$region, lat=rich$lat, lon=rich$lon), FUN=lu) # only ever has 1 value of richness per grid cell. that's a problem
		summary(nrich) # up to five values: good
	
	i<- rich$lat == 52.875 & rich$lon == 170.625
	plot(rich$period[i], rich$rich[i])
	
	regs <- unique(rich$region)
	allgrids <- paste(rich$lat, rich$lon)
	col.ln <- rgb(0.1, 0.1, 0.1, 0.5)
	quartz(width=5,height=5)
	# pdf(width=5, height=5, file='figures/richness_proj_by_grid.pdf')
	par(mfrow=c(3,4), mai=c(0.4, 0.4, 0.2, 0.1), mgp=c(2, 0.6, 0), las=1)
	for(i in 1:length(regs)){
		print(i)
		inds1 <- rich$region == regs[i]
		mxr <- max(rich$rich[inds1])
		thesegrids <- unique(allgrids[inds1])

		inds2 <- inds1 & allgrids == thesegrids[1]
		plot(as.numeric(rich$period)[inds2], rich$rich[inds2], type='l', xlab='Period', ylab='Richness', ylim=c(0,mxr), main=regs[i], cex.main=0.8, cex.axis=0.8, col=col.ln, lwd=0.1)

		for(j in 2:length(thesegrids)){
			inds2 <- inds1 & allgrids == thesegrids[j]
			lines(as.numeric(rich$period)[inds2], rich$rich[inds2], col=col.ln, lwd=0.1)
		}
	}
	
	dev.off()
	
	# write out
	save(rich, file='data/rich.RData')
	
# load in richness calcs
load('data/rich.RData') # loads rich data.frame

	
# calculate trend in richness by grid cell
richtrend <- Hmisc::summarize(X=rich[,c('period', 'rich')], by=list(region=rich$region, lat=rich$lat, lon=rich$lon, dummy=rich$lon), FUN=calcrichtrend, stat.name='trend') # calculate richness (# taxa) by grid cell in each time period. for an unknown resion, the last column in the by list gets NAs, so padded with a dummy column here
	richtrend <- richtrend[,-grep('dummy', names(richtrend))]
	
	# examine
	hist(richtrend$trend) # nicely centered around 0. perhaps a slightly longer tail to the left (negative)
	
	# write out
	save(richtrend, file='data/richtrend.RData')



# calculate turnover metrics for each grid cell from start to end period
turn <- rich[rich$period == '2006-2020',c('region', 'lat', 'lon', 'rich')] # initial richness
	names(turn)[names(turn)=='rich'] <- 'rich_init'
turn <- merge(turn, rich[rich$period == '2081-2100',c('region', 'lat', 'lon', 'rich')])
	names(turn)[names(turn)=='rich'] <- 'rich_final'

	# efficient, but not working
#presmap <- presmap[order(presmap$region, presmap$period, presmap$lat, presmap$lon, presmap$sppocean),]
#inds <- presmap$period %in% c('2006-2020', '2081-2100')
#temp <- Hmisc::summarize(X=presmap[inds,c('region', 'lat', 'lon', 'period', 'sppocean', 'pres')], by=list(region=presmap$region[inds], lat=presmap$lat[inds], lon=presmap$lon[inds]), FUN=turnover, stat.name=NULL)

	# slow: takes an hour
turn$nstart <- turn$nend <- turn$nlost <- turn$ngained <- turn$flost <- turn$fgained <- turn$beta_sor <- NA
for(i in 1:nrow(turn)){
	if(i %% 100 == 0) print(paste(i, 'of', nrow(turn), Sys.time()))
	x<-presmap[presmap$region==turn$region[i] & presmap$period %in% c('2006-2020', '2081-2100') & presmap$lat==turn$lat[i] & presmap$lon==turn$lon[i],c('region', 'lat', 'lon', 'period', 'sppocean', 'pres')]
	turn[i,c('nstart', 'nend', 'nlost', 'ngained', 'flost', 'fgained', 'beta_sor')] <- turnover(x)
}
	
	# write out
	save(turn, file='data/turn.RData')
	
	
################
## Plots
################
load('data/rich.RData')
load('data/turn.RData')

# plot maps of richness (all of North America)
	colfun <- colorRamp(rev(brewer.pal(11, 'Spectral')))
	periods <- sort(unique(rich$period))
	quartz(width=7, height=5)
	# pdf(width=7, height=5, file='figures/richness_proj_map.pdf')
	for(i in 1:length(periods)){
		par(mai=c(0.7,0.7,0.3, 0.1), las=1, mgp=c(2,1,0))
		j = rich$period == periods[i]
		plot(rich$lon[j], rich$lat[j], col=rgb(colfun(norm01(rich$rich[j])), maxColorValue=255), pch=16, cex=0.3, xlab='Longitude', ylab='Latitude', main=paste('Projected species richness', periods[i]))
		legend('bottomleft', legend=round(seq(min(rich$rich[j]), max(rich$rich[j]), length.out=10),2), col=rgb(colfun(norm01(seq(min(rich$rich[j]), max(rich$rich[j]), length.out=10))), maxColorValue=255), pch=16, cex=0.8, title='Species', bty='n')
	}
	
	dev.off()

# plot maps of richness (by region)
	colfun <- colorRamp(rev(brewer.pal(11, 'Spectral')))
	regs <- sort(unique(rich$region))
	periods <- sort(unique(rich$period))
	cexs = c(0.7, 1.2, 0.5, 0.8, 0.5, 1.1, 1, 1.5, 0.5, 0.8, 1) # to adjust for each region
	quartz(width=14, height=3)
	# pdf(width=14, height=3, file='figures/richness_proj_map_byregion.pdf')
	for(k in 1:length(regs)){
		par(mfrow=c(1,length(periods)), mai=c(0.7,0.7,0.3, 0.1), las=1, mgp=c(2,1,0))
		for(i in 1:length(periods)){
			j = rich$period == periods[i] & rich$region == regs[k]
			plot(rich$lon[j], rich$lat[j], col=rgb(colfun(norm01(rich$rich[j])), maxColorValue=255), pch=16, cex=cexs[k], xlab='Longitude', ylab='Latitude', main=paste(regs[k], periods[i]))
			legend('bottomleft', legend=round(seq(min(rich$rich[j]), max(rich$rich[j]), length.out=10),2), col=rgb(colfun(norm01(seq(min(rich$rich[j]), max(rich$rich[j]), length.out=10))), maxColorValue=255), pch=16, cex=0.6, title='Species', bty='n')
		}
	}	
	dev.off()

	
# plot map of richness trend
	colfun <- colorRamp(brewer.pal(11, 'Spectral'))
	quartz(width=7, height=5)
	# pdf(width=7, height=5, file='figures/richness_proj_trend_map.pdf')
	par(mai=c(0.7,0.7,0.3, 0.1), las=1, mgp=c(2,1,0))
	plot(richtrend$lon, richtrend$lat, col=rgb(colfun(norm01(richtrend$trend)), maxColorValue=255), pch=16, cex=0.5, xlab='Longitude', ylab='Latitude', main='Change in species richness')
	legend('bottomleft', legend=round(seq(min(richtrend$trend), max(richtrend$trend), length.out=10),2), col=rgb(colfun(norm01(seq(min(richtrend$trend), max(richtrend$trend), length.out=10))), maxColorValue=255), pch=16, cex=0.8, title='Species/year', bty='n')

	dev.off()
	

# plot maps of turnover (by region)
	colfun <- colorRamp(rev(brewer.pal(11, 'Spectral')))
	regs <- sort(unique(turn$region))
	cexs = c(0.7, 1.2, 0.5, 0.8, 0.5, 1.1, 1, 1.5, 0.5, 0.8, 1) # to adjust for each region


	#quartz(width=4, height=3)
	pdf(width=4, height=3, file='figures/turnover_proj_map_byregion_beta_sor.pdf')
	for(k in 1:length(regs)){
		par(mai=c(0.7,0.7,0.3, 0.1), las=1, mgp=c(2,1,0))
		j = turn$region == regs[k]
		plot(turn$lon[j], turn$lat[j], col=rgb(colfun(turn$beta_sor[j]), maxColorValue=255), pch=16, cex=cexs[k], xlab='Longitude', ylab='Latitude', main=regs[k])
		legend('bottomleft', legend=seq(0,1, length.out=11), col=rgb(colfun(seq(0, 1, length.out=11)), maxColorValue=255), pch=16, cex=0.6, title='Sorenson', bty='n')
	}	
	dev.off()

	#quartz(width=4, height=3)
	pdf(width=4, height=3, file='figures/turnover_proj_map_byregion_flost.pdf')
	for(k in 1:length(regs)){
		par(mai=c(0.7,0.7,0.3, 0.1), las=1, mgp=c(2,1,0))
		j = turn$region == regs[k]
		plot(turn$lon[j], turn$lat[j], col=rgb(colfun(turn$flost[j]), maxColorValue=255), pch=16, cex=cexs[k], xlab='Longitude', ylab='Latitude', main=regs[k])
		legend('bottomleft', legend=seq(0,1, length.out=11), col=rgb(colfun(seq(0, 1, length.out=11)), maxColorValue=255), pch=16, cex=0.6, title='Fraction lost', bty='n')
	}	
	dev.off()


	#quartz(width=4, height=3)
	pdf(width=4, height=3, file='figures/turnover_proj_map_byregion_fgained.pdf')
	y<-turn$fgained
	y[y>1] <- 1
	for(k in 1:length(regs)){
		par(mai=c(0.7,0.7,0.3, 0.1), las=1, mgp=c(2,1,0))
		j = turn$region == regs[k]
		plot(turn$lon[j], turn$lat[j], col=rgb(colfun(y[j]), maxColorValue=255), pch=16, cex=cexs[k], xlab='Longitude', ylab='Latitude', main=regs[k])
		legend('bottomleft', legend=c(seq(0,0.9, length.out=10),'>=1'), col=rgb(colfun(seq(0, 1, length.out=11)), maxColorValue=255), pch=16, cex=0.6, title='Fraction gained', bty='n')
	}	
	dev.off()
	
