require(Hmisc) # for summarize()
require(RColorBrewer)

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

############################################
## Summarize community change by grid cell
############################################

load('data/biomassavemap_testK6noSeas.RData') # loads biomassavemap data.frame
	dim(biomassavemap)

# summarize across seasons within each region
biomassavemap.maxseas <- aggregate(list(wtcpue.proj = biomassavemap$wtcpue.proj), by=list(sppocean=biomassavemap$sppocean, region=biomassavemap$region, period=biomassavemap$period, lat=biomassavemap$lat, lon=biomassavemap$lon), FUN=max) # max across any season. SLOW (a few minutes)
	dim(biomassavemap.maxseas)

biomassavemap.meanseas <- aggregate(list(wtcpue.proj = biomassavemap$wtcpue.proj), by=list(sppocean=biomassavemap$sppocean, region=biomassavemap$region, period=biomassavemap$period, lat=biomassavemap$lat, lon=biomassavemap$lon), FUN=mean) # mean across all seasons. SLOW (a few minutes)
	dim(biomassavemap.meanseas)

# find abundance threshold to count as present for each taxon
# use cumulative 5% of wtcpue from earliest timeperiod, by region
i <- !duplicated(biomassavemap[,c('region', 'sppocean')])
taxthresh <- biomassavemap[i,c('region', 'sppocean')]
	taxthresh$thresh<-NA
for(i in 1:nrow(taxthresh)){
	print(paste(i, 'of', nrow(taxthresh)))
	wts <- sort(biomassavemap.meanseas$wtcpue.proj[biomassavemap.meanseas$sppocean == taxthresh$sppocean[i] & biomassavemap.meanseas$region == taxthresh$region[i] & biomassavemap.meanseas$period == '2006-2020'])
	ind<-max(which(cumsum(wts)/sum(wts)<0.05))
	taxthresh$thresh[i] <- wts[ind]
}

# apply threshold to projections to determine presence/absence
presmap <- merge(biomassavemap.meanseas, taxthresh) # surprisingly fast
presmap$pres <- presmap$wtcpue.proj > presmap$thresh
	
	# examine the results
	table(presmap$sppocean, presmap$pres) # how many marked present vs. absent

	i = presmap$sppocean == 'gadus morhua_Atl' & presmap$period == '2006-2020'
	i = presmap$sppocean == 'gadus morhua_Atl' & presmap$period == '2081-2100'
	i = presmap$sppocean == 'theragra chalcogramma_Pac' & presmap$period == '2006-2020'
	i = presmap$sppocean == 'theragra chalcogramma_Pac' & presmap$period == '2081-2100'
	plot(presmap$lon[i], presmap$lat[i], col=c('black', 'red')[presmap$pres[i]+1])
	
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
	
	
################
## Plots
################
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