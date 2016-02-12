#######################################################
# Read in results and analyze against ensemble mean
#######################################################
consplan1 <- read.csv(paste(marxfolder, 'output', runname, '/', runname1, '_best.csv', sep=''))
consplan2 <- read.csv(paste(marxfolder, 'output', runname2, '/', runname2, '_best.csv', sep=''))
load(paste(inputfolder, '/spps.Rdata', sep=''))

load(paste('data/rich_', runtype, projtype, '_rcp', rcp, '.RData', sep='')) # loads rich data.frame with presence/absence information
load(paste('data/presmap_', runtype, projtype, '_rcp', rcp, '.RData', sep='')) # loads presmap data.frame with presence/absence information

# add zone to pus
pusplan <- merge(pus, consplan, by.x='id', by.y='planning_unit')
	dim(pus)
	dim(pusplan)

pusplan2 <- merge(pus2, consplan2, by.x='id', by.y='planning_unit')
	dim(pus2)
	dim(pusplan2)

# plot map of selected grids, on top of richness (#1)
	colfun <- colorRamp(rev(brewer.pal(11, 'Spectral')))
	cexs = 0.5 # to adjust
	periods <- sort(unique(rich$period))
	# quartz(width=10, height=3)
	pdf(width=10, height=3, file=paste('figures/MarZone_NEUSSpring_on_richness_', runname, '.pdf', sep=''))
	par(mfrow=c(1,length(periods)), mai=c(0.5,0.5,0.3, 0.1), las=1, mgp=c(2,1,0))
	j <- rich$region == 'NEFSC_NEUSSpring'
	for(i in 1:length(periods)){
		j2 <- rich$period == periods[i] & j
		plot(rich$lon[j2], rich$lat[j2], col=rgb(colfun(pnorm01(rich$rich[j2], rich$rich[j])), maxColorValue=255), pch=16, cex=cexs, xlab='Longitude', ylab='Latitude', main=paste('NEUSSpring', periods[i]), cex.main=0.9)

		i <- pusplan$zone == 2
		points(pusplan$lon[i], pusplan$lat[i], pch=1, col='black', cex=cexs)
	}
	legend('bottomright', legend=round(seq(min(rich$rich[j]), max(rich$rich[j]), length.out=10),2), col=rgb(colfun(norm01(seq(min(rich$rich[j]), max(rich$rich[j]), length.out=10))), maxColorValue=255), pch=16, cex=0.7, title='Species', bty='n')


	dev.off()

# plot map of selected grids, on top of richness (#1 and #2)
	colfun <- colorRamp(rev(brewer.pal(11, 'Spectral')))
	cexs = 0.5 # to adjust
	periods <- sort(unique(rich$period))
	# quartz(width=10, height=3)
	pdf(width=10, height=3, file=paste('figures/MarZone_NEUSSpring_on_richness_', runname, '&', runname2, '.pdf', sep=''))
	par(mfrow=c(1,length(periods)), mai=c(0.5,0.5,0.3, 0.1), las=1, mgp=c(2,1,0))
	j <- rich$region == 'NEFSC_NEUSSpring'
	for(i in 1:length(periods)){
		j2 <- rich$period == periods[i] & j
		plot(rich$lon[j2], rich$lat[j2], col=rgb(colfun(pnorm01(rich$rich[j2], rich$rich[j])), maxColorValue=255), pch=16, cex=cexs, xlab='Longitude', ylab='Latitude', main=paste('NEUSSpring', periods[i]), cex.main=0.9)

		i <- pusplan$zone == 2
		i2 <- pusplan2$zone == 2
		points(pusplan$lon[i], pusplan$lat[i], pch=1, col='black', cex=cexs)
		points(pusplan2$lon[i2], pusplan2$lat[i2], pch=16, col='black', cex=0.3*cexs)
	}
	legend('bottomright', legend=round(seq(min(rich$rich[j]), max(rich$rich[j]), length.out=10),2), col=rgb(colfun(norm01(seq(min(rich$rich[j]), max(rich$rich[j]), length.out=10))), maxColorValue=255), pch=16, cex=0.7, title='Species', bty='n')


	dev.off()


# evaluate # targets met in each time period
	pinds <- presmap$region == myreg # trim to this region
	totals <- aggregate(list(total = presmap$pres[pinds]), by=list(sppocean=presmap$sppocean[pinds], period=presmap$period[pinds]), FUN=sum) # how many spp presences in each period
	totals <- totals[totals$sppocean %in% spps$sppocean,]
		length(unique(totals$sppocean))

	temp <- merge(presmap[pinds, ], pusplan[pusplan$zone==2,]) # only keep the conserved zones
	temp2 <- merge(presmap[pinds, ], pusplan2[pusplan2$zone==2,]) # only keep the conserved zones
	consabund <- aggregate(list(conserved = temp$pres), by=list(sppocean=temp$sppocean, period=temp$period), FUN=sum)
		dim(consabund)
	consabund2 <- aggregate(list(conserved = temp2$pres), by=list(sppocean=temp2$sppocean, period=temp2$period), FUN=sum)
		dim(consabund2)

	consabund2.1 <- merge(consabund, totals)
		dim(consabund2.1)
	consabund2.2 <- merge(consabund2, totals)
		dim(consabund2.2)

	consabund2.1$prop <- consabund2.1$conserved/consabund2.1$total
	consabund2.2$prop <- consabund2.2$conserved/consabund2.2$total

	goalsmet <- aggregate(list(nmet=consabund2.1$prop>=goal), by=list(period=consabund2.1$period), FUN=sum)
	goalsmet$mid <- sapply(strsplit(as.character(goalsmet$period), split='-'), FUN=function(x) mean(as.numeric(x)))
	goalsmet2 <- aggregate(list(nmet=consabund2.2$prop>=goal), by=list(period=consabund2.2$period), FUN=sum)
	goalsmet2$mid <- sapply(strsplit(as.character(goalsmet2$period), split='-'), FUN=function(x) mean(as.numeric(x)))

# plot goals met (solution #1)
	# quartz(width=4, height=4)
	cols = brewer.pal(4, 'Paired')
	pdf(width=4, height=4, file=paste('figures/MarZone_NEUSSpring_goalsmet_', runname, '.pdf', sep=''))

	plot(goalsmet$mid, goalsmet$nmet, xlab='Year', ylab='# Goals met', ylim=c(0, 125), type='o', pch=16, las=1, col=cols[2])
	
	dev.off()
	
# plot goals met (solution #1 and #2)
	# quartz(width=4, height=4)
	cols = brewer.pal(4, 'Paired')
	pdf(width=4, height=4, file=paste('figures/MarZone_NEUSSpring_goalsmet_', runname, '&', runname2, '.pdf', sep=''))

	plot(goalsmet$mid, goalsmet$nmet, xlab='Year', ylab='# Goals met', ylim=c(0, 125), type='o', pch=16, las=1, col=cols[2])
	points(goalsmet2$mid, goalsmet2$nmet, type='o', pch=16, col=cols[4])
	
	dev.off()
	
	

######################################################################
# Read in results and analyze against each climate model projection  #
# across each rcp                                                    #
######################################################################
runname <- 'conservationtest'
runname2 <- 'conservationtest2per'
goal <- 0.2
consplan <- read.csv(paste(marxfolder, 'output', runname, '/', runname, '_best.csv', sep=''))
consplan2 <- read.csv(paste(marxfolder, 'output', runname2, '/', runname2, '_best.csv', sep=''))
load(paste(marxfolder, 'input', runname, '/pus.Rdata', sep='')) # pus
load(paste(marxfolder, 'input', runname2, '/pus.Rdata', sep='')) # pus2
load(paste(marxfolder, 'input', runname, '/spps.Rdata', sep='')) # spps
load(paste(marxfolder, 'input', runname2, '/spps.Rdata', sep='')) # spps2

load('data/presmapbymod.RData') # loads presmap data.frame with presence/absence information from each model (slow to load)

# add zone to pus
pusplan <- merge(pus, consplan, by.x='id', by.y='planning_unit')
	dim(pus)
	dim(pusplan)

pusplan2 <- merge(pus2, consplan2, by.x='id', by.y='planning_unit')
	dim(pus2)
	dim(pusplan2)



# evaluate # targets met in each time period in each model
	pinds <- presmapbymod$region == 'NEFSC_NEUSSpring' # trim to this region
	totals <- aggregate(list(total = presmapbymod$pres[pinds]), by=list(sppocean=presmapbymod$sppocean[pinds], period=presmapbymod$period[pinds], model=presmapbymod$model[pinds]), FUN=sum) # how many presences in each period in each model
	totals <- totals[totals$sppocean %in% spps$sppocean,]
		length(unique(totals$sppocean))

	temp <- merge(presmapbymod[pinds, ], pusplan[pusplan$zone==2,]) # only keep the conserved zones
	temp2 <- merge(presmapbymod[pinds, ], pusplan2[pusplan2$zone==2,]) # only keep the conserved zones
	consabund <- aggregate(list(conserved = temp$pres), by=list(sppocean=temp$sppocean, period=temp$period, model=temp$model), FUN=sum)
		dim(consabund)
	consabund2 <- aggregate(list(conserved = temp2$pres), by=list(sppocean=temp2$sppocean, period=temp2$period, model=temp2$model), FUN=sum)
		dim(consabund2)

	intersect(names(consabund), names(totals)) # check before merging
	consabund2.1 <- merge(consabund, totals)
		dim(consabund2.1)
	consabund2.2 <- merge(consabund2, totals)
		dim(consabund2.2)

	consabund2.1$prop <- consabund2.1$conserved/consabund2.1$total
	consabund2.2$prop <- consabund2.2$conserved/consabund2.2$total

	goalsmet <- aggregate(list(nmet=consabund2.1$prop>=goal), by=list(period=consabund2.1$period, model=consabund2.1$model), FUN=sum, na.rm=TRUE) # remove NAs for species with 0 abundance
	goalsmet$mid <- sapply(strsplit(as.character(goalsmet$period), split='-'), FUN=function(x) mean(as.numeric(x)))
	goalsmet2 <- aggregate(list(nmet=consabund2.2$prop>=goal), by=list(period=consabund2.2$period, model=consabund2.2$model), FUN=sum, na.rm=TRUE)
	goalsmet2$mid <- sapply(strsplit(as.character(goalsmet2$period), split='-'), FUN=function(x) mean(as.numeric(x)))

	write.csv(consabund2.1, file=paste('output/consabund_', runname, '.csv', sep=''))
	write.csv(consabund2.2, file=paste('output/consabund_', runname2, '.csv', sep=''))
	write.csv(goalsmet, file=paste('output/goalsmet_', runname, '.csv', sep=''))
	write.csv(goalsmet2, file=paste('output/goalsmet_', runname2, '.csv', sep=''))

# compare goals met
	t.test(goalsmet$nmet[goalsmet$period=='2081-2100'], goalsmet2$nmet[goalsmet2$period=='2081-2100'])

# plot goals met (solution #1)
	# quartz(width=4, height=4)
	cols = brewer.pal(4, 'Paired')
	mods <- sort(unique(goalsmet$model))
	pdf(width=4, height=4, file=paste('figures/MarZone_NEUSSpring_goalsmetbymod_', runname, '.pdf', sep=''))

	inds <- goalsmet$model == 1
	plot(goalsmet$mid[inds], goalsmet$nmet[inds], xlab='Year', ylab='# Goals met', ylim=c(0, 125), type='o', pch=16, las=1, col=cols[1])
	for(i in 2:length(mods)){
		inds <- goalsmet$model == i
		points(goalsmet$mid[inds], goalsmet$nmet[inds], type='o', pch=16, col=cols[1])
	
	}
	ensmean <- aggregate(list(nmet=goalsmet$nmet), by=list(mid=goalsmet$mid), FUN=mean)
	lines(ensmean$mid, ensmean$nmet, col=cols[2], lwd=2)
	
	dev.off()
	
# plot goals met (solution #1 and #2)
	# quartz(width=4, height=4)
	cols = brewer.pal(4, 'Paired')
	mods <- sort(unique(goalsmet$model))
	pdf(width=4, height=4, file=paste('figures/MarZone_NEUSSpring_goalsmetbymod_', runname, '&', runname2, '.pdf', sep=''))

	inds <- goalsmet$model == 1
	plot(goalsmet$mid[inds], goalsmet$nmet[inds], xlab='Year', ylab='# Goals met', ylim=c(0, 125), type='o', pch=16, las=1, col=cols[1])
	for(i in 2:length(mods)){
		inds <- goalsmet$model == i
		points(goalsmet$mid[inds], goalsmet$nmet[inds], type='o', pch=16, col=cols[1])	
	}
	ensmean <- aggregate(list(nmet=goalsmet$nmet), by=list(mid=goalsmet$mid), FUN=mean)
	lines(ensmean$mid, ensmean$nmet, col=cols[2], lwd=2)

	inds <- goalsmet2$model == 1
	points(goalsmet2$mid[inds], goalsmet2$nmet[inds], type='o', pch=16, col=cols[3])
	for(i in 2:length(mods)){
		inds <- goalsmet2$model == i
		points(goalsmet2$mid[inds], goalsmet2$nmet[inds], type='o', pch=16, col=cols[3])
	
	}
	ensmean2 <- aggregate(list(nmet=goalsmet2$nmet), by=list(mid=goalsmet2$mid), FUN=mean)
	lines(ensmean2$mid, ensmean2$nmet, col=cols[4], lwd=2)

	
	dev.off()
	
