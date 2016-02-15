# Set up a Marxan with Zones run for CMSP

## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	marxfolder <- '../MarZone_runs/'
	presmapbymodfolder <- '../data/'
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	.libPaths(new='~/R/x86_64-redhat-linux-gnu-library/3.1/') # so that it can find my old packages
	marxfolder <- 'MarZone_runs/'
	presmapbymodfolder <- 'data/'
	}
# could add code for Lauren's working directory here


############
## Flags
############

## choose which runs to use
## runtype refers to how the Species Distribution Models (SDMs) were fit
## projtype refers to how the SDM projections were done
#runtype <- 'test'; projtype=''
#runtype <- ''; projtype=''`
#runtype <- 'testK6noSeas'; projtype='_xreg'
runtype <- 'fitallreg'; projtype='_xreg'

# choose the rcp (for runs using just one)
rcp <- 85
otherrcp <- 45

# CMSP goals
consgoal <- 0.1 # proportion of presences to capture in conservation
energygoal <- 0.1 # proportion of NPV
fishgoal <- 0.5 # proportion of biomass
conscolnm <- 'proppres'
fishcolnm <- 'proppres' # which column to use for fish goal (pres = occurrences, propwtcpue = biomass)

# choose region and name these runs
myreg <- 'NEFSC_NEUSSpring'; runname1 <- 'cmsphistonly'; runname2 <- 'cmsp2per'


# folders
inputfolder1 <- paste(marxfolder, runname1, '_input', sep='')
inputfolder2 <- paste(marxfolder, runname2, '_input', sep='')
outputfolder1 <- paste(marxfolder, runname1, '_output', sep='')
outputfolder2 <- paste(marxfolder, runname2, '_output', sep='')

######################
## Helper functions
######################
require(RColorBrewer)


#######################################################
# Read in results and simple prep
#######################################################
consplan1 <- read.csv(paste(outputfolder1, '/', runname1, '_best.csv', sep=''))
consplan2 <- read.csv(paste(outputfolder2, '/', runname2, '_best.csv', sep=''))

load(paste(inputfolder1, '/pus.Rdata', sep='')) # pus
load(paste(inputfolder2, '/pus.Rdata', sep='')) # pus2
load(paste(inputfolder1, '/spps.Rdata', sep='')) # spps
load(paste(inputfolder2, '/spps.Rdata', sep='')) # spps2
fisheryspps <- read.csv('cmsp_data/fishery_spps.csv', row.names=1) # which spp to include in fishery goal in each region

load(paste('data/rich_', runtype, projtype, '_rcp', rcp, '.RData', sep='')) # loads rich data.frame with presence/absence information
load(paste('data/presmap_', runtype, projtype, '_rcp', rcp, '.RData', sep='')) # loads presmap data.frame with presence/absence information

# load all model runs... very slow
load(paste(presmapbymodfolder, 'presmapbymod_', runtype, projtype, '_rcp', rcp, '.RData', sep='')) # loads presmap data.frame with presence/absence information from each model (slow to load).
	presmapbymod.1 <- presmapbymod[presmapbymod$region==myreg,]
	presmapbymod.1$rcp <- rcp
	dim(presmapbymod.1)
load(paste(presmapbymodfolder, 'presmapbymod_', runtype, projtype, '_rcp', otherrcp, '.RData', sep='')) # loads presmap data.frame with presence/absence information from each model (slow to load) (for the other rcp, the one not used for planning)
	presmapbymod$rcp <- otherrcp
	presmapbymod <- presmapbymod[presmapbymod$region==myreg,]
	dim(presmapbymod)
	presmapbymod <- rbind(presmapbymod.1, presmapbymod)
	dim(presmapbymod)
	rm(presmapbymod.1)

# add zone to pus
pusplan1 <- merge(pus, consplan1, by.x='id', by.y='planning_unit')
	dim(pus)
	dim(pusplan1)
	table(pusplan1$zone)

pusplan2 <- merge(pus2, consplan2, by.x='id', by.y='planning_unit')
	dim(pus2)
	dim(pusplan2)
	table(pusplan2$zone)


##############################
## Basic maps of solutions
##############################

# plot map of selected grids, on top of richness (consplan #1)
#	colfun <- colorRamp(rev(brewer.pal(11, 'Spectral')))
#	cexs = 1 # to adjust
#	periods <- sort(unique(rich$period))
#	# quartz(width=20, height=5)
#	pdf(width=20, height=5, file=paste('figures/MarZone_', myreg, '_on_richness_', runname1, '_', runtype, projtype, '_rcp', rcp, '.pdf', sep=''))
#	par(mfrow=c(1,length(periods)), mai=c(0.5,0.5,0.3, 0.1), las=1, mgp=c(2,1,0))
#	j <- rich$region == 'NEFSC_NEUSSpring'
#	for(i in 1:length(periods)){
#		j2 <- rich$period == periods[i] & j
#		plot(rich$lon[j2], rich$lat[j2], col=rgb(colfun(pnorm01(rich$rich[j2], rich$rich[j])), maxColorValue=255), pch=16, cex=cexs, xlab='Longitude', ylab='Latitude', main=paste(myreg, periods[i]), cex.main=0.9)
#
#		i <- pusplan1$zone == 2
#		points(pusplan1$lon[i], pusplan1$lat[i], pch=1, col='black', cex=cexs) # conservation
#
#		i <- pusplan1$zone == 3
#		points(pusplan1$lon[i], pusplan1$lat[i], pch=1, col='blue', cex=cexs) # fishery
#
#		i <- pusplan1$zone == 4
#		points(pusplan1$lon[i], pusplan1$lat[i], pch=1, col='purple', cex=cexs) # energy
#	}
#	legend('bottomright', legend=c(round(seq(min(rich$rich[j]), max(rich$rich[j]), length.out=10),2), 'Conservation', 'Fishery', 'Energy'), col=c(rgb(colfun(norm01(seq(min(rich$rich[j]), max(rich$rich[j]), length.out=10))), maxColorValue=255), 'black', 'blue', 'purple'), pch=c(rep(16,10), rep(1,3)), cex=0.7, title='Species', bty='n')
#
#
#	dev.off()

# plot map of selected grids, on top of richness (consplan #1 and #2)
#	colfun <- colorRamp(rev(brewer.pal(11, 'Spectral')))
#	cexs = 1 # to adjust
#	periods <- sort(unique(rich$period))
#	# quartz(width=20, height=10)
#	pdf(width=20, height=10, file=paste('figures/MarZone_', myreg, '_on_richness_', runname1, '&', runname2, '_', runtype, projtype, '_rcp', rcp, '.pdf', sep=''))
#	par(mfrow=c(2,length(periods)), mai=c(0.5,0.5,0.3, 0.1), las=1, mgp=c(2,1,0))
#	j <- rich$region == myreg
#	for(i in 1:length(periods)){ # for plan1
#		j2 <- rich$period == periods[i] & j
#		plot(rich$lon[j2], rich$lat[j2], col=rgb(colfun(pnorm01(rich$rich[j2], rich$rich[j])), maxColorValue=255), pch=16, cex=cexs, xlab='Longitude', ylab='Latitude', main=paste(myreg, periods[i]), cex.main=0.9)
#
#		i <- pusplan1$zone == 2
#		points(pusplan1$lon[i], pusplan1$lat[i], pch=1, col='black', cex=cexs) # conservation
#
#		i <- pusplan1$zone == 3
#		points(pusplan1$lon[i], pusplan1$lat[i], pch=1, col='blue', cex=cexs) # fishery
#
#		i <- pusplan1$zone == 4
#		points(pusplan1$lon[i], pusplan1$lat[i], pch=1, col='purple', cex=cexs) # energy
#	}
#	for(i in 1:length(periods)){ # for plan2
#		j2 <- rich$period == periods[i] & j
#		plot(rich$lon[j2], rich$lat[j2], col=rgb(colfun(pnorm01(rich$rich[j2], rich$rich[j])), maxColorValue=255), pch=16, cex=cexs, xlab='Longitude', ylab='Latitude', main=paste(myreg, periods[i]), cex.main=0.9)
#
#		i <- pusplan2$zone == 2
#		points(pusplan2$lon[i], pusplan2$lat[i], pch=1, col='black', cex=cexs) # conservation
#
#		i <- pusplan2$zone == 3
#		points(pusplan2$lon[i], pusplan2$lat[i], pch=1, col='blue', cex=cexs) # fishery
#
#		i <- pusplan2$zone == 4
#		points(pusplan2$lon[i], pusplan2$lat[i], pch=1, col='purple', cex=cexs) # energy
#	}
#	legend('bottomright', legend=c(round(seq(min(rich$rich[j]), max(rich$rich[j]), length.out=10),2), 'Conservation', 'Fishery', 'Energy'), col=c(rgb(colfun(norm01(seq(min(rich$rich[j]), max(rich$rich[j]), length.out=10))), maxColorValue=255), 'black', 'blue', 'purple'), pch=c(rep(16,10), rep(1,3)), cex=0.7, title='Species', bty='n')
#
#
#	dev.off()
	
	
# plot map of selected grids (consplan #1 and #2)
#	cexs = 0.5 # to adjust
#	periods <- sort(unique(rich$period))
#	# quartz(width=5, height=3)
#	pdf(width=5, height=3, file=paste('figures/MarZone_', myreg, '_', runname1, '&', runname2, '_', runtype, projtype, '_rcp', rcp, '.pdf', sep=''))
#	par(mfrow=c(1,2), mai=c(0.5,0.5,0.3, 0.1), las=1, mgp=c(2,1,0))
#	j <- rich$region == myreg
#	i <- pusplan1$zone == 2
#	plot(pusplan1$lon[i], pusplan1$lat[i], pch=16, col='black', cex=cexs, main='Historical only', xlab='', ylab='') # conservation
#	i <- pusplan1$zone == 3
#	points(pusplan1$lon[i], pusplan1$lat[i], pch=16, col='blue', cex=cexs) # fishery
#	i <- pusplan1$zone == 4
#	points(pusplan1$lon[i], pusplan1$lat[i], pch=16, col='red', cex=cexs) # energy
#
#	i <- pusplan2$zone == 2
#	plot(pusplan2$lon[i], pusplan2$lat[i], pch=16, col='black', cex=cexs, main='Two periods', xlab='', ylab='') # conservation
#	i <- pusplan2$zone == 3
#	points(pusplan2$lon[i], pusplan2$lat[i], pch=16, col='blue', cex=cexs) # fishery
#	i <- pusplan2$zone == 4
#	points(pusplan2$lon[i], pusplan2$lat[i], pch=16, col='red', cex=cexs) # energy
#	
#	legend('bottomright', legend=c('Conservation', 'Fishery', 'Energy'), col=c('black', 'blue', 'purple'), pch=rep(16,3), cex=0.7, title='Zones', bty='n')
#
#
#	dev.off()


#######################################################
# Analyze plan against the ensemble mean
#######################################################
pinds <- presmap$region == myreg & presmap$pres # trim to presences in this region
finds <- fisheryspps$region == myreg # trim to this region

# evaluate # targets met in each time period (only biogical targets)

	# calculate total presences and wtcpue by species and time-period
	totals <- expand.grid(sppocean=spps$sppocean[spps$name != 'energy'], period=sort(unique(presmap$period))) # grid of all spp and time periods
	totcalcs <- aggregate(list(totalpres = presmap$pres[pinds], totalwtcpue=presmap$wtcpue.proj[pinds]), by=list(sppocean=presmap$sppocean[pinds], period=presmap$period[pinds]), FUN=sum) # how many spp presences and amount in each period
	totals2 <- merge(totals, totcalcs, all.x=TRUE) # make sure all spp and time period represented
	totals2$totalpres[is.na(totals2$totalpres)] <- 0
	totals2$totalwtcpue[is.na(totals2$totalwtcpue)] <- 0
	totals2$zone <- 2 # for conservation

	# add totals for the fishery zones
	temp <- totals2
	temp$zone <- 3
	temp <- temp[temp$sppocean %in% fisheryspps$projname[finds],]
		nrow(temp)/5 # 8 species
	totals2 <- rbind(totals2, temp) # for conservation

		length(unique(totals2$sppocean)) # number of species
		length(unique(paste(totals2$sppocean, totals2$zone))) # number of goals: 75
		nrow(totals2)

	# merge in plan zones
	temp1 <- merge(presmap[pinds, ], pusplan1[pusplan1$zone %in% c(2,3),]) # only keep the conserved & fishery zones for goals
	temp2 <- merge(presmap[pinds, ], pusplan2[pusplan2$zone %in% c(2,3),])
	abundbyzone1 <- aggregate(list(npres = temp1$pres, sumwtcpue = temp1$wtcpue.proj), by=list(sppocean=temp1$sppocean, period=temp1$period, zone=temp1$zone), FUN=sum)
		dim(abundbyzone1)
	abundbyzone2 <- aggregate(list(npres = temp2$pres, sumwtcpue = temp2$wtcpue.proj), by=list(sppocean=temp2$sppocean, period=temp2$period, zone=temp2$zone), FUN=sum)
		dim(abundbyzone2)

	# add totals for pres and wtcpue across all zones
	# intersect(names(abundbyzone1), names(totals2))
	abundbyzone1.2 <- merge(abundbyzone1, totals2, all.y=TRUE)
		abundbyzone1.2$npres[is.na(abundbyzone1.2$npres)] <- 0
		abundbyzone1.2$sumwtcpue[is.na(abundbyzone1.2$sumwtcpue)] <- 0
		dim(abundbyzone1.2) # 375
		length(unique(abundbyzone1.2$sppocean)) # 67 species
		sort(table(as.character(abundbyzone1.2$sppocean))) # entries per species, from low to high. 8 species appear in both conservation and fishery zones.
	abundbyzone2.2 <- merge(abundbyzone2, totals2, all.y=TRUE)
		abundbyzone2.2$npres[is.na(abundbyzone2.2$npres)] <- 0
		abundbyzone2.2$sumwtcpue[is.na(abundbyzone2.2$sumwtcpue)] <- 0
		dim(abundbyzone2.2) # 375
		length(unique(abundbyzone2.2$sppocean)) # 67 species
		sort(table(as.character(abundbyzone2.2$sppocean))) # entries per species, from low to high. 8 species appear in both conservation and fishery zones.

	# calculate proportion of presences and wtcpue
	abundbyzone1.2$proppres <- abundbyzone1.2$npres/abundbyzone1.2$totalpres
	abundbyzone2.2$proppres <- abundbyzone2.2$npres/abundbyzone2.2$totalpres
	abundbyzone1.2$propwtcpue <- abundbyzone1.2$sumwtcpue/abundbyzone1.2$totalwtcpue
	abundbyzone2.2$propwtcpue <- abundbyzone2.2$sumwtcpue/abundbyzone2.2$totalwtcpue

	# force 0/0 to 1 so that it counts as a goal met
	abundbyzone1.2$proppres[abundbyzone1.2$npres==0 & abundbyzone1.2$totalpres==0] <- 1
	abundbyzone2.2$proppres[abundbyzone2.2$npres==0 & abundbyzone2.2$totalpres==0] <- 1
	abundbyzone1.2$propwtcpue[abundbyzone1.2$sumwtcpue==0 & abundbyzone1.2$totalwtcpue==0] <- 1
	abundbyzone2.2$propwtcpue[abundbyzone2.2$sumwtcpue==0 & abundbyzone2.2$totalwtcpue==0] <- 1
	
	# mark where goals met for conservation or fishery
	abundbyzone1.2$metgoal <- FALSE
	abundbyzone2.2$metgoal <- FALSE
	abundbyzone1.2$metgoal[abundbyzone1.2$zone==2 & abundbyzone1.2[[conscolnm]] >= consgoal] <- TRUE
	abundbyzone2.2$metgoal[abundbyzone2.2$zone==2 & abundbyzone2.2[[conscolnm]] >= consgoal] <- TRUE
	abundbyzone1.2$metgoal[abundbyzone1.2$zone==3 & abundbyzone1.2[[fishcolnm]] >= fishgoal] <- TRUE
	abundbyzone2.2$metgoal[abundbyzone2.2$zone==3 & abundbyzone2.2[[fishcolnm]] >= fishgoal] <- TRUE

	# calculate number goals met in each timeperiod
	goalsmet1 <- aggregate(list(nmet=abundbyzone1.2$metgoal), by=list(period=abundbyzone1.2$period), FUN=sum)
	goalsmet1$mid <- sapply(strsplit(as.character(goalsmet1$period), split='-'), FUN=function(x) mean(as.numeric(x)))
	goalsmet2 <- aggregate(list(nmet=abundbyzone2.2$metgoal), by=list(period=abundbyzone2.2$period), FUN=sum)
	goalsmet2$mid <- sapply(strsplit(as.character(goalsmet2$period), split='-'), FUN=function(x) mean(as.numeric(x)))

# write out
	write.csv(abundbyzone1.2, file=paste('output/abundbyzone_', runtype, projtype, '_', runname1, '.csv', sep=''))
	write.csv(abundbyzone2.2, file=paste('output/abundbyzone_', runtype, projtype, '_', runname2, '.csv', sep=''))
	write.csv(goalsmet1, file=paste('output/goalsmet_', runtype, projtype, '_', runname1, '.csv', sep=''))
	write.csv(goalsmet2, file=paste('output/goalsmet_', runtype, projtype, '_', runname2, '.csv', sep=''))


# plot goals met (solution #1)
	# quartz(width=4, height=4)
#	cols = brewer.pal(4, 'Paired')
#	ylims <- c(0, max(goalsmet1$nmet))
#	pdf(width=4, height=4, file=paste('figures/MarZone_NEUSSpring_goalsmet_', runname, '.pdf', sep=''))
#
#	plot(goalsmet1$mid, goalsmet1$nmet, xlab='Year', ylab='# Goals met', ylim=ylims, type='o', pch=16, las=1, col=cols[2])
#	
#	dev.off()
	
# plot goals met (solution #1 and #2)
	cols = brewer.pal(4, 'Paired')
	ylims <- c(0, max(c(max(goalsmet1$nmet), max(goalsmet2$nmet))))
	# quartz(width=4, height=4)
	pdf(width=4, height=4, file=paste('figures/MarZone_', myreg, '_goalsmet_', runname1, '&', runname2, '_', runtype, projtype, '_rcp', rcp, '.pdf', sep=''))

	plot(goalsmet1$mid, goalsmet1$nmet, xlab='Year', ylab='# Goals met', ylim=ylims, type='o', pch=16, las=1, col=cols[2])
	points(goalsmet2$mid, goalsmet2$nmet, type='o', pch=16, col=cols[4])
	
	dev.off()
	
	

######################################################################
# Analyze against each climate model projection                      #
# across each rcp                                                    #
######################################################################
# better to run this on Amphiprion: takes 30 GB memory

pinds <- presmapbymod$region == myreg # trim to this region
finds <- fisheryspps$region == myreg # trim to this region


# evaluate # targets met in each time period in each model
	# calculate total presences and wtcpue by species and time-period
	totals <- expand.grid(sppocean=spps$sppocean[spps$name != 'energy'], period=sort(unique(presmapbymod$period)), model=sort(unique(presmapbymod$model)), rcp=c(rcp, otherrcp)) # grid of all spp and time periods and models and rcps
	totcalcs <- aggregate(list(totalpres = presmapbymod$pres[pinds], totalwtcpue=presmapbymod$wtcpue.proj[pinds]), by=list(sppocean=presmapbymod$sppocean[pinds], period=presmapbymod$period[pinds], model=presmapbymod$model[pinds], rcp=presmapbymod$rcp[pinds]), FUN=sum) # how many spp presences and amount in each period in each model for focal rcp

	totals2 <- merge(totals, totcalcs, all.x=TRUE) # make sure all spp and time periods and models and rcps represented
	totals2$totalpres[is.na(totals2$totalpres)] <- 0
	totals2$totalwtcpue[is.na(totals2$totalwtcpue)] <- 0
	totals2$zone <- 2 # for conservation

	# add totals for the fishery zones
	temp <- totals2
	temp$zone <- 3
	temp <- temp[temp$sppocean %in% fisheryspps$projname[finds],]
		nrow(temp)/5/13/2 # 8 species: good
	totals2 <- rbind(totals2, temp) # for fishery

		length(unique(totals2$sppocean)) # number of species: 67
		numgoals <- length(unique(paste(totals2$sppocean, totals2$zone))) # number of goals: 75
			numgoals
		nrow(totals2)

	# merge in plan zones
	temp1 <- merge(presmapbymod[pinds, ], pusplan1[pusplan1$zone %in% c(2,3),]) # only keep the conserved & fishery zones for goals
	temp2 <- merge(presmapbymod[pinds, ], pusplan2[pusplan2$zone %in% c(2,3),])
	abundbyzone1 <- aggregate(list(npres = temp1$pres, sumwtcpue = temp1$wtcpue.proj), by=list(sppocean=temp1$sppocean, period=temp1$period, zone=temp1$zone, model=temp1$model, rcp=temp1$rcp), FUN=sum)
		dim(abundbyzone1) # 67080
	abundbyzone2 <- aggregate(list(npres = temp2$pres, sumwtcpue = temp2$wtcpue.proj), by=list(sppocean=temp2$sppocean, period=temp2$period, zone=temp2$zone, model=temp2$model, rcp=temp2$rcp), FUN=sum)
		dim(abundbyzone2) # 67080

	# add totals for pres and wtcpue across all zones
	# intersect(names(abundbyzone1), names(totals2))
	abundbyzonebymod1.2 <- merge(abundbyzone1, totals2, all.y=TRUE)
		abundbyzonebymod1.2$npres[is.na(abundbyzonebymod1.2$npres)] <- 0
		abundbyzonebymod1.2$sumwtcpue[is.na(abundbyzonebymod1.2$sumwtcpue)] <- 0
		dim(abundbyzonebymod1.2) # 9750
		length(unique(abundbyzonebymod1.2$sppocean)) # 67 species
		sort(table(as.character(abundbyzonebymod1.2$sppocean))) # entries per species, from low to high. 8 species appear in both conservation and fishery zones.
	abundbyzonebymod2.2 <- merge(abundbyzone2, totals2, all.y=TRUE)
		abundbyzonebymod2.2$npres[is.na(abundbyzonebymod2.2$npres)] <- 0
		abundbyzonebymod2.2$sumwtcpue[is.na(abundbyzonebymod2.2$sumwtcpue)] <- 0
		dim(abundbyzonebymod2.2) # 9750
		length(unique(abundbyzonebymod2.2$sppocean)) # 67 species
		sort(table(as.character(abundbyzonebymod2.2$sppocean))) # entries per species, from low to high. 8 species appear in both conservation and fishery zones.

	# calculate proportion of presences and wtcpue
	abundbyzonebymod1.2$proppres <- abundbyzonebymod1.2$npres/abundbyzonebymod1.2$totalpres
	abundbyzonebymod2.2$proppres <- abundbyzonebymod2.2$npres/abundbyzonebymod2.2$totalpres
	abundbyzonebymod1.2$propwtcpue <- abundbyzonebymod1.2$sumwtcpue/abundbyzonebymod1.2$totalwtcpue
	abundbyzonebymod2.2$propwtcpue <- abundbyzonebymod2.2$sumwtcpue/abundbyzonebymod2.2$totalwtcpue

	# force 0/0 to 1 so that it counts as a goal met
	abundbyzonebymod1.2$proppres[abundbyzonebymod1.2$npres==0 & abundbyzonebymod1.2$totalpres==0] <- 1
	abundbyzonebymod2.2$proppres[abundbyzonebymod2.2$npres==0 & abundbyzonebymod2.2$totalpres==0] <- 1
	abundbyzonebymod1.2$propwtcpue[abundbyzonebymod1.2$sumwtcpue==0 & abundbyzonebymod1.2$totalwtcpue==0] <- 1
	abundbyzonebymod2.2$propwtcpue[abundbyzonebymod2.2$sumwtcpue==0 & abundbyzonebymod2.2$totalwtcpue==0] <- 1
	
	# mark where goals met for conservation or fishery
	abundbyzonebymod1.2$metgoal <- FALSE
	abundbyzonebymod2.2$metgoal <- FALSE
	abundbyzonebymod1.2$metgoal[abundbyzonebymod1.2$zone==2 & abundbyzonebymod1.2[[conscolnm]] >= consgoal] <- TRUE
	abundbyzonebymod2.2$metgoal[abundbyzonebymod2.2$zone==2 & abundbyzonebymod2.2[[conscolnm]] >= consgoal] <- TRUE
	abundbyzonebymod1.2$metgoal[abundbyzonebymod1.2$zone==3 & abundbyzonebymod1.2[[fishcolnm]] >= fishgoal] <- TRUE
	abundbyzonebymod2.2$metgoal[abundbyzonebymod2.2$zone==3 & abundbyzonebymod2.2[[fishcolnm]] >= fishgoal] <- TRUE

	# calculate number goals met in each timeperiod
	goalsmetbymod1 <- aggregate(list(nmet=abundbyzonebymod1.2$metgoal), by=list(period=abundbyzonebymod1.2$period, model=abundbyzonebymod1.2$model, rcp=abundbyzonebymod1.2$rcp), FUN=sum)
	goalsmetbymod1$mid <- sapply(strsplit(as.character(goalsmetbymod1$period), split='-'), FUN=function(x) mean(as.numeric(x)))
	goalsmetbymod1$pmet <- goalsmetbymod1$nmet/numgoals

	goalsmetbymod2 <- aggregate(list(nmet=abundbyzonebymod2.2$metgoal), by=list(period=abundbyzonebymod2.2$period, model=abundbyzonebymod1.2$model, rcp=abundbyzonebymod1.2$rcp), FUN=sum)
	goalsmetbymod2$mid <- sapply(strsplit(as.character(goalsmetbymod2$period), split='-'), FUN=function(x) mean(as.numeric(x)))
	goalsmetbymod2$pmet <- goalsmetbymod2$nmet/numgoals

	# write out
	write.csv(abundbyzonebymod1.2, file=paste('output/abundbyzonebymod_', runtype, projtype, '_', runname1, '.csv', sep=''))
	write.csv(abundbyzonebymod2.2, file=paste('output/abundbyzonebymod_', runtype, projtype, '_', runname2, '.csv', sep=''))
	write.csv(goalsmetbymod1, file=paste('output/goalsmetbymod_', runtype, projtype, '_', runname1, '.csv', sep=''))
	write.csv(goalsmetbymod2, file=paste('output/goalsmetbymod_', runtype, projtype, '_', runname2, '.csv', sep=''))

# compare goals met
	t.test(goalsmetbymod1$nmet[goalsmetbymod1$period=='2081-2100'], goalsmetbymod2$nmet[goalsmetbymod2$period=='2081-2100'])


# plot goals met (solution #1)
	# quartz(width=4, height=4)
#	cols = brewer.pal(4, 'Paired')
#	mods <- sort(unique(goalsmetbymod1$model))
#	pdf(width=4, height=4, file=paste('figures/MarZone_NEUSSpring_goalsmetbymod_', runname1, '.pdf', sep=''))
#
#	inds <- goalsmetbymod1$model == 1
#	plot(goalsmetbymod1$mid[inds], goalsmetbymod1$nmet[inds], xlab='Year', ylab='# Goals met', ylim=c(0, 125), type='o', pch=16, las=1, col=cols[1])
#	for(i in 2:length(mods)){
#		inds <- goalsmetbymod1$model == i
#		points(goalsmetbymod1$mid[inds], goalsmetbymod1$nmet[inds], type='o', pch=16, col=cols[1])
#	
#	}
#	ensmean <- aggregate(list(nmet=goalsmetbymod1$nmet), by=list(mid=goalsmetbymod1$mid), FUN=mean)
#	lines(ensmean$mid, ensmean$nmet, col=cols[2], lwd=2)
#	
#	dev.off()
	
# plot goals met (solution #1 and #2)
	# quartz(width=4, height=4)
	colmat <- t(col2rgb(brewer.pal(4, 'Paired')))
	cols <- rgb(red=colmat[,1], green=colmat[,2], blue=colmat[,3], alpha=c(80,255,80,255), maxColorValue=255)
	mods <- sort(unique(goalsmetbymod1$model))
	rcps <- sort(unique(goalsmetbymod1$rcp))
	pdf(width=4, height=4, file=paste('figures/MarZone_', myreg, '_goalsmetbymod_', runtype, projtype, '_', runname1, '&', runname2, '.pdf', sep=''))

	inds <- goalsmetbymod1$model == mods[1] & goalsmetbymod1$rcp == rcps[1]
	plot(goalsmetbymod1$mid[inds], goalsmetbymod1$pmet[inds], xlab='Year', ylab='# Goals met', ylim=c(0, 1), type='l', pch=16, las=1, col=cols[1])
	for(i in 1:length(mods)){
		for(j in 1:length(rcps)){
			if(!(i==1 & j==1)){
				inds <- goalsmetbymod1$model == mods[i] & goalsmetbymod1$rcp == rcps[j]
				points(goalsmetbymod1$mid[inds], goalsmetbymod1$pmet[inds], type='l', pch=16, col=cols[1])	
			}
		}
	}
	ensmean <- aggregate(list(nmet=goalsmetbymod1$nmet, pmet=goalsmetbymod1$pmet), by=list(mid=goalsmetbymod1$mid), FUN=mean)
	lines(ensmean$mid, ensmean$pmet, col=cols[2], lwd=2)

	for(i in 1:length(mods)){
		for(j in 1:length(rcps)){
			inds <- goalsmetbymod2$model == mods[i] & goalsmetbymod2$rcp == rcps[j]
			points(goalsmetbymod2$mid[inds], goalsmetbymod2$pmet[inds], type='l', pch=16, col=cols[3])
		}	
	}
	ensmean2 <- aggregate(list(nmet=goalsmetbymod2$nmet, pmet=goalsmetbymod2$pmet), by=list(mid=goalsmetbymod2$mid), FUN=mean)
	lines(ensmean2$mid, ensmean2$pmet, col=cols[4], lwd=2)

	
	dev.off()
	
