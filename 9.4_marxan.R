# Set up a Marxan with Zones run

## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	marxfolder <- '../MarZone_runs/'
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	}
# could add code for Lauren's working directory here

############################################
## Set up a Marxan run
############################################
load('data/biomassavemap_testK6noSeas.RData') # loads biomassavemap data.frame
	dim(biomassavemap)
load('data/presmap.RData') # loads presmap data.frame with presence/absence information

# Just one region for now: NEUS Spring spring
	# pu.dat
	# planning features are each 1/4 deg square
	pus <- presmap[presmap$region=='NEFSC_NEUS',c('lat', 'lon')]
	pus <- pus[!duplicated(pus),]
		dim(pus)
	pus <- pus[order(pus$lat, pus$lon),]
	pus$id <- 1:nrow(pus)
	pus$dummycost <- 1
		
	pu.dat<-pus[,c('id', 'dummycost')]
	
	write.csv(pu.dat, file=paste(marxfolder, 'data/pu.dat', sep=''), row.names=FALSE)
	
	# feat.dat
	# goal for every species present
	feats <- presmap[presmap$region=='NEFSC_NEUS',c('sppocean')]
	feats <- feats[!duplicated(feats)]
		length(feats) # 123
	feats <- data.frame(id=1:length(feats), name=gsub(' ', '_', feats)) #  fill spaces in species names.
	
	write.csv(feats, file=paste(marxfolder, 'data/feat.dat', sep=''), row.names=FALSE)

	# puvfeat.dat
	# features by planning units
	puvfeat <- presmap[presmap$region=='NEFSC_NEUS' & presmap$period=='2006-2020',c('lat', 'lon', 'sppocean', 'pres')]
		dim(puvfeat)
	puvfeat <- merge(puvfeat, pus[,c('lat', 'lon', 'id')]) # add pu id
		dim(puvfeat)
		names(puvfeat)[names(puvfeat)=='id'] <- 'puid'
	puvfeat$name <- gsub(' ', '_', puvfeat$sppocean)
	puvfeat <- merge(puvfeat, feats[,c('id', 'name')])
		dim(puvfeat)
		names(puvfeat)[names(puvfeat)=='id'] <- 'featureid'
	puvfeat$amount <- as.numeric(puvfeat$pres)
	
	puvfeat.dat <- puvfeat[,c('featureid', 'puid', 'amount')]
	write.csv(puvfeat.dat, file=paste(marxfolder, 'data/puvfeat.dat', sep=''), row.names=FALSE)

	# zones
	zones <- data.frame(zoneid=1:2, zonename=c('open', 'conservation'))
	
	write.csv(zones, file=paste(marxfolder, 'data/zones.dat', sep=''), row.names=FALSE)

	#costs
	costs <- data.frame(costid=1, costname='dummycost')
	
	write.csv(costs, file=paste(marxfolder, 'data/costs.dat', sep=''), row.names=FALSE)
	
	#zone cost
	zonecost <- expand.grid(list(zoneid=zones$zoneid, costid=costs$costid))
	zonecost$multiplier <- 1
	
	write.csv(zonecost, file=paste(marxfolder, 'data/zonecost.dat', sep=''), row.names=FALSE)
	
	#boundary length
	# optional
	
	#zone boundary cost
	# optional
	
	#planning unit zone
	#optional
	
	#planning unit lock
	#optional
	
	#zone target
	# set zone-specific targets
	zonetarget <- expand.grid(list(zoneid=zones$zoneid, featureid=feats$id))
	zonetarget$target <- 0
	zonetarget$target[zonetarget$zoneid==2] <- 0.1 # 10% in conservation zones
	zonetarget$targettype <- 3 # proportion of total occurrences
	
	write.csv(zonetarget, file=paste(marxfolder, 'data/zonetarget.dat', sep=''), row.names=FALSE)


	#input parameters
	input <- data.frame(parameter=c('BLM', 'PROP', 'RANDSEED', 'NUMREPS', 'AVAILABLEZONE', 'NUMITNS', 'VERBOSITY'), value=c(0,0.5,-1, 100, 1, 1000000, 0))

	write.csv(input, file=paste(marxfolder, 'data/input.dat', sep=''), row.names=FALSE)