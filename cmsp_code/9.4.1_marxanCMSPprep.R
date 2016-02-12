# Set up a Marxan with Zones run for CMSP

## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	marxfolder <- '../MarZone_runs/'
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	.libPaths(new='~/R/x86_64-redhat-linux-gnu-library/3.1/') # so that it can find my old packages
	marxfolder <- 'MarZone_runs/'
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

# CMSP goals
consgoal <- 0.3 # proportion of presences to capture in conservation
energygoal <- 0.1 # proportion of NPV
fishgoal <- 0.7 # proportion of biomass

# name this run
runname1 <- 'cmsphistonly' # only historical period
runname2 <- 'cmsp2per' # 2 periods (historical and end of century)

# which region to run this for
myreg <- 'NEFSC_NEUSSpring'

# folders
inputfolder <- paste(marxfolder, 'input', runname1, sep='')
outputfolder <- paste(marxfolder, 'output', runname1, sep='')



####################
## helper functions
####################
require(RColorBrewer)

# normalize to 0-1
norm01 <- function(x){
	mn <- min(x)
	mx <- max(x)
	return((x-mn)/(mx-mn))
}

# parallel normalize to 0-1
# use y to guide the scale
pnorm01 <- function(x, y){
	mn <- min(y)
	mx <- max(y)
	return((x-mn)/(mx-mn))
}


############################################
## Set up a Marxan run just on 2006-2020
## use just one rcp
############################################
load(paste('data/presmap_', runtype, projtype, '_rcp', rcp, '.RData', sep='')) # loads presmap data.frame with presence/absence information
windnpv <- read.csv('cmsp_data/wind_npv.csv', row.names=1)
wavenpv <- read.csv('cmsp_data/wave_npv.csv', row.names=1)
fisheryspps <- read.csv('cmsp_data/fishery_spps.csv', row.names=1) # which spp to include in fishery goal in each region

# Create directory for input and output if missing
if(!dir.exists(inputfolder)){
	dir.create(inputfolder)
}
if(!dir.exists(outputfolder)){
	dir.create(outputfolder)
}

# Trim to just one region
presmap <- presmap[presmap$region==myreg,]
fisheryspps <- fisheryspps[fisheryspps$region==myreg,]

# Just one region for now: NEUS Spring spring
	# pu.dat
	# planning features are each 1/4 deg square
		pus <- presmap[,c('lat', 'lon')]
		pus <- pus[!duplicated(pus),]
			dim(pus)
		pus <- pus[order(pus$lat, pus$lon),]
		pus$id <- 1:nrow(pus)
	#	pus$dummycost <- runif(nrow(pus), 0, 1)
		pus$dummycost <- rep(0.1, nrow(pus)) # set the same cost in each planning unit. can add separate costs for each zone.
	
	#	pus <- pus[1:50,] # trim as a test
		
		pu.dat<-pus[,c('id', 'dummycost')] # the version to write out for Marxan

		save(pus, file=paste(inputfolder, '/pus.Rdata', sep=''))	
		write.csv(pu.dat, file=paste(inputfolder, '/pu.dat', sep=''), row.names=FALSE, quote=FALSE)
	
	# spec.dat (takes the place of feat.dat?)
	# goal for every planning feature present
	# not documented in 1.0.1 manual. copying format from unix example
		spps <- presmap[presmap$pres,c('sppocean')]
		spps <- spps[!duplicated(spps)]
			length(spps) # 102

		sppstokeep <- presmap[presmap$period=='2006-2020' & presmap$pres,c('lat', 'lon', 'sppocean', 'pres')]
			dim(sppstokeep)
		sppstokeep <- merge(sppstokeep, pus[,c('lat', 'lon', 'id')]) # add pu id (and trim to focal pus)
			names(sppstokeep)[names(sppstokeep)=='id'] <- 'pu'
			dim(sppstokeep)
	
			ngrid <- aggregate(list(ngrid=sppstokeep$pu), by=list(sppocean=sppstokeep$sppocean), FUN=function(x) length(unique(x)))
			sppstokeep <- merge(sppstokeep, ngrid)
			summary(sppstokeep$ngrid) # 1 to 543 (squalus acanthias is 543)

			nspps <- aggregate(list(nspp=sppstokeep$sppocean), by=list(pu=sppstokeep$pu), FUN=function(x) length(unique(x)))
			sppstokeep <- merge(sppstokeep, nspps)
			summary(sppstokeep$nspp) # 56 to 79

			sppstokeep <- sppstokeep[sppstokeep$ngrid> (nrow(pus)*0.05),] # trim to species found in at least 10% of grids

			length(unique(sppstokeep$sppocean)) # 67

		spps <- spps[spps %in% sppstokeep$sppocean]
			length(spps) # 67

		spps <- data.frame(id=1:length(spps), name=gsub(' |_', '', spps), sppocean=spps) #  fill spaces in species names.

	#	spps <- spps[1:20,] # trim as a test

		# add wind and wave energy feature
		spps <- rbind(spps, data.frame(id=max(spps$id)+1, name=c('energy'), sppocean=c(NA)))

		spps.dat <- spps[,c('id', 'name')]
	
		save(spps, file=paste(inputfolder, '/spps.Rdata', sep=''))		
		write.csv(spps.dat, file=paste(inputfolder, '/spec.dat', sep=''), row.names=FALSE, quote=FALSE)

	# puvsp.dat (takes the place of puvfeat.dat?)
	# planning features by planning units
	# not documented in 1.0.1 manual. copying format from unix example
		# Format species data
		puvsp <- presmap[presmap$period=='2006-2020',c('lat', 'lon', 'sppocean', 'wtcpue.proj', 'pres')]
			dim(puvsp)
		puvsp <- merge(puvsp, pus[,c('lat', 'lon', 'id')]) # add pu id (and trim to focal pus)
			dim(puvsp)
			names(puvsp)[names(puvsp)=='id'] <- 'pu'


		puvsp$name <- gsub(' |_', '', puvsp$sppocean) # trim out spaces on species names
		puvsp <- merge(puvsp, spps[,c('id', 'name')]) # merge in species IDs and trim to focal species
			dim(puvsp)
			names(puvsp)[names(puvsp)=='id'] <- 'species'
		puvsp$amount <- as.numeric(puvsp$wtcpue.proj) # use projected biomass as amount
		puvsp$amount[!puvsp$pres] <- 0 # set amount to zero where our cutoff says not present
	
		# Format wind and wave data
		puvenergy <- merge(windnpv, pus[,c('lat', 'lon', 'id')], all.y=TRUE)
		puvenergy <- merge(puvenergy, wavenpv)
			names(puvenergy)[names(puvenergy)=='id'] <- 'pu'
			head(puvenergy)
			dim(windnpv)
			dim(wavenpv)
			dim(puvenergy)
		puvenergy$wind_npv[puvenergy$wind_npv<0 | is.na(puvenergy$wind_npv)] <- 0 # set negative or NA NPV to 0
		puvenergy$wave_npv[puvenergy$wave_npv<0 | is.na(puvenergy$wave_npv)] <- 0
		puvenergy$amount <- puvenergy$wind_npv + puvenergy$wave_npv
		
		puvenergy$species <- spps$id[spps$name=='energy']
				
		length(unique(puvsp$pu))
		length(unique(puvenergy$pu))
		length(unique(puvsp$species))

		# Combine species, wind, and wave data for output
		puvsp.dat <- rbind(puvsp[,c('species', 'pu', 'amount')], puvenergy[,c('species', 'pu', 'amount')])
		puvsp.dat <- puvsp.dat[order(puvsp.dat$pu, puvsp.dat$species),]
		puvsp.dat <- puvsp.dat[puvsp.dat$amount>0,] # trim only to presences

			table(puvsp.dat$species) # make sure all species show up in some planning units
			table(puvsp.dat$pu) # make sure all planning units have some species

			sort(unique(table(puvsp.dat$pu, puvsp.dat$species))) # should be all 0s and 1s

		write.csv(puvsp.dat, file=paste(inputfolder, '/puvsp.dat', sep=''), row.names=FALSE, quote=FALSE)

	# zones
		zones <- data.frame(zoneid=1:4, zonename=c('available', 'conservation', 'fishery', 'energy'))
	
		write.csv(zones, file=paste(inputfolder, '/zones.dat', sep=''), row.names=FALSE, quote=FALSE)

	#costs
		costs <- data.frame(costid=1, costname='dummycost')
	
		write.csv(costs, file=paste(inputfolder, '/costs.dat', sep=''), row.names=FALSE, quote=FALSE)
	
	#zone cost
	#to adjust the importance of each cost in each zone
		zonecost <- expand.grid(list(zoneid=zones$zoneid, costid=costs$costid))
		zonecost$multiplier <- 0
		zonecost$multiplier[zonecost$zoneid %in% c(2,3,4) & zonecost$costid==1] <- 1 # set multiplier to non-zero for zones with a use
	
		zonecost <- zonecost[zonecost$multiplier>0,] # trim out zeros
	
		write.csv(zonecost, file=paste(inputfolder, '/zonecost.dat', sep=''), row.names=FALSE, quote=FALSE)
	
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
		zonetarget <- expand.grid(list(zoneid=zones$zoneid, speciesid=spps$id))
		zonetarget$target <- 0

		# set conservation zone target
		consinds <- zonetarget$zoneid==zones$zoneid[zones$zonename=='conservation']
		zonetarget$target[consinds] <- consgoal # XX proportion
		zonetarget$targettype[consinds] <- 3 # 3: proportion of total occurrences. 1: proportion of total amount

		# set fishing zone target
		fishinds <- zonetarget$zoneid==zones$zoneid[zones$zonename=='fishery'] & zonetarget$speciesid %in% spps$id[spps$sppocean %in% fisheryspps$projname]
		zonetarget$target[fishinds] <- fishgoal # XX proportion
		zonetarget$targettype[fishinds] <- 1 # 3: proportion of total occurrences. 1: proportion of total amount
	
		# set energy goal target
		energyinds <- zonetarget$zoneid==zones$zoneid[zones$zonename=='energy'] & zonetarget$speciesid %in% spps$id[spps$name == 'energy']
		zonetarget$target[energyinds] <- energygoal # XX proportion
		zonetarget$targettype[energyinds] <- 1 # 3: proportion of total occurrences. 1: proportion of total amount

		# format for output
		zonetarget.dat <- zonetarget[zonetarget$target>0,] # trim to only positive targets
		zonetarget.dat <- zonetarget.dat[order(zonetarget.dat$zoneid, zonetarget.dat$speciesid),] # order

		# write out	
		write.csv(zonetarget.dat, file=paste(inputfolder, '/zonetarget.dat', sep=''), row.names=FALSE, quote=FALSE)


	#input parameters
	input <- data.frame(BLM=0, PROP=0.5, RANDSEED=-1, NUMREPS=100, AVAILABLEZONE=1, NUMITNS='1000000', STARTTEMP=-1, NUMTEMP='10000', COSTTHRESH=0, THRESHPEN1=14, THRESHPEN2=10, INPUTDIR=paste('input', runname1, sep=''), PUNAME='pu.dat', SPECNAME='spec.dat', PUVSPRNAME='puvsp.dat', ZONESNAME='zones.dat', COSTSNAME='costs.dat', ZONECOSTNAME='zonecost.dat', ZONETARGETNAME='zonetarget.dat', SCENNAME=runname, SAVERUN=3, SAVEBEST=3, SAVESUMMARY=3, SAVESCEN=3, SAVETARGMET=3, SAVESUMSOLN=3, SAVEPENALTY=3, SAVELOG=3, OUTPUTDIR=paste('output', runname, sep=''),  RUNMODE=1, MISSLEVEL=1, ITIMPTYPE=0, HEURTYPE=-1, CLUMPTYPE=0, VERBOSITY=3, SAVESOLUTIONSMATRIX=3, SAVEANNEALINGTRACE=0, ANNEALINGTRACEROWS=1000)

	write.table(t(input), file=paste(marxfolder, 'input.dat', sep=''), row.names=TRUE, quote=FALSE, sep=' ', col.names=FALSE)
	write.table(t(input), file=paste(inputfolder, '/input.dat', sep=''), row.names=TRUE, quote=FALSE, sep=' ', col.names=FALSE)
	
	
# Go run MarZone!
# cd /Users/mpinsky/Documents/Rutgers/Range\ projections/MarZone_runs 
# ./MarZone_v201_Mac32 # will read in input.dat and write to the output folder

##################################################################
## Set up a Marxan run on 2006-2020 and ensemble mean 2081-2100
##################################################################
load('data/biomassavemap_testK6noSeas.RData') # loads biomassavemap data.frame
	dim(biomassavemap)
load('data/presmap.RData') # loads presmap data.frame with presence/absence information

goal <- 0.2 # proportion to capture in conservation
runname <- 'conservationtest2per'

# Just one region for now: NEUS Spring spring
	# pu.dat
	# planning features are each 1/4 deg square
	pus2 <- presmap[presmap$region=='NEFSC_NEUSSpring',c('lat', 'lon')]
	pus2 <- pus2[!duplicated(pus2),]
		dim(pus2)
	pus2 <- pus2[order(pus2$lat, pus2$lon),]
	pus2$id <- 1:nrow(pus2)
#	pus2$dummycost <- runif(nrow(pus2), 0, 1)
	pus2$dummycost <- rep(0.1, nrow(pus2))
	
#	pus2 <- pus2[1:50,] # trim as a test
		
	pu2.dat<-pus2[,c('id', 'dummycost')]
	
	save(pus2, file=paste(marxfolder, 'input', runname, '/pus.Rdata', sep=''))	
	write.csv(pu2.dat, file=paste(marxfolder, 'input', runname, '/pu.dat', sep=''), row.names=FALSE, quote=FALSE)
	
	# spec.dat (takes the place of feat.dat?)
	# goal for every species present
	# not documented in 1.0.1 manual. copying format from unix example
	inds <- presmap$region=='NEFSC_NEUSSpring' & presmap$period %in% c('2006-2020', '2081-2100')
	spps2 <- paste(presmap$sppocean[inds], presmap$period[inds])
	spps2 <- spps2[!duplicated(spps2)]
		length(spps2) # 246

	sppstokeep2 <- presmap[presmap$region=='NEFSC_NEUSSpring' & presmap$period %in% c('2006-2020', '2081-2100') & presmap$pres==TRUE,c('lat', 'lon', 'period', 'sppocean', 'pres')]
		dim(sppstokeep2)
	sppstokeep2 <- merge(sppstokeep2, pus2[,c('lat', 'lon', 'id')]) # add pu id (and trim to focal pus2)
		names(sppstokeep2)[names(sppstokeep2)=='id'] <- 'pu'
		dim(sppstokeep2)
	
		ngrid2 <- aggregate(list(ngrid=sppstokeep2$pu), by=list(sppocean=sppstokeep2$sppocean, period=sppstokeep2$period), FUN=function(x) length(unique(x)))
		sppstokeep2 <- merge(sppstokeep2, ngrid2)
		summary(sppstokeep2$ngrid) # 3 to 552

		nspps2 <- aggregate(list(nspp=sppstokeep2$sppocean), by=list(pu=sppstokeep2$pu, period=sppstokeep2$period), FUN=function(x) length(unique(x)))
		sppstokeep2 <- merge(sppstokeep2, nspps2)
		summary(sppstokeep2$nspp) # 4 to 85

		sppstokeep2 <- sppstokeep2[sppstokeep2$ngrid>10,] # trim to species at least minimally common

		length(unique(paste(sppstokeep2$sppocean, sppstokeep2$period))) # 242

	spps2 <- spps2[spps2 %in% paste(sppstokeep2$sppocean, sppstokeep2$period)]
		length(spps2) # 242

#	spps2 <- data.frame(id=1:length(spps2), prop=rep(goal, length(spps2)), spf=rep(10000, length(spps2)), name=gsub(' |_', '', spps2)) #  fill spaces in species names.
	spps2 <- data.frame(id=1:length(spps2), name=gsub(' |_|-', '', spps2), sppocean=spps2) #  remove spaces and dashes in species names.

#	spps2 <- spps2[1:20,] # trim as a test

	spps2.dat <- spps2[,c('id', 'name')]
	
	save(spps2, file=paste(marxfolder, 'input', runname, '/spps.Rdata', sep=''))		
	write.csv(spps2.dat, file=paste(marxfolder, 'input', runname, '/spec.dat', sep=''), row.names=FALSE, quote=FALSE)

	# puvsp.dat (takes the place of puvfeat.dat?)
	# species by planning units
	# not documented in 1.0.1 manual. copying format from unix example
	puvsp2 <- presmap[presmap$region=='NEFSC_NEUSSpring' & presmap$period %in% c('2006-2020', '2081-2100'),c('lat', 'lon', 'period', 'sppocean', 'pres')]
		dim(puvsp2)
	puvsp2 <- merge(puvsp2, pus2[,c('lat', 'lon', 'id')]) # add pu id (and trim to focal pus2)
		dim(puvsp2)
		names(puvsp2)[names(puvsp2)=='id'] <- 'pu'


	puvsp2$name <- gsub(' |_|-', '', paste(puvsp2$sppocean, puvsp2$period))
		dim(puvsp2)
	puvsp2 <- merge(puvsp2, spps2[,c('id', 'name')]) # add species id and trim to focal species
		dim(puvsp2)
		names(puvsp2)[names(puvsp2)=='id'] <- 'species'
	puvsp2$amount <- as.numeric(puvsp2$pres)
	
	length(unique(puvsp2$pu)) # 552
	length(unique(puvsp2$species)) # 242

	puvsp2.dat <- puvsp2[,c('species', 'pu', 'amount')]
		puvsp2.dat <- puvsp2.dat[order(puvsp2.dat$pu, puvsp2.dat$species),]
		puvsp2.dat <- puvsp2.dat[puvsp2.dat$amount>0,] # trim only to presences

		table(puvsp2.dat$species) # make sure all species show up in some planning units
			range(table(puvsp2.dat$species)) # make sure all species show up in some planning units
		table(puvsp2.dat$pu) # make sure all planning units have some species
			range(table(puvsp2.dat$pu)) # make sure all planning units have some species

#		table(puvsp2.dat$pu, puvsp2.dat$species) # giant matrix of who is where

	write.csv(puvsp2.dat, file=paste(marxfolder, 'input', runname, '/puvsp.dat', sep=''), row.names=FALSE, quote=FALSE)

	# zones
	zones2 <- data.frame(zoneid=1:2, zonename=c('available', 'conservation'))
	
	write.csv(zones2, file=paste(marxfolder, 'input', runname, '/zones.dat', sep=''), row.names=FALSE, quote=FALSE)

	#costs
	costs2 <- data.frame(costid=1, costname='dummycost')
	
	write.csv(costs2, file=paste(marxfolder, 'input', runname, '/costs.dat', sep=''), row.names=FALSE, quote=FALSE)
	
	#zone cost
	zonecost2 <- expand.grid(list(zoneid=zones$zoneid, costid=costs$costid))
	zonecost2$multiplier <- 0
	zonecost2$multiplier[zonecost2$zoneid==2 & zonecost2$costid==1] <- 1
	
	zonecost2 <- zonecost2[zonecost2$multiplier>0,]
	
	write.csv(zonecost2, file=paste(marxfolder, 'input', runname, '/zonecost.dat', sep=''), row.names=FALSE, quote=FALSE)
	
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
	zonetarget2 <- expand.grid(list(zoneid=zones$zoneid, speciesid=spps2$id))
	zonetarget2$target <- 0
#	zonetarget2$target[zonetarget2$zoneid==1 & zonetarget2$speciesid==1] <- 0.2
	zonetarget2$target[zonetarget2$zoneid==2] <- goal # XX proportion in conservation zones
	zonetarget2$targettype <- 1 # 3: proportion of total occurrences. 1: proportion of total amount
	
	zonetarget2.dat <- zonetarget2[zonetarget2$target>0,]
	zonetarget2.dat <- zonetarget2.dat[order(zonetarget2.dat$zoneid, zonetarget2.dat$speciesid),]
	
	write.csv(zonetarget2.dat, file=paste(marxfolder, 'input', runname, '/zonetarget.dat', sep=''), row.names=FALSE, quote=FALSE)


	#input parameters
	input2 <- data.frame(BLM=0, PROP=0.5, RANDSEED=-1, NUMREPS=100, AVAILABLEZONE=1, NUMITNS='1000000', STARTTEMP=-1, NUMTEMP='10000', COSTTHRESH=0, THRESHPEN1=14, THRESHPEN2=10, INPUTDIR=paste('input', runname, sep=''), PUNAME='pu.dat', SPECNAME='spec.dat', PUVSPRNAME='puvsp.dat', ZONESNAME='zones.dat', COSTSNAME='costs.dat', ZONECOSTNAME='zonecost.dat', ZONETARGETNAME='zonetarget.dat', SCENNAME=runname, SAVERUN=3, SAVEBEST=3, SAVESUMMARY=3, SAVESCEN=3, SAVETARGMET=3, SAVESUMSOLN=3, SAVEPENALTY=3, SAVELOG=3, OUTPUTDIR=paste('output', runname, sep=''),  RUNMODE=1, MISSLEVEL=1, ITIMPTYPE=0, HEURTYPE=-1, CLUMPTYPE=0, VERBOSITY=3, SAVESOLUTIONSMATRIX=3, SAVEANNEALINGTRACE=0, ANNEALINGTRACEROWS=1000)


	write.table(t(input2), file=paste(marxfolder, 'input.dat', sep=''), row.names=TRUE, quote=FALSE, sep=' ', col.names=FALSE)
	write.table(t(input2), file=paste(marxfolder, 'input', runname, '/input.dat', sep=''), row.names=TRUE, quote=FALSE, sep=' ', col.names=FALSE)
	
	
# Go run MarZone!	


