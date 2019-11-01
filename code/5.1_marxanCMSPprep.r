# Set up a Marxan with Zones run for CMSP


############
## Flags
############

# choose the rcp (for runs using just one)
rcp <- 85

# CMSP goals
consgoal <- 0.1 # proportion of presences to capture in conservation
energygoal <- 0.2 # proportion of NPV
fishgoal <- 0.5 # proportion of biomass
cost <- 0.01 # basic cost of including each planning unit in a given zone
zonecosts <- c(0.1, 1, 1, 1) # cost multipliers for available, conservation, fishery, and energy zones
fpf <- 10 # feature penalty factor

# poccur threshold: how high does the probability of occurrence in the projections need to be to consider the species "present"?
poccurthresh <- 0.1

# choose region and name these runs
myreg <- 'ebs'; runname1 <- 'hist_ebs'; runname2 <- '2per_ebs'
#myreg <- 'goa'; runname1 <- 'hist_goa'; runname2 <- '2per_goa'
#myreg <- 'bc'; runname1 <- 'hist_bc'; runname2 <- '2per_bc'
#myreg <- 'wc'; runname1 <- 'hist_wc'; runname2 <- '2per_wc'
#myreg <- 'gmex'; runname1 <- 'hist_gmex'; runname2 <- '2per_gmex'
#myreg <- 'seus'; runname1 <- 'hist_seus'; runname2 <- '2per_seus'
#myreg <- 'neus'; runname1 <- 'hist_neus'; runname2 <- '2per_neus'
#myreg <- 'maritime'; runname1 <- 'hist_maritime'; runname2 <- '2per_maritime'
#myreg <- 'newf'; runname1 <- 'hist_newf'; runname2 <- '2per_newf'

# which time periods to use in the multi-period planning
planningperiods <- c('2007-2020', '2081-2100')

# folders
marxfolder <- 'marzone_runs/'

######################
# Functions
######################
require(data.table)
library(prioritizr) # only runs in R 3.5.3 for now (Gurobi 8.1.1)


#####################
## Load data
#####################
# define ocean
if (myreg %in% c('wc', 'bc', 'goa', 'ebs')) myocean <- 'Pac'
if (myreg %in% c('gmex', 'seus', 'neus', 'maritime', 'newf')) myocean <- 'Atl'
if (!(myreg %in% c('gmex', 'seus', 'neus', 'maritime', 'newf', 'wc', 'bc', 'goa', 'ebs'))) stop('myreg not an existing region')

# loads presence/absence and biomass data
presmap <- fread(cmd = paste0('gunzip -c temp/presmap_', myocean, '_rcp', rcp, '.csv.gz'), drop = 1) 
biomassmap <- fread(cmd = paste0('gunzip -c temp/biomassmap_', myocean, '_rcp', rcp, '.csv.gz'), drop = 1) 

# load NatCap calculations
windnpv <- fread(cmd = 'gunzip -c output/wind_npv.csv.gz', drop = 1)
wavenpv <- fread(cmd = 'gunzip -c output/wave_npv.csv.gz', drop = 1)

setnames(windnpv, c('lat', 'lon'), c('latgrid', 'longrid'))
setnames(wavenpv, c('lat', 'lon'), c('latgrid', 'longrid'))

# definition of fishery species by region
fisheryspps <- fread('output/fishery_spps.csv', drop = 1) # which spp to include in fishery goal in each region

# region definitions
regiongrid <- fread(cmd = 'gunzip -c output/region_grid.csv.gz', drop = 1)

################
## Set up data
################

# Fix lon in regiongrid to match presmap (-360 to 0)
regiongrid[longrid > 0, longrid := longrid - 360] 

# Add region information to presmap
setkey(presmap, latgrid, longrid)
setkey(regiongrid, latgrid, longrid)
presmap <- merge(presmap, regiongrid[, .(latgrid, longrid, region)], all.x = TRUE) # add region information
    presmap[is.na(region) & !duplicated(presmap[,.(latgrid, longrid)]), .N] # 0 missing region: good!
    # presmap[is.na(region) & !duplicated(presmap[,.(latgrid, longrid)]), ]
    # presmap[is.na(region) & !duplicated(presmap[,.(latgrid, longrid)]), plot(longrid, latgrid)]
    
# Add region information to biomassmap
setkey(biomassmap, latgrid, longrid)
setkey(regiongrid, latgrid, longrid)
biomassmap <- merge(biomassmap, regiongrid[, .(latgrid, longrid, region)], all.x = TRUE) # add region information
    biomassmap[is.na(region) & !duplicated(biomassmap[,.(latgrid, longrid)]), .N] # 0 missing region: good!

# Trim to just one region
presmap <- presmap[region == myreg, ]
biomassmap <- biomassmap[region == myreg, ]
fisheryspps <- fisheryspps[region == myreg, ]
	dim(presmap)
	dim(biomassmap)
	dim(fisheryspps)

# Fix a species name
presmap[spp == 'theragra chalcogramma', spp := 'gadus chalcogrammus']

############################################
## Run prioritizr just on 2006-2020
## use just one rcp
############################################
# pu.dat
# planning features are each 1/4 deg square
	pus <- presmap[,c('latgrid', 'longrid')]
	pus <- pus[!duplicated(pus),]
		dim(pus) # 2195 (ebs), 795 (goa), (bc), 229 (wc), 245 (gmex), (seus), 552 (neus), 396 (maritime), 1342 (newf)
	pus <- pus[order(pus$latgrid, pus$longrid),]
	pus$id <- 1:nrow(pus)
	pus$dummycost <- rep(cost, nrow(pus)) # set the same cost in each planning unit. can add separate costs for each zone.

	pu.dat <- pus[,c('id', 'dummycost')] # the version to write out for Marxan

	# write.csv(pus, file = paste0(folder1, '/pus.csv'))	
	# write.csv(pu.dat, file = paste0(folder1, '/pu.dat'), row.names = FALSE, quote = FALSE)

# spec.dat (takes the place of feat.dat?)
# id and name for each species
# fishery features entered separately from conservation features, even if same species
	sppstokeep <- presmap[year_range == '2007-2020' & rcp == rcp & model %in% c(1:11, 13:15, 17:18), .(poccur = mean(poccur)), by = c('latgrid', 'longrid', 'spp')] # average across models
		dim(sppstokeep)
	sppstokeep <- sppstokeep[poccur >= poccurthresh, ]
	sppstokeep <- merge(sppstokeep, pus[, .(latgrid, longrid, id)], by = c('latgrid', 'longrid')) # add pu id (and trim to focal pus)
		setnames(sppstokeep, 'id', 'pu')
		dim(sppstokeep)

		ngrid <- sppstokeep[ , .(ngrid = length(unique(pu))), by = 'spp']
		sppstokeep <- merge(sppstokeep, ngrid, by = 'spp')
		sppstokeep[ , summary(ngrid)] #

		nspps <- sppstokeep[ , .(nspp = length(unique(spp))), by = 'pu']
		sppstokeep <- merge(sppstokeep, nspps, by = 'pu')
		sppstokeep[, summary(nspp)] #

		# sppstokeep <- sppstokeep[ngrid > (nrow(pus)*0.05),] # trim to species found at poccur > poccurthresh in at least 5% of grids

		sppstokeep[ , length(unique(spp))] # 176 (ebs), 115 (goa), (bc), 115 (wc), 140 (gmex), (seus), 67 (neus), 48 (maritime), 70 (newf)

	spps <- data.table(id = 1:length(unique(sppstokeep$spp)), name = gsub(' |_', '', sort(unique(sppstokeep$spp))), spp = sort(unique(sppstokeep$spp))) #  fill spaces in species names.

	# add fishery features
	spps <- rbind(spps, data.table(id = max(spps$id) + 1:length(fisheryspps$projname), 
	                               name = paste0(gsub(' |_', '', fisheryspps$projname), '_fishery'),
	                               spp = fisheryspps$projname))
	
	# add wind and wave energy feature
	spps <- rbind(spps, data.table(id = max(spps$id) + 1, name = c('energy'), spp = c(NA)))


	spec.dat <- spps[, .(id, name)]

	# write.csv(spps, file = paste0(folder1, '/spps.csv'))
	# write.csv(spec.dat, file = paste0(folder1, '/spec.dat'), row.names = FALSE, quote = FALSE)

# puvsp.dat (takes the place of puvfeat.dat)
# which features are in each planning unit
# not documented in 1.0.1 manual. copying format from unix example
# use the average poccur and biomass from the first 16 climate models
	# Format conservation data
	puvsppa <- presmap[year_range == '2007-2020' & rcp == rcp & model %in% c(1:11, 13:15, 17:18), .(poccur = mean(poccur)), by = c('latgrid', 'longrid', 'spp')] # pres/abs data. climate models are bcc-csm1-1-m, bcc-csm1-1, CanESM2, CCSM4, CESM1-CAM5, CNRM-CM5, GFDL-CM3, GFDL-ESM2M, GFDL-ESM2G, GISS-E2-R, GISS-E2-H, HadGEM2-ES, IPSL-CM5A-LR, IPSL-CM5A-MR, MIROC-ESM, MIROC5, MPI-ESM-LR, NorESM1-ME. only use those that match biomass projections
		dim(puvsppa)
	
	puvsppa[, amount := as.numeric(poccur >= poccurthresh)] # use pres/abs as conservation amount. should this instead be left as poccur?
	    puvsppa[, summary(amount)]
	    puvsppa[, sort(unique(amount))]
	puvsppa[, poccur := NULL]

	puvsppa[ , name := gsub(' |_', '', spp)] # trim out spaces from species names
	
    # Format fishery data
    puvspbio <- biomassmap[year_range == '2007-2020' & rcp == rcp & model %in% 1:16 & spp %in% fisheryspps$projname, .(biomass = mean(biomass)), by = c('latgrid', 'longrid', 'spp')] # biomass data. models are bcc-csm1-1-m, bcc-csm1-1, CanESM2, CCSM4, CESM1-CAM5, CNRM-CM5, GFDL-CM3, GFDL-ESM2M, GFDL-ESM2G, GISS-E2-R, GISS-E2-H, IPSL-CM5A-LR, IPSL-CM5A-MR, MIROC-ESM, MPI-ESM-LR, NorESM1-ME
        dim(puvspbio)
        puvspbio[, length(unique(spp))] # should be 10

    puvspbio[, amount := biomass/100] # use biomass as amount for fishery targets. Scale down to pass presolve checks.

    puvspbio[ , name := paste0(gsub(' |_', '', spp), '_fishery')] # trim out spaces from species names
 
       
	# Format wind and wave data
	setnames(windnpv, 'npv', 'wind_npv')
	setnames(wavenpv, 'npv', 'wave_npv')
	puvenergy <- merge(windnpv, wavenpv, by = c('latgrid', 'longrid'), all = TRUE)
		head(puvenergy)
		dim(windnpv)
		dim(wavenpv)
		dim(puvenergy)
	puvenergy[wind_npv < 0 | is.na(wind_npv), wind_npv := 0] # set negative or NA NPV to 0
	puvenergy[wave_npv < 0 | is.na(wave_npv), wave_npv := 0]
	puvenergy[, amount := (wind_npv + wave_npv)/10000] # scale down to pass presolve checks
	puvenergy[, name := 'energy']
	

	# combine
	puvsp.dat <- rbind(puvsppa[, .(name, latgrid, longrid, amount, zone = 1)], 
	               puvspbio[, .(name, latgrid, longrid, amount, zone = 2)], 
	               puvenergy[, .(name, latgrid, longrid, amount, zone = 3)])
	
	# Add species ids
	nrow(puvsp.dat)
	puvsp.dat <- merge(puvsp.dat, spps[, .(id, name)], by = 'name') # merge in species IDs and trim to focal species
	nrow(puvsp.dat)
	setnames(puvsp.dat, 'id', 'species')

		
    # Add planning units
    puvsp.dat <- merge(puvsp.dat, pus[, .(latgrid, longrid, id)], by = c('latgrid', 'longrid')) # add pu id (and trim to focal pus)
	nrow(puvsp.dat)
	setnames(puvsp.dat, 'id', 'pu')
	
	# Sort and trim columns and rows
	setkey(puvsp.dat, pu, species) # order by pu then species
	puvsp.dat <- puvsp.dat[amount > 0, ] # trim only to presences
    puvsp.dat[, c('latgrid', 'longrid', 'name') := NULL]
	
    # checks
	length(unique(puvsp.dat$pu)) # planning units for species + NatCap: 2195 (ebs), 549 (neus), 1342 (newf)
	length(unique(puvsp.dat$species)) # features that are species + fishery + NatCap: 187 (ebs), 68 (neus), 71 (newf)
	puvsp.dat[, summary(amount)] # range of amounts
	
	sort(unique(table(puvsp.dat$species))) # make sure all species show up in some planning units (shouldn't see any 0s)
        
	sort(unique(table(puvsp.dat$pu))) # make sure all planning units have some species (shouldn't see any 0s)

	sort(unique(table(puvsp.dat$pu, puvsp.dat$species))) # should be all 0s and 1s

	# write.csv(puvsp.dat, file = paste0(folder1, '/puvsp.dat'), row.names = FALSE, quote = FALSE)


# zones
# id and names for each zone
	zones <- data.frame(id = 1:3, name = c('conservation', 'fishery', 'energy'))

	# write.csv(zones, file = paste0(folder1, '/zones.dat'), row.names = FALSE, quote = FALSE)


#costs
	# costs <- data.frame(costid = 1, costname = 'dummycost')
	# 
	# write.csv(costs, file = paste0(folder1, '/costs.dat'), row.names = FALSE, quote = FALSE)


#zone cost
#to adjust the importance of each cost in each zone
	# zonecost <- expand.grid(list(zoneid = zones$zoneid, costid = costs$costid))
	# zonecost$multiplier <- 0
	# zonecost$multiplier[zonecost$zoneid == 2 & zonecost$costid == 1] <- zonecosts[2] # set multiplier for conservation
	# zonecost$multiplier[zonecost$zoneid == 3 & zonecost$costid == 1] <- zonecosts[3] # set multiplier for fishery
	# zonecost$multiplier[zonecost$zoneid == 4 & zonecost$costid == 1] <- zonecosts[4] # set multiplier for energy
	# 
	# zonecost <- zonecost[zonecost$multiplier > 0,] # trim out zeros
	# 
	# write.csv(zonecost, file = paste0(folder1, '/zonecost.dat'), row.names = FALSE, quote = FALSE)

#boundary length
# optional

#zone boundary cost
# optional

#planning unit zone
#optional

#planning unit lock
#optional

#zone target
# set zone-specific targets: rows are features, columns are zones
	zonetarget <- matrix(0, nrow = nrow(spec.dat), ncol = nrow(zones), dimnames = list(spec.dat$name, zones$name))
	                     
	# set conservation zone target
	zonetarget[!grepl('energy|fishery', rownames(zonetarget)), 'conservation'] <- consgoal # XX proportion

	# set fishing zone target
	zonetarget[grepl('fishery', rownames(zonetarget)), 'fishery'] <- fishgoal # XX proportion
	
	# set energy goal target
	zonetarget[grepl('energy', rownames(zonetarget)), 'energy'] <- energygoal # XX proportion

	# write out	
	# write.csv(zonetarget.dat, file = paste0(folder1, '/zonetarget.dat'), row.names = FALSE, quote = FALSE)


#input parameters
# input <- data.frame(BLM = 0, PROP = 0.8, RANDSEED = -1, NUMREPS = 100, AVAILABLEZONE = 1, NUMITNS = '1000000', 
#                     STARTTEMP = -1, NUMTEMP = '10000', COSTTHRESH = 0, THRESHPEN1 = 0, THRESHPEN2 = 0, 
#                     INPUTDIR = './', PUNAME = 'pu.dat', SPECNAME = 'spec.dat', 
#                     PUVSPRNAME = 'puvsp.dat', ZONESNAME = 'zones.dat', COSTSNAME = 'costs.dat', 
#                     ZONECOSTNAME = 'zonecost.dat', ZONETARGETNAME = 'zonetarget.dat', SCENNAME = runname1, 
#                     SAVERUN = 3, SAVEBEST = 3, SAVESUMMARY = 3, SAVESCEN = 3, SAVETARGMET = 3, SAVESUMSOLN = 3, 
#                     SAVEPENALTY = 3, SAVELOG = 3, OUTPUTDIR = './output', RUNMODE = 1, MISSLEVEL = 1, 
#                     ITIMPTYPE = 0, HEURTYPE = -1, CLUMPTYPE = 0, VERBOSITY = 3, SAVESOLUTIONSMATRIX = 3, 
#                     SAVEANNEALINGTRACE = 0, ANNEALINGTRACEROWS = 1000)
# 
# write.table(t(input), file = paste0(folder1, '/input.dat'), row.names = TRUE, quote = FALSE, sep = ' ', col.names = FALSE)
	
	
# basic checks
colSums(zonetarget) # reasonable targets?
dim(zonetarget) # should match # features
nrow(spec.dat) # the # of features
all(rownames(zonetarget) == spec.dat$name) # target order matches species?
sum(!(puvsp.dat$pu %in% pu.dat$id)) # any planning units not in pu.dat?
sum(!(puvsp.dat$species %in% spec.dat$id)) # any species not in spec.dat?
sum(!(pu.dat$id %in% puvsp.dat$pu)) # any pus not in puvsp.dat?
sum(!(spec.dat$id %in% puvsp.dat$species)) # any species not in puvsp.dat?
	

# Define the problem in prioritzr format
# rij may need a zone column
p1 <- problem(pu.dat, spec.dat, cost_column = c('dummycost', 'dummycost', 'dummycost'), 
              rij = puvsp.dat, zones = zones) %>%
    add_min_set_objective() %>%
    add_relative_targets(zonetarget) %>%
    add_binary_decisions() %>%
    add_gurobi_solver(gap = 0)

# solve it
s1 <- solve(p1)

# force solution if fails presolve checks
s1 <- p1 %>%
    add_gurobi_solver(numeric_focus = TRUE) %>%
    solve(force = TRUE)

# check solution in various ways
s1[,.(navail = sum(solution_1_conservation == 0 & solution_1_fishery == 0 & solution_1_energy == 0), 
      ncons = sum(solution_1_conservation), 
      nfish = sum(solution_1_fishery), 
      nenergy = sum(solution_1_energy))]  # number of cells in each zone
r1 <- feature_representation(p1, s1[, .(solution_1_conservation, solution_1_fishery, solution_1_energy)]) # representation of each feature in each zone
r1dt <- as.data.table(r1)
r1dt[zone=='conservation',] # meeting conservation targets?
r1dt[zone=='conservation' & !grepl('fishery|energy', feature), .(summary(absolute_held), summary(relative_held))] # meeting conservation targets?
r1dt[zone=='fishery' & grepl('fishery', feature),] # fishery targets
r1dt[zone=='energy',]


# write out
write.csv(s1, file = paste0(marxfolder, 'solution_', runname1, '.csv'))

##################################################################
## Set up a Marxan run on 2006-2020 and ensemble mean 2081-2100
## This assumes that the historical-only code (previous section) has been run and is loaded in memory
##################################################################
# Create directory for input and output if missing
if(!dir.exists(inputfolder2)){
	dir.create(inputfolder2)
}
if(!dir.exists(outputfolder2)){
	dir.create(outputfolder2)
}


# pu.dat
# planning features are each 1/4 deg square
	pus2 <- presmap[,c('lat', 'lon')]
	pus2 <- pus2[!duplicated(pus2),]
		dim(pus2)
	pus2 <- pus2[order(pus2$lat, pus2$lon),]
	pus2$id <- 1:nrow(pus2)
	pus2$dummycost <- rep(cost, nrow(pus2))  # set the same cost in each planning unit. can add separate costs for each zone.

	pu2.dat<-pus2[,c('id', 'dummycost')]

	save(pus2, file=paste(inputfolder2, '/pus.Rdata', sep=''))	
	write.csv(pu2.dat, file=paste(inputfolder2, '/pu.dat', sep=''), row.names=FALSE, quote=FALSE)

# spec.dat (takes the place of feat.dat?)
# document every species present
# not documented in 1.0.1 manual. copying format from unix example
	spps2 <- spps #  use the same species as in the historical-only run
		spps2$name <- as.character(spps2$name)

	sppinds <- !grepl('energy', spps$name) # don't include energy in each time period
	temp1 <- spps2[sppinds,]
	spps2$name[sppinds] <- paste(spps2$name[sppinds], gsub('-', '', planningperiods[1]), sep='')
	temp1$name <- paste(temp1$name, gsub('-', '', planningperiods[2]), sep='')
	temp1$id = temp1$id + max(spps2$id) # make sure the ids don't overlap
	spps2 <- rbind(spps2, temp1)

	# set feature penalty factor
	spps2$fpf <- fpf


	spps2.dat <- spps2[,c('id', 'fpf', 'name')]

	save(spps2, file=paste(inputfolder2, '/spps.Rdata', sep=''))		
	write.csv(spps2.dat, file=paste(inputfolder2, '/spec.dat', sep=''), row.names=FALSE, quote=FALSE)

# puvsp.dat (takes the place of puvfeat.dat?)
# table of species by planning units
# not documented in 1.0.1 manual. copying format from unix example
	puvsp2 <- presmap[presmap$period %in% c('2006-2020', '2081-2100'),c('lat', 'lon', 'period', 'sppocean', 'wtcpue.proj', 'pres')]
		dim(puvsp2)
	puvsp2 <- merge(puvsp2, pus2[,c('lat', 'lon', 'id')]) # add pu id (and trim to focal pus2)
		dim(puvsp2)
		names(puvsp2)[names(puvsp2)=='id'] <- 'pu'

	puvsp2$name <- gsub(' |_|-', '', paste(puvsp2$sppocean, puvsp2$period)) # create names, and remove spaces and dashes
		dim(puvsp2)
	puvsp2 <- merge(puvsp2, spps2[,c('id', 'name')]) # add species id and trim to focal species
		dim(puvsp2)
		names(puvsp2)[names(puvsp2)=='id'] <- 'species'
	puvsp2$amount <- as.numeric(puvsp2[[amountcolnm]]) # use projected biomass or pres as amount
	puvsp2$amount[!puvsp2$pres] <- 0 # set amount to zero where our cutoff says not present

	# Combine species, wind, and wave data for output
	puvsp2.dat <- rbind(puvsp2[,c('species', 'pu', 'amount')], puvenergy[,c('species', 'pu', 'amount')])
	puvsp2.dat <- puvsp2.dat[order(puvsp2.dat$pu, puvsp2.dat$species),]
	puvsp2.dat <- puvsp2.dat[puvsp2.dat$amount>0,] # trim only to presences

	length(unique(puvsp2$pu)) # 552 (neus), 677 (ebs), 
	length(unique(puvsp2$species)) # 134 (neus), 252 (ebs)
	length(unique(puvsp2.dat$pu)) # 552 (neus), 677 (ebs)
	length(unique(puvsp2.dat$species)) # 129 (neus), 235 (ebs)
	range(table(puvsp2.dat$species)) # make sure all species show up in some planning units
	range(table(puvsp2.dat$pu)) # make sure all planning units have some species

	write.csv(puvsp2.dat, file=paste(inputfolder2, '/puvsp.dat', sep=''), row.names=FALSE, quote=FALSE)


# zones
# id and names for each zone
	zones2 <- data.frame(zoneid=1:4, zonename=c('available', 'conservation', 'fishery', 'energy'))
	write.csv(zones2, file=paste(inputfolder2, '/zones.dat', sep=''), row.names=FALSE, quote=FALSE)

#costs
	costs2 <- data.frame(costid=1, costname='dummycost')	
	write.csv(costs2, file=paste(inputfolder2, '/costs.dat', sep=''), row.names=FALSE, quote=FALSE)


#zone cost
	zonecost2 <- expand.grid(list(zoneid=zones2$zoneid, costid=costs2$costid))
	zonecost2$multiplier <- 0
	zonecost2$multiplier[zonecost2$zoneid==1 & zonecost2$costid==1] <- zonecosts[1] # set multiplier for available
	zonecost2$multiplier[zonecost2$zoneid==2 & zonecost2$costid==1] <- zonecosts[2] # set multiplier for conservation
	zonecost2$multiplier[zonecost2$zoneid==3 & zonecost2$costid==1] <- zonecosts[3] # set multiplier for fishery
	zonecost2$multiplier[zonecost2$zoneid==4 & zonecost2$costid==1] <- zonecosts[4] # set multiplier for energy

	zonecost2 <- zonecost2[zonecost2$multiplier>0,]

	write.csv(zonecost2, file=paste(inputfolder2, '/zonecost.dat', sep=''), row.names=FALSE, quote=FALSE)

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
	zonetarget2 <- expand.grid(list(zoneid=zones2$zoneid, speciesid=spps2$id))
	zonetarget2$target <- 0

	# set conservation zone target
	consinds <- zonetarget2$zoneid==zones2$zoneid[zones2$zonename=='conservation'] & zonetarget2$speciesid %in% spps2$id[spps$name != 'energy']
	zonetarget2$target[consinds] <- consgoal # XX proportion
	if(amountcolnm=='pres'){
		zonetarget2$targettype[consinds] <- 1 # 3: proportion of total occurrences. 1: proportion of total amount
	}
	if(amountcolnm=='wtcpue.proj'){
		zonetarget2$targettype[consinds] <- 3 # 3: proportion of total occurrences. 1: proportion of total amount
	}

	# set fishing zone target
	fishinds <- zonetarget2$zoneid==zones2$zoneid[zones2$zonename=='fishery'] & zonetarget2$speciesid %in% spps2$id[spps2$sppocean %in% fisheryspps$projname]
	zonetarget2$target[fishinds] <- fishgoal # XX proportion
	zonetarget2$targettype[fishinds] <- 1 # 3: proportion of total occurrences. 1: proportion of total amount

	# set energy goal target
	energyinds <- zonetarget2$zoneid==zones2$zoneid[zones2$zonename=='energy'] & zonetarget2$speciesid %in% spps2$id[spps2$name == 'energy']
	zonetarget2$target[energyinds] <- energygoal # XX proportion
	zonetarget2$targettype[energyinds] <- 1 # 3: proportion of total occurrences. 1: proportion of total amount

	# format for output	
	zonetarget2.dat <- zonetarget2[zonetarget2$target>0,] # trim to only positive targets
	zonetarget2.dat <- zonetarget2.dat[order(zonetarget2.dat$zoneid, zonetarget2.dat$speciesid),]

	# write out	
	write.csv(zonetarget2.dat, file=paste(inputfolder2, '/zonetarget.dat', sep=''), row.names=FALSE, quote=FALSE)



#input parameters
input2 <- data.frame(BLM=0, PROP=0.5, RANDSEED=-1, NUMREPS=100, AVAILABLEZONE=1, NUMITNS='1000000', STARTTEMP=-1, NUMTEMP='10000', COSTTHRESH=0, THRESHPEN1=0, THRESHPEN2=0, INPUTDIR=paste(runname2, '_input', sep=''), PUNAME='pu.dat', SPECNAME='spec.dat', PUVSPRNAME='puvsp.dat', ZONESNAME='zones.dat', COSTSNAME='costs.dat', ZONECOSTNAME='zonecost.dat', ZONETARGETNAME='zonetarget.dat', SCENNAME=runname2, SAVERUN=3, SAVEBEST=3, SAVESUMMARY=3, SAVESCEN=3, SAVETARGMET=3, SAVESUMSOLN=3, SAVEPENALTY=3, SAVELOG=3, OUTPUTDIR=paste(runname2, '_output', sep=''), RUNMODE=1, MISSLEVEL=1, ITIMPTYPE=0, HEURTYPE=-1, CLUMPTYPE=0, VERBOSITY=3, SAVESOLUTIONSMATRIX=3, SAVEANNEALINGTRACE=0, ANNEALINGTRACEROWS=1000)


write.table(t(input2), file=paste(marxfolder, 'input.dat', sep=''), row.names=TRUE, quote=FALSE, sep=' ', col.names=FALSE)
write.table(t(input2), file=paste(inputfolder2, '/input.dat', sep=''), row.names=TRUE, quote=FALSE, sep=' ', col.names=FALSE)
	

